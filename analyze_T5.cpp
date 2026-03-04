#include <ostream>
#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>

#include <ROOT/RVec.hxx>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TF2.h>
#include <TLatex.h>

#include <nlohmann/json.hpp>

#include "./utils.h"
#include "./return_TOF_position.h"
#include "TString.h"
#include "buffer.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::cerr;

using namespace ROOT;

int main(int argc, char** argv){
	
	if (argc < 2){
		cerr << "One argument expected, but was not provided. Usage: ./analyze_T5 <run_number>" << endl;
		return -1;
	}
//	TApplication app("app", &argc, argv);

	string arg = argv[1];
	auto run_number = std::stoi(arg);
	RUN_NUMBER = run_number;
	TString filename = "WCTE_data/charged_particle/WCTE_offline_R" + arg + "S0_VME_matched.root";
	auto file = TFile::Open(filename, "READ");
	
	if (!file || file-> IsZombie()){
		cerr << "ERROR, file did not open" << endl;
		return -1;
	}

	std::ifstream config_file("/home/frantisek/scripts/config.json");
	nlohmann::json config = nlohmann::json:: parse(config_file);
	auto& run_config = config[std::to_string(run_number)];
	BEAM_MOMENTUM = run_config.value("Beam momentum (MeV/c)", 0);
	cout << "Config file loaded, beam momentum is " << BEAM_MOMENTUM << "MeV/c" << endl ;


	auto tree = file->Get<TTree>("WCTEReadoutWindows");


	vector<float>* arr_bm_times = nullptr;
	vector<float>* arr_bm_charges = nullptr;
	vector<int>* arr_bm_time_ids = nullptr;
	vector<int>* arr_bm_charge_ids = nullptr;

	vector<double>* arr_pmt_times = nullptr;
	vector<int>* arr_pmt_ids = nullptr;
	vector<int>* arr_mpmt_ids = nullptr;

	tree-> SetBranchStatus("*", 0);
	tree-> SetBranchStatus("beamline_pmt_qdc_ids", 1);
	tree-> SetBranchStatus("beamline_pmt_tdc_ids", 1);
	tree-> SetBranchStatus("beamline_pmt_tdc_times", 1);
	tree-> SetBranchStatus("beamline_pmt_qdc_charges", 1);
	tree-> SetBranchStatus("hit_pmt_times", 1);
	tree-> SetBranchStatus("hit_mpmt_card_ids", 1);
	tree-> SetBranchStatus("hit_pmt_channel_ids", 1);

	tree-> SetBranchAddress("beamline_pmt_qdc_ids", &arr_bm_charge_ids);
	tree-> SetBranchAddress("beamline_pmt_tdc_ids", &arr_bm_time_ids);
	tree-> SetBranchAddress("beamline_pmt_tdc_times", &arr_bm_times);
	tree-> SetBranchAddress("beamline_pmt_qdc_charges", &arr_bm_charges);
	tree-> SetBranchAddress("hit_pmt_times", &arr_pmt_times);
	tree-> SetBranchAddress("hit_mpmt_card_ids", &arr_mpmt_ids);
	tree-> SetBranchAddress("hit_pmt_channel_ids", &arr_pmt_ids);


	Cuts cut;
	TOF_reconstructor recon;
	Histograms hists;
	setup_histograms(hists, recon);
	int n_pass_cut = 0;
	int n_T5_valid_events = 0;
	int n_events = tree->GetEntries();
	int verb = 1000;
	int n_events_with_multiple_valid_hits = 0;
	int n_events_with_valid_hits_in_expected_window = 0;
	int n_events_with_multiple_scint_hits = 0;
	int n_events_with_multiple_valid_hits_had_one_in_expected_window = 0;
	int n_invalid_hits = 0;
	int n_events_out_of_bounds = 0;

	vector<event_T5_detection> invalid_T5_hits;
	vector<event_T5_detection> multivalidhits_events;
	

	for(size_t i = 0; i < n_events; i++){
		tree->GetEntry(i);
		//Print progress

		if (i % verb == 0) cout << "\rAnalyzed " << i << " of " << n_events << std::setprecision(2) << std::fixed <<
			" events (" << static_cast<float>(i)/n_events * 100 << " %)" << std::flush << endl; 

		RVecI bm_time_ids(arr_bm_time_ids->data(), arr_bm_time_ids->size());
		RVecI bm_charge_ids(arr_bm_charge_ids->data(), arr_bm_charge_ids->size());
		RVecF bm_times(arr_bm_times->data(), arr_bm_times->size());
		RVecF bm_charges(arr_bm_charges->data(), arr_bm_charges->size());

		RVecD pmt_times(arr_pmt_times->data(), arr_pmt_times->size());
		RVecI pmt_ids(arr_pmt_ids->data(), arr_pmt_ids->size());
		RVecI mpmt_ids(arr_mpmt_ids->data(), arr_mpmt_ids->size());

		
		if (!cut.hit_T0_T1(bm_time_ids, bm_charge_ids) ||
		    !cut.hit_T4(bm_charge_ids, bm_charges) ||
		    !cut.did_not_hit_HC(bm_charge_ids, bm_charges) ||
		    !cut.hit_T5(mpmt_ids, pmt_ids))
			continue;

		n_pass_cut++;

		auto mask_T5_board = (mpmt_ids == cut.get_T5_board());
		auto T5_board_ids = pmt_ids[mask_T5_board];
		auto T5_board_times = pmt_times[mask_T5_board];

		auto detections = recon.Return_position(i, *arr_mpmt_ids, *arr_pmt_ids, *arr_pmt_times);

		if (detections.HasValidHit) n_T5_valid_events++;
		if (detections.HasMultipleValidHits){
			n_events_with_multiple_valid_hits++;
			multivalidhits_events.push_back(detections);
			if (detections.HasInTimeWindow) n_events_with_multiple_valid_hits_had_one_in_expected_window++;
			if (detections.HasMultipleScintillatorsHit) n_events_with_multiple_scint_hits++;
		}
		if (detections.HasInTimeWindow) n_events_with_valid_hits_in_expected_window++;
		if (detections.HasHit && !detections.HasValidHit) {invalid_T5_hits.push_back(detections); n_invalid_hits++;}
		if (detections.HasOutOfBounds) n_events_out_of_bounds++;


		int n_hits_in_T5_in_single_event = 0;
		for (int i = 0; i < cut.Get_T5_ids().size(); i++){
			auto T5_id = cut.Get_T5_ids().at(i);
			int sum_hits_T5_i = VecOps::Sum(T5_board_ids == T5_id);
			hists.fill(Form("T5_number_of_hits_%i", i), sum_hits_T5_i);
			n_hits_in_T5_in_single_event += sum_hits_T5_i;
		}
		hists.fill("n_event_hits", n_hits_in_T5_in_single_event);	
		
		for (const auto& hit : detections.T5_hits){
			if (!hit.is_valid_hit || hit.quality != HitQuality::Perfect){
				continue;
			}
			hists.fill("positions", hit.position->first, hit.position->second);
		}

	}
	auto hist = hists.get_histogram_2D("positions");
	TF2* gaus_2D = new TF2("gaus_2D", "bigaus", recon.Get_scint_xmin(3), recon.Get_scint_xmax(3), recon.Get_ymin(), recon.Get_ymax());
	gaus_2D->SetParameters(130, 0, 40, 0, 40, 0);
	hist->Fit(gaus_2D, "RQ");

	gaus_2D = (TF2*)hist->GetFunction("gaus_2D");

	double volume = gaus_2D->GetParameter(0);
	double sig_x  = gaus_2D->GetParameter(2);
	double sig_y  = gaus_2D->GetParameter(4);
	double rho    = gaus_2D->GetParameter(5); // Correlation factor

	// 3. Calculate the TRUE mathematical peak height of the bigaus function
	double denominator = 2.0 * TMath::Pi() * sig_x * sig_y * std::sqrt(1.0 - rho*rho);
	double peak_amplitude = volume / denominator;

	// 4. Define your contour levels! 
	// 1-sigma drops to e^(-0.5)
	// 2-sigma drops to e^(-2.0)
	// 3-sigma drops to e^(-4.5)

	// Let's draw all 3 levels to make it look incredibly professional:
	double contours[3];
	contours[0] = peak_amplitude * std::exp(-4.5); // 3-sigma (widest, lowest)
	contours[1] = peak_amplitude * std::exp(-2.0); // 2-sigma
	contours[2] = peak_amplitude * std::exp(-0.5); // 1-sigma (tightest, highest)

	// 5. Apply the contours to your TF2
	// The arguments are: (number_of_levels, array_of_levels)
	gaus_2D->SetContour(3, contours);

	// 6. Make the contour lines stand out against the color map
	gaus_2D->SetLineColor(kRed);
	gaus_2D->SetLineWidth(2);
	gaus_2D->SetLineStyle(1); // Solid lines

	double sigma_x = gaus_2D->GetParameter(2);
	double sigma_y = gaus_2D->GetParameter(4);

	TString txt_sigma_x = Form("#sigma_{x} = %.2f mm", sigma_x);
	TString txt_sigma_y = Form("#sigma_{y} = %.2f mm", sigma_y);

	double offset = 10.0;
	TLatex* ltx_sigX = new TLatex(sigma_x - offset, sigma_y, txt_sigma_x);
	TLatex* ltx_sigY = new TLatex(sigma_x - offset, sigma_y - 5.0, txt_sigma_y);

	ltx_sigX->SetTextSize(0.04); ltx_sigX->SetTextColor(kBlack);
	ltx_sigY->SetTextSize(0.04); ltx_sigY->SetTextColor(kBlack);

	TLatex* contour_sigma = new TLatex(sigma_x - offset, -sigma_y + offset, "1#sigma");
	contour_sigma->SetTextSize(0.04); contour_sigma->SetTextColor(kBlack);

	hist->GetListOfFunctions()->Add(ltx_sigX);
	hist->GetListOfFunctions()->Add(ltx_sigY);
	hist->GetListOfFunctions()->Add(contour_sigma);


	for (int i = 0; i < 8; i++){
		TString h_name = "positions_" + std::to_string(i); 
		hists. hist_projectX("positions", h_name.Data(), i+1, i+1);
		hists. get_histogram(h_name.Data())->Fit("gaus", "QR", "", recon.Get_scint_xmin(i), recon. Get_scint_xmax(i));
	}
	TString plots_directory = "plots/Run_" + std::to_string(run_number);
	gSystem->Exec("mkdir -p " + plots_directory);
	gSystem->cd(plots_directory);
	hists.print_exclusive("positions", 1000, 900);
	hists.print_all();
	hists.save_all("hists");

	cout << "Printing out invalid hits: " << endl;
	for (const auto& event : invalid_T5_hits){
		cout << "Event " << event.event_nr << ":" << endl;
		for(const auto& hit : event.T5_hits){
			cout << "Scintillator ID: " << hit.scintillator_id.value() << "\t"
			     << "Average time: " << hit.hit_time.value() << "\t"
			     << "First SiPM time: " << hit.sipm_time_a.value() << "\t"
			     << "Second SiPM time: " << hit.sipm_time_b.value() << "\t"
			     << endl;
		}	
	}

	cout << endl;
	cout << "Printing out events with multiple valid hits: " << endl;
	for (const auto& event : multivalidhits_events){
		cout << "Event " << event.event_nr << ": " << endl;
		for (const auto& hit : event.T5_hits){
			cout << "Scintillator ID: " << hit.scintillator_id.value() << "\t"
			     << "Average time: " << hit.hit_time.value() << "\t"
			     << "First SiPM time: " << hit.sipm_time_a.value() << "\t"
			     << "Second SiPM time: " << hit.sipm_time_b.value() << "\t"
			     << "Valid hit?: " << hit.is_valid_hit << "\t"
			     << endl;
		}
	}

	cout << endl << n_pass_cut << " events out of " << n_events << " passed cuts" << endl;
	cout << n_T5_valid_events << " events got a valid reconstruction -- "
	     << n_pass_cut - n_T5_valid_events << " were mismatched events?" << endl
	     << n_invalid_hits << " events were invalid" << endl
	     << n_events_with_valid_hits_in_expected_window << " events of them had a hit in the expected time window" << endl
	     << n_events_out_of_bounds << " events had a reconstruction out of bounds" << endl
	     << endl
		

	     << n_events_with_multiple_valid_hits << " events had multiple valid hits -- " 
	     << n_events_with_multiple_valid_hits_had_one_in_expected_window << " of those had at least one hit in the expected time window" << endl
	     << n_events_with_multiple_scint_hits << " events had hits in multiple scintillators -- in the expected time window" << endl
	     << endl; 
	cout << endl;

	std::ofstream file_out;
	file_out.open("Beam_profile_widths.dat", std::ofstream::app);
	file_out << run_number << "\t"
		<< BEAM_MOMENTUM << "\t"
		<< sigma_x << "\t"
		<< sigma_y << "\t" << endl;


	//	app.Run();

	return 0;
}
