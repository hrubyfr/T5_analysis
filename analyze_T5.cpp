#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

#include <ROOT/RVec.hxx>
#include <TApplication.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// #include <nlohmann/json.hpp>

#include "return_TOF_position.h"
#include "utils.h"
// #include "buffer.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace ROOT;

int main(int argc, char **argv) {

  int run_number = -1;
  vector<string> input_paths;
  TString output_path = "output.root";
  bool debug = false;

  int opt;
  while ((opt = getopt(argc, argv, "r:i:o:d")) != -1) {
    switch (opt) {
    case 'r':
      run_number = std::stoi(optarg);
      break;
    case 'i':
      input_paths.push_back(optarg);
      break;
    case 'o':
      output_path = optarg;
      break;
    case 'd':
      debug = true;
      break;
    default:
      cerr << "Usage: " << argv[0]
           << " -r <run_number> [-i <input_file>] [-o <output_file>] [-d]"
           << endl;
      return -1;
    }
  }
  // debug print the input paths
  cout << "Input paths: " << endl;
  for (size_t input_idx = 0; input_idx < input_paths.size(); ++input_idx) {
    cout << input_paths[input_idx] << endl;
  }

  RUN_NUMBER = run_number;

  if (input_paths.empty()) {
    throw std::runtime_error("No input files provided");
  }

  for (size_t input_idx = 0; input_idx < input_paths.size(); ++input_idx) {
    TString filename = input_paths[input_idx];
    TString current_output_path = output_path;

    TString base_name = gSystem->BaseName(filename);
    base_name.ReplaceAll(".root", "_T5.root");

    current_output_path = output_path + "/" + base_name;

    auto file = TFile::Open(filename, "READ");

    if (!file || file->IsZombie()) {
      cerr << "ERROR, file did not open: " << filename << endl;
      continue;
    }

    // std::ifstream config_file("/home/frantisek/scripts/config.json");
    // nlohmann::json config = nlohmann::json:: parse(config_file);
    // auto& run_config = config[std::to_string(run_number)];
    // BEAM_MOMENTUM = run_config.value("Beam momentum (MeV/c)", 0);
    // cout << "Config file loaded, beam momentum is " << BEAM_MOMENTUM <<
    // "MeV/c" << endl ;

    bool hardware_processed_data = false;
    string tree_name = "";
    if (file->GetListOfKeys()->Contains("ProcessedWaveforms")) {
      tree_name = "ProcessedWaveforms";
      hardware_processed_data = true;
    } else if (file->GetListOfKeys()->Contains("WCTEReadoutWindows")) {
      tree_name = "WCTEReadoutWindows";
    } else {
      cerr << "ERROR: Input file does not contain recognized tree "
              "(WCTEReadoutWindows or ProcessedWaveforms)"
           << endl;
      return -1;
    }
    TTreeReader tree(tree_name.c_str(), file);

    // vector<float>* arr_bm_times = nullptr;
    // vector<float>* arr_bm_charges = nullptr;
    // vector<int>* arr_bm_time_ids = nullptr;
    // vector<int>* arr_bm_charge_ids = nullptr;

    std::unique_ptr<TTreeReaderValue<vector<double>>> arr_pmt_times;
    std::unique_ptr<TTreeReaderValue<vector<int>>> arr_pmt_ids;
    std::unique_ptr<TTreeReaderValue<vector<int>>> arr_mpmt_ids;

    // needed for array input of the processed waveforms file

    const int MAX_HITS = 4000;

    std::unique_ptr<TTreeReaderValue<int>> n_pmt_times;
    std::unique_ptr<TTreeReaderValue<int>> n_pmt_chans;
    std::unique_ptr<TTreeReaderValue<int>> n_pmt_cards;

    std::unique_ptr<TTreeReaderArray<double>> processed_hit_times;
    std::unique_ptr<TTreeReaderArray<int>> processed_hit_cards;
    std::unique_ptr<TTreeReaderArray<int>> processed_hit_chans;

    if (hardware_processed_data) {
      n_pmt_times = std::make_unique<TTreeReaderValue<int>>(tree, "nhit_time");
      n_pmt_chans = std::make_unique<TTreeReaderValue<int>>(tree, "nhit_chan");
      n_pmt_cards = std::make_unique<TTreeReaderValue<int>>(tree, "nhit_card");

      processed_hit_times =
          std::make_unique<TTreeReaderArray<double>>(tree, "hit_time");
      processed_hit_cards =
          std::make_unique<TTreeReaderArray<int>>(tree, "hit_card");
      processed_hit_chans =
          std::make_unique<TTreeReaderArray<int>>(tree, "hit_chan");

    } else {
      arr_pmt_times = std::make_unique<TTreeReaderValue<vector<double>>>(
          tree, "hit_pmt_times");
      arr_pmt_ids = std::make_unique<TTreeReaderValue<vector<int>>>(
          tree, "hit_pmt_channel_ids");
      arr_mpmt_ids = std::make_unique<TTreeReaderValue<vector<int>>>(
          tree, "hit_mpmt_card_ids");
    }

    // tree-> SetBranchAddress("beamline_pmt_qdc_ids", &arr_bm_charge_ids);
    // tree-> SetBranchAddress("beamline_pmt_tdc_ids", &arr_bm_time_ids);
    // tree-> SetBranchAddress("beamline_pmt_tdc_times", &arr_bm_times);
    // tree-> SetBranchAddress("beamline_pmt_qdc_charges", &arr_bm_charges);

    Cuts cut;
    TOF_reconstructor recon;
    // Histograms hists;
    // setup_histograms(hists, recon);
    int n_pass_cut = 0;
    int n_T5_valid_events = 0;
    auto n_events = tree.GetEntries();
    if (debug) {
      cout << "Debug mode enabled: limiting to 5000 events" << endl;
      n_events = std::min(n_events, 5000LL);
    }
    int verb = 1000;
    int n_events_with_multiple_valid_hits = 0;
    int n_events_with_valid_hits_in_expected_window = 0;
    int n_events_with_multiple_scint_hits = 0;
    int n_events_with_multiple_valid_hits_had_one_in_expected_window = 0;
    int n_invalid_hits = 0;
    int n_events_out_of_bounds = 0;

    vector<event_T5_detection> all_T5_hits;

    while (tree.Next()) {
      auto i = tree.GetCurrentEntry();

      vector<double> vec_hit_time;
      vector<int> vec_hit_chan;
      vector<int> vec_hit_card;
      // Print progress
      event_T5_detection detections;
      if (hardware_processed_data) {
        if (**n_pmt_times >= MAX_HITS) {
          // continue and push back dummy detections
          std::cerr << "ERROR: Too many hits in event "
                    << tree.GetCurrentEntry() << " nhits " << **n_pmt_times
                    << " (max " << MAX_HITS << ")" << std::endl;
          detections.event_nr = tree.GetCurrentEntry();
          all_T5_hits.push_back(detections);
          continue;
        }
        // change array output into std vectors
        vec_hit_time.assign((*processed_hit_times).begin(),
                            (*processed_hit_times).end());
        vec_hit_card.assign((*processed_hit_cards).begin(),
                            (*processed_hit_cards).end());
        vec_hit_chan.assign((*processed_hit_chans).begin(),
                            (*processed_hit_chans).end());
      } else {
        vec_hit_time.assign((**arr_pmt_times).begin(), (**arr_pmt_times).end());
        vec_hit_card.assign((**arr_mpmt_ids).begin(), (**arr_mpmt_ids).end());
        vec_hit_chan.assign((**arr_pmt_ids).begin(), (**arr_pmt_ids).end());
      }

      if (i % verb == 0)
        cout << "\rAnalyzed " << i << " of " << n_events << std::setprecision(2)
             << std::fixed << " events ("
             << static_cast<float>(i) / n_events * 100 << " %)" << std::flush
             << endl;

      // RVecI bm_time_ids(arr_bm_time_ids->data(), arr_bm_time_ids->size());
      // RVecI bm_charge_ids(arr_bm_charge_ids->data(),
      // arr_bm_charge_ids->size()); RVecF bm_times(arr_bm_times->data(),
      // arr_bm_times->size()); RVecF bm_charges(arr_bm_charges->data(),
      // arr_bm_charges->size());

      RVecD pmt_times(vec_hit_time.data(), vec_hit_time.size());
      RVecI pmt_ids(vec_hit_chan.data(), vec_hit_chan.size());
      RVecI mpmt_ids(vec_hit_card.data(), vec_hit_card.size());

      if (!cut.hit_T5(mpmt_ids, pmt_ids)) {
        detections.event_nr = i;
        all_T5_hits.push_back(detections);
        continue;
      }

      n_pass_cut++;

      auto mask_T5_board = (mpmt_ids == cut.get_T5_board());
      auto T5_board_ids = pmt_ids[mask_T5_board];
      auto T5_board_times = pmt_times[mask_T5_board];

      detections =
          recon.Return_position(i, vec_hit_card, vec_hit_chan, vec_hit_time);

      if (detections.HasValidHit)
        n_T5_valid_events++;
      if (detections.HasMultipleValidHits) {
        n_events_with_multiple_valid_hits++;
        if (detections.HasInTimeWindow)
          n_events_with_multiple_valid_hits_had_one_in_expected_window++;
        if (detections.HasMultipleScintillatorsHit)
          n_events_with_multiple_scint_hits++;
      }
      if (detections.HasInTimeWindow)
        n_events_with_valid_hits_in_expected_window++;
      if (detections.HasHit && !detections.HasValidHit) {
        n_invalid_hits++;
      }
      if (detections.HasOutOfBounds)
        n_events_out_of_bounds++;

      int n_hits_in_T5_in_single_event = 0;
      for (size_t j = 0; j < cut.Get_T5_ids().size(); j++) {
        auto T5_id = cut.Get_T5_ids().at(j);
        int sum_hits_T5_i = VecOps::Sum(T5_board_ids == T5_id);
        // hists.fill(Form("T5_number_of_hits_%i", j), sum_hits_T5_i);
        n_hits_in_T5_in_single_event += sum_hits_T5_i;
      }
      // hists.fill("n_event_hits", n_hits_in_T5_in_single_event);

      // for (const auto& hit : detections.T5_hits){
      // 	if (!hit.is_valid_hit || hit.quality != HitQuality::Perfect){
      // 		continue;
      // 	}
      // 	// hists.fill("positions", hit.position_x, hit.position_y);
      // }
      all_T5_hits.push_back(detections);
    }
    // auto hist = hists.get_histogram_2D("positions");
    // TF2* gaus_2D = new TF2("gaus_2D", "bigaus", recon.Get_scint_xmin(3),
    // recon.Get_scint_xmax(3), recon.Get_ymin(), recon.Get_ymax());
    // gaus_2D->SetParameters(130, 0, 40, 0, 40, 0);
    // hist->Fit(gaus_2D, "R");
    //
    // gaus_2D = (TF2*)hist->GetFunction("gaus_2D");
    //
    // double volume = gaus_2D->GetParameter(0);
    // double sig_x  = gaus_2D->GetParameter(2);
    // double sig_y  = gaus_2D->GetParameter(4);
    // double rho    = gaus_2D->GetParameter(5); // Correlation factor
    //
    // // 3. Calculate the TRUE mathematical peak height of the bigaus function
    // double denominator = 2.0 * TMath::Pi() * sig_x * sig_y * std::sqrt(1.0 -
    // rho*rho); double peak_amplitude = volume / denominator;
    //
    // // 4. Define your contour levels!
    // // 1-sigma drops to e^(-0.5)
    // // 2-sigma drops to e^(-2.0)
    // // 3-sigma drops to e^(-4.5)
    //
    // // Let's draw all 3 levels to make it look incredibly professional:
    // double contours[3];
    // contours[0] = peak_amplitude * std::exp(-4.5); // 3-sigma (widest,
    // lowest) contours[1] = peak_amplitude * std::exp(-2.0); // 2-sigma
    // contours[2] = peak_amplitude * std::exp(-0.5); // 1-sigma (tightest,
    // highest)
    //
    // // 5. Apply the contours to your TF2
    // // The arguments are: (number_of_levels, array_of_levels)
    // gaus_2D->SetContour(3, contours);
    //
    // // 6. Make the contour lines stand out against the color map
    // gaus_2D->SetLineColor(kRed);
    // gaus_2D->SetLineWidth(2);
    // gaus_2D->SetLineStyle(1); // Solid lines
    //
    // double sigma_x = gaus_2D->GetParameter(2);
    // double sigma_y = gaus_2D->GetParameter(4);
    //
    // TString txt_sigma_x = Form("#sigma_{x} = %.2f mm", sigma_x);
    // TString txt_sigma_y = Form("#sigma_{y} = %.2f mm", sigma_y);
    //
    // double offset = 10.0;
    // TLatex* ltx_sigX = new TLatex(sigma_x - offset, sigma_y, txt_sigma_x);
    // TLatex* ltx_sigY = new TLatex(sigma_x - offset, sigma_y - 5.0,
    // txt_sigma_y);
    //
    // ltx_sigX->SetTextSize(0.04); ltx_sigX->SetTextColor(kBlack);
    // ltx_sigY->SetTextSize(0.04); ltx_sigY->SetTextColor(kBlack);
    //
    // TLatex* contour_sigma = new TLatex(sigma_x - offset, -sigma_y + offset,
    // "1#sigma"); contour_sigma->SetTextSize(0.04);
    // contour_sigma->SetTextColor(kBlack);
    //
    // hist->GetListOfFunctions()->Add(ltx_sigX);
    // hist->GetListOfFunctions()->Add(ltx_sigY);
    // hist->GetListOfFunctions()->Add(contour_sigma);

    // for (int i = 0; i < 8; i++){
    // 	TString h_name = "positions_" + std::to_string(i);
    // 	hists. hist_projectX("positions", h_name.Data(), i+1, i+1);
    // 	hists. get_histogram(h_name.Data())->Fit("gaus", "QR", "",
    // recon.Get_scint_xmin(i), recon. Get_scint_xmax(i));
    // }
    // TString plots_directory = "plots/Run_" + std::to_string(run_number);
    // gSystem->Exec("mkdir -p " + plots_directory);
    // gSystem->cd(plots_directory);
    // hists.print_exclusive("positions", 1000, 900);
    // hists.print_all();
    // hists.save_all("hists");

    // cout << "Printing out invalid hits: " << endl;
    // for (const auto& event : invalid_T5_hits){
    // 	cout << "Event " << event.event_nr << ":" << endl;
    // 	for(const auto& hit : event.T5_hits){
    // 		cout << "Scintillator ID: " << hit.scintillator_id.value() <<
    // "\t"
    // 		     << "Average time: " << hit.hit_time.value() << "\t"
    // 		     << "First SiPM time: " << hit.sipm_time_a.value() << "\t"
    // 		     << "Second SiPM time: " << hit.sipm_time_b.value() << "\t"
    // 		     << endl;
    // 	}
    // }
    //
    // cout << endl;
    // cout << "Printing out events with multiple valid hits: " << endl;
    // for (const auto& event : multivalidhits_events){
    // 	cout << "Event " << event.event_nr << ": " << endl;
    // 	for (const auto& hit : event.T5_hits){
    // 		cout << "Scintillator ID: " << hit.scintillator_id.value() <<
    // "\t"
    // 		     << "Average time: " << hit.hit_time.value() << "\t"
    // 		     << "First SiPM time: " << hit.sipm_time_a.value() << "\t"
    // 		     << "Second SiPM time: " << hit.sipm_time_b.value() << "\t"
    // 		     << "Valid hit?: " << hit.is_valid_hit << "\t"
    // 		     << "Out of Bounds? " << (hit.quality ==
    // HitQuality::OutOfBounds)
    // 		     << endl;
    // 	}
    // }

    cout << endl
         << n_pass_cut << " events out of " << n_events << " passed cuts"
         << endl;
    cout << n_T5_valid_events << " events got a valid reconstruction -- "
         << n_pass_cut - n_T5_valid_events << " were mismatched events?" << endl
         << n_invalid_hits
         << " events were invalid -- mismatched events (the only paired SiPM "
            "hits were at totally different times)"
         << endl
         << n_events_with_valid_hits_in_expected_window
         << " events of them had a hit in the expected time window" << endl
         << n_events_out_of_bounds
         << " events had a reconstruction out of bounds" << endl
         << endl

         << n_events_with_multiple_valid_hits
         << " events had multiple valid hits -- "
         << n_events_with_multiple_valid_hits_had_one_in_expected_window
         << " of those had at least one hit in the expected time window" << endl
         << n_events_with_multiple_scint_hits
         << " events had hits in multiple scintillators -- in the expected "
            "time window"
         << endl
         << endl;
    cout << endl;

    file->Close();

    TFile *output_file = TFile::Open(current_output_path, "RECREATE");
    if (!output_file || output_file->IsZombie()) {
      cerr << "ERROR: Did not open output file: " << current_output_path
           << endl;
      continue;
    }

    output_file->cd();
    TTree *output_tree = new TTree("T5_Events", "Reconstructed T5 events");
    // --- Single-value branches (per event) ---
    int b_n_particles = 0;
    int b_event_nr = 0;
    bool b_HasValidHit;
    bool b_HasMultipleScintillatorsHit;
    bool b_HasOutOfTimeWindow;
    bool b_HasInTimeWindow;

    output_tree->Branch("event_nr", &b_event_nr, "event_nr/I");
    output_tree->Branch("T5_particle_nr", &b_n_particles, "T5_particle_nr/I");
    output_tree->Branch("T5_HasValidHit", &b_HasValidHit, "T5_HasValidHit/O");
    output_tree->Branch("T5_HasMultipleScintillatorsHit",
                        &b_HasMultipleScintillatorsHit,
                        "T5_HasMultipleScintillatorsHit/O");
    output_tree->Branch("T5_HasOutOfTimeWindow", &b_HasOutOfTimeWindow,
                        "T5_HasOutOfTimeWindow/O");
    output_tree->Branch("T5_HasInTimeWindow", &b_HasInTimeWindow,
                        "T5_HasInTimeWindow/O");

    // --- Vector branches (multiple hits per event) ---
    // Primary hits -- hits in the expected timeframe
    std::vector<int> *b_hit_is_in_bounds = new std::vector<int>();
    std::vector<double> *b_hit_pos_x = new std::vector<double>();
    std::vector<double> *b_hit_pos_y = new std::vector<double>();
    std::vector<double> *b_hit_time = new std::vector<double>();

    output_tree->Branch("T5_hit_is_in_bounds", &b_hit_is_in_bounds);
    output_tree->Branch("T5_hit_pos_x", &b_hit_pos_x);
    output_tree->Branch("T5_hit_pos_y", &b_hit_pos_y);
    output_tree->Branch("T5_hit_time", &b_hit_time);

    // Secondary hits -- hits outside of the main bunch
    std::vector<bool> *b_secondary_hit_is_in_bounds = new std::vector<bool>();
    std::vector<double> *b_secondary_hit_pos_x = new std::vector<double>();
    std::vector<double> *b_secondary_hit_pos_y = new std::vector<double>();
    std::vector<double> *b_secondary_hit_time = new std::vector<double>();

    output_tree->Branch("T5_secondary_hit_is_in_bounds",
                        &b_secondary_hit_is_in_bounds);
    output_tree->Branch("T5_secondary_hit_pos_x", &b_secondary_hit_pos_x);
    output_tree->Branch("T5_secondary_hit_pos_y", &b_secondary_hit_pos_y);
    output_tree->Branch("T5_secondary_hit_time", &b_secondary_hit_time);

    for (const auto &event : all_T5_hits) {
      b_n_particles = 0;
      b_event_nr = event.event_nr;
      b_HasValidHit = false;
      b_HasMultipleScintillatorsHit = false;
      b_HasOutOfTimeWindow = false;
      b_HasInTimeWindow = false;

      b_hit_is_in_bounds->clear();
      b_hit_pos_x->clear();
      b_hit_pos_y->clear();
      b_hit_time->clear();

      b_secondary_hit_is_in_bounds->clear();
      b_secondary_hit_pos_x->clear();
      b_secondary_hit_pos_y->clear();
      b_secondary_hit_time->clear();

      if (!event.HasValidHit) {
        output_tree->Fill();
        continue;
      }
      b_HasValidHit = event.HasValidHit;
      b_HasMultipleScintillatorsHit = event.HasMultipleScintillatorsHit;
      b_HasInTimeWindow = event.HasInTimeWindow;
      b_HasOutOfTimeWindow = event.HasOutOfTimeWindow;

      for (const auto &hit : event.T5_hits) {
        if (hit.quality == HitQuality::AccidentalCoincidence)
          continue;
        if (hit.is_in_time_window) {
          b_hit_time->push_back(hit.hit_time);
          b_hit_pos_x->push_back(hit.position_x);
          b_hit_pos_y->push_back(hit.position_y);
          if (hit.quality == HitQuality::Perfect)
            b_hit_is_in_bounds->push_back(true);
          else
            b_hit_is_in_bounds->push_back(false);
        } else {
          b_secondary_hit_time->push_back(hit.hit_time);
          b_secondary_hit_pos_x->push_back(hit.position_x);
          b_secondary_hit_pos_y->push_back(hit.position_y);
          if (hit.quality == HitQuality::Perfect)
            b_secondary_hit_is_in_bounds->push_back(true);
          else
            b_secondary_hit_is_in_bounds->push_back(false);
        }
        b_n_particles++;
      }

      output_tree->Fill();
    }
    output_tree->Write();

    delete b_hit_time;
    delete b_hit_pos_x;
    delete b_hit_pos_y;
    delete b_hit_is_in_bounds;
    delete b_secondary_hit_time;
    delete b_secondary_hit_pos_x;
    delete b_secondary_hit_pos_y;
    delete b_secondary_hit_is_in_bounds;

    output_file->Close();
  }

  return 0;
}
