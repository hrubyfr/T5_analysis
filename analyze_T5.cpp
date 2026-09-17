#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <ostream>
#include <stdexcept>
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

#include <nlohmann/json.hpp>

#include "FitParameters.h"
#include "RtypesCore.h"
#include "buffer.h"
#include "return_TOF_position.h"
#include "utils.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace ROOT;

struct run_config {
    int run_number = -1;
    string input_path = "";
    TString output_path = "";
    TString plots_directory = "";

    string config_file_path = "";
    bool is_processed = false;
    bool is_tagged_gamma = false;
    bool is_minimum_bias = false;
    bool is_charged_hadron = false;
    bool is_self_trigger = false;
    bool is_hardware_trigger = false;
};

struct Fit_results {
    bool valid = false;
    double mean_x;
    double mean_y;
    double mean_x_error;
    double mean_y_error;
    double sigma_x;
    double sigma_y;
    double sigma_x_error;
    double sigma_y_error;
    vector<int> SiPM_id;
    vector<double> SiPM_means;
    vector<double> SiPM_mean_errors;
    vector<double> SiPM_sigmas;
    vector<double> SiPM_sigma_errors;
};

TFile *open_file(run_config &run) {
    TString base_path = run.input_path;
    if (!base_path.EndsWith("/") && base_path.Length() > 0) {
        base_path += "/";
    }
    bool is_production =
        (run.input_path.find("production_v1_0") != std::string::npos);

    TString file_name;
    if (is_production) {
        file_name = TString::Format(
            "%s%d/processed_waveforms/"
            "WCTE_offline_R%dS0_VME_matched_processed_waveforms.root",
            base_path.Data(), run.run_number, run.run_number);
        run.is_processed = true;
    } else {
        file_name = TString::Format("%sWCTE_offline_R%dS0_VME_matched.root",
                                    base_path.Data(), run.run_number);
    }
    TFile *file = TFile::Open(file_name, "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: File " << file_name
                  << " could not be opened or is corrupted: " << endl;
        if (file) {
            file->Close();
        }
        throw std::runtime_error(
            Form("Could not open ROOT file: %s", file_name.Data()));
    }
    return file;
}

std::ofstream open_out_file(const run_config &run_directory) {

    TString out_file_name = run_directory.plots_directory;
    if (!out_file_name.EndsWith("/"))
        out_file_name.Append("/");
    gSystem->mkdir(out_file_name, true);
    out_file_name.Append("output.csv");

    cout << "Opening file " << out_file_name << " to output times to CSV"
         << endl;
    std::ofstream file_out(out_file_name.Data());
    if (!file_out.is_open()) {
        cerr << "Error opening csv file" << endl;
        throw std::runtime_error("Failed to open CSV file");
    }
    file_out << "RunNumber,HitTime,RawTime\n";
    return file_out;
}
void save_1D_fit_result(const Fit_results &fit_results, TString filename) {
    std::ofstream file_out_1D(filename.Data());
    if (!file_out_1D.is_open()) {
        throw std::runtime_error(
            "ERROR: file for 1D fit results could not be opened");
    }
    file_out_1D << "scint_id,mean,mean_error,sigma,sigma_error\n";
    for (int i = 0; i < fit_results.SiPM_id.size(); i++) {
        file_out_1D << fit_results.SiPM_id[i] << ","
                    << fit_results.SiPM_means[i] << ","
                    << fit_results.SiPM_mean_errors[i] << ","
                    << fit_results.SiPM_sigmas[i] << ","
                    << fit_results.SiPM_sigma_errors[i] << "\n";
    }
    file_out_1D.close();
}

void fit_positions_2D(Fit_results &fit_results, Histograms &hists,
                      const TOF_reconstructor &recon) {
    auto hist = hists.get_histogram_2D("positions");
    TF2 *gaus_2D = new TF2("gaus_2D", "bigaus", recon.Get_scint_xmin(3) * 2,
                           recon.Get_scint_xmax(3) * 2, recon.Get_ymin(),
                           recon.Get_ymax());
    gaus_2D->SetParameters(130, 0, 40, 0, 40, 0);
    hist->Fit(gaus_2D, "R");

    gaus_2D = (TF2 *)hist->GetFunction("gaus_2D");

    double volume = gaus_2D->GetParameter(0);
    double rho = gaus_2D->GetParameter(5); // Correlation factor

    auto mean_x = gaus_2D->GetParameter(1);
    auto mean_y = gaus_2D->GetParameter(3);
    auto mean_x_error = gaus_2D->GetParError(1);
    auto mean_y_error = gaus_2D->GetParError(3);

    double sig_x = gaus_2D->GetParameter(2);
    double sig_y = gaus_2D->GetParameter(4);
    auto sig_x_error = gaus_2D->GetParError(2);
    auto sig_y_error = gaus_2D->GetParError(4);

    fit_results.mean_x = mean_x;
    fit_results.mean_y = mean_y;
    fit_results.mean_x_error = mean_x_error;
    fit_results.mean_y_error = mean_y_error;

    // errors

    fit_results.sigma_x = sig_x;
    fit_results.sigma_y = sig_y;
    fit_results.sigma_x_error = sig_x_error;
    fit_results.sigma_y_error = sig_y_error;

    double chi2 = gaus_2D->GetChisquare();
    int ndf = gaus_2D->GetNDF();

    // 3. Calculate the TRUE mathematical peak height of the bigaus function
    double denominator =
        2.0 * TMath::Pi() * sig_x * sig_y * std::sqrt(1.0 - rho * rho);
    double peak_amplitude = volume / denominator;

    // 4. Define your contour levels!
    // 1-sigma drops to e^(-0.5)
    // 2-sigma drops to e^(-2.0)
    // 3-sigma drops to e^(-4.5)

    // Let's draw all 3 levels to make it look incredibly professional:
    double contours[3];
    contours[0] = peak_amplitude * std::exp(-4.5); // 3-sigma (widest, lowest)
    contours[1] = peak_amplitude * std::exp(-2.0); // 2-sigma
    contours[2] =
        peak_amplitude * std::exp(-0.5); // 1-sigma (tightest, highest)

    // 5. Apply the contours to your TF2
    // The arguments are: (number_of_levels, array_of_levels)
    gaus_2D->SetContour(3, contours);

    // 6. Make the contour lines stand out against the color map
    gaus_2D->SetLineColor(kRed);
    gaus_2D->SetLineWidth(2);
    gaus_2D->SetLineStyle(1); // Solid lines

    TString txt_sig_x =
        Form("#sigma_{x} [mm] = \n (%.2f #pm %.2f)", sig_x, sig_x_error);
    TString txt_sig_y =
        Form("#sigma_{y} [mm] = \n (%.2f #pm %.2f)", sig_y, sig_y_error);
    TString txt_rho = Form("#rho = %.2f", rho);
    TString txt_chi2_ndf = Form("#chi^{2}/NDF = %.2f/%i", chi2, ndf);

    double offset = 10.0;
    TLatex *ltx_sigX = new TLatex(sig_x - offset, sig_y, txt_sig_x);
    TLatex *ltx_sigY = new TLatex(sig_x - offset, sig_y - 10.0, txt_sig_y);
    double placement_rho_x = -40.0;
    double placement_rho_y = 45.0;
    TLatex *ltx_rho = new TLatex(placement_rho_x, placement_rho_y, txt_rho);
    TLatex *ltx_chi2_ndf =
        new TLatex(placement_rho_x, placement_rho_y - 7.5, txt_chi2_ndf);

    ltx_sigX->SetTextSize(0.04);
    ltx_sigX->SetTextColor(kBlack);
    ltx_sigY->SetTextSize(0.04);
    ltx_sigY->SetTextColor(kBlack);
    ltx_rho->SetTextSize(0.04);
    ltx_rho->SetTextColor(kBlack);
    ltx_chi2_ndf->SetTextSize(0.04);
    ltx_chi2_ndf->SetTextColor(kBlack);

    TLatex *contour_sigma = new TLatex(mean_x, -sig_y, "1#sigma");
    contour_sigma->SetTextSize(0.04);
    contour_sigma->SetTextColor(kBlack);
    TLatex *contour_sigma_2 =
        new TLatex(mean_x, mean_y + (-2 * sig_y), "2#sigma");
    contour_sigma_2->SetTextSize(0.04);
    contour_sigma_2->SetTextColor(kBlack);
    TLatex *contour_sigma_3 =
        new TLatex(mean_x, mean_y + (-3 * sig_y), "3#sigma");
    contour_sigma_3->SetTextSize(0.04);
    contour_sigma_3->SetTextColor(kBlack);

    hist->GetListOfFunctions()->Add(ltx_sigX);
    hist->GetListOfFunctions()->Add(ltx_sigY);
    hist->GetListOfFunctions()->Add(ltx_rho);
    hist->GetListOfFunctions()->Add(ltx_chi2_ndf);
    hist->GetListOfFunctions()->Add(contour_sigma);
    hist->GetListOfFunctions()->Add(contour_sigma_2);
    hist->GetListOfFunctions()->Add(contour_sigma_3);
}

void fit_positions_1D(Fit_results &fit_results, Histograms &hists,
                      TOF_reconstructor recon) {
    for (int i = 0; i < 8; i++) {
        TString h_name = "positions_" + std::to_string(i);
        auto fit_fun = new TF1("gaussian", "gaus", 2 * recon.Get_scint_xmin(i),
                               2 * recon.Get_scint_xmax(i));
        hists.hist_projectX("positions", h_name.Data(), i + 1, i + 1);
        hists.get_histogram(h_name.Data())->Fit(fit_fun, "R");
        auto mean = fit_fun->GetParameter(1);
        auto mean_error = fit_fun->GetParError(1);
        auto sigma = fit_fun->GetParameter(2);
        auto sigma_error = fit_fun->GetParError(2);
        fit_results.SiPM_means.push_back(mean);
        fit_results.SiPM_mean_errors.push_back(mean_error);
        fit_results.SiPM_sigmas.push_back(sigma);
        fit_results.SiPM_sigma_errors.push_back(sigma_error);
    }
}

void analyze_tagged_gamma(TFile *file, run_config run_directory,
                          Fit_results fit_results) {}

void analyze_processed(TFile *file, run_config run_directory,
                       Fit_results &fit_results) {

    cout << "Loading Tree..." << endl;

    auto tree = file->Get<TTree>("ProcessedWaveforms");

    cout << "tree loaded, printing structure  ..." << endl;

    tree->Print();

    cout << "Finished printing tree structure, setting branch addresses ..."
         << endl;

    // vector<float> *arr_bm_times = nullptr;
    // vector<float> *arr_bm_charges = nullptr;
    // vector<int> *arr_bm_time_ids = nullptr;
    // vector<int> *arr_bm_charge_ids = nullptr;
    Int_t nhits;
    const int MAX_HITS = 5000;
    Double_t hit_time_carray[MAX_HITS];
    Double_t hit_charge_carray[MAX_HITS];
    Int_t hit_pmt_carray[MAX_HITS];
    Int_t hit_card_carray[MAX_HITS];

    vector<double> *arr_pmt_times = nullptr;
    vector<double> *arr_pmt_charges = nullptr;
    vector<int> *arr_pmt_ids = nullptr;
    vector<int> *arr_mpmt_ids = nullptr;

    // vector<vector<double>> *pmt_waveforms = nullptr;
    // vector<double> *pmt_waveform_times = nullptr;
    // vector<int> *pmt_waveform_card_ids = nullptr;
    // vector<int> *pmt_waveform_pmt_ids = nullptr;

    TString b_hit_times_name, b_hit_charges_name, b_hit_channel_name,
        b_hit_card_name;
    if (run_directory.is_processed) {
        b_hit_times_name = "hit_time";
        b_hit_charges_name = "hit_charge";
        b_hit_channel_name = "hit_chan";
        b_hit_card_name = "hit_card";
    } else {
        b_hit_times_name = "hit_pmt_times";
        b_hit_charges_name = "hit_pmt_charges";
        b_hit_channel_name = "hit_pmt_channel_ids";
        b_hit_card_name = "hit_mpmt_card_ids";
    }
    TBranch *b_nhits = nullptr;

    tree->SetBranchStatus("*", 0);
    // tree->SetBranchStatus("beamline_pmt_qdc_ids", 1);
    // tree->SetBranchStatus("beamline_pmt_tdc_ids", 1);
    // tree->SetBranchStatus("beamline_pmt_tdc_times", 1);
    // tree->SetBranchStatus("beamline_pmt_qdc_charges", 1);
    tree->SetBranchStatus(b_hit_times_name, 1);
    tree->SetBranchStatus(b_hit_channel_name, 1);
    tree->SetBranchStatus(b_hit_charges_name, 1);
    tree->SetBranchStatus(b_hit_card_name, 1);

    if (run_directory.is_processed) {
        tree->SetBranchStatus("nhit_time", 1);
        tree->SetBranchAddress("nhit_time", &nhits, &b_nhits);
    }
    // tree->SetBranchStatus("pmt_waveforms", 1);
    // tree->SetBranchStatus("pmt_waveform_times", 1);
    // tree->SetBranchStatus("pmt_waveform_mpmt_card_ids", 1);
    // tree->SetBranchStatus("pmt_waveform_pmt_channel_ids", 1);

    // tree->SetBranchAddress("beamline_pmt_qdc_ids", &arr_bm_charge_ids);
    // tree->SetBranchAddress("beamline_pmt_tdc_ids", &arr_bm_time_ids);
    // tree->SetBranchAddress("beamline_pmt_tdc_times", &arr_bm_times);
    // tree->SetBranchAddress("beamline_pmt_qdc_charges", &arr_bm_charges);

    if (run_directory.is_processed) {

        tree->SetBranchAddress(b_hit_times_name, &hit_time_carray);
        tree->SetBranchAddress(b_hit_charges_name, &hit_charge_carray);
        tree->SetBranchAddress(b_hit_card_name, &hit_card_carray);
        tree->SetBranchAddress(b_hit_channel_name, &hit_pmt_carray);
    } else {
        tree->SetBranchAddress(b_hit_times_name, &arr_pmt_times);
        tree->SetBranchAddress(b_hit_charges_name, &arr_pmt_charges);
        tree->SetBranchAddress(b_hit_card_name, &arr_mpmt_ids);
        tree->SetBranchAddress(b_hit_channel_name, &arr_pmt_ids);
    }
    // tree->SetBranchAddress("pmt_waveforms", &pmt_waveforms);
    // tree->SetBranchAddress("pmt_waveform_times", &pmt_waveform_times);
    // tree->SetBranchAddress("pmt_waveform_mpmt_card_ids",
    //                        &pmt_waveform_card_ids);
    // tree->SetBranchAddress("pmt_waveform_pmt_channel_ids",
    //                        &pmt_waveform_pmt_ids);

    Cuts cut;
    TOF_reconstructor recon;
    Histograms hists;
    setup_histograms(hists, recon);
    int n_pass_cut = 0;
    int n_T5_valid_events = 0;
    auto n_events = tree->GetEntries();
    int verb = 1000;
    int n_events_with_multiple_valid_hits = 0;
    int n_events_with_valid_hits_in_expected_window = 0;
    int n_events_with_multiple_scint_hits = 0;
    int n_events_with_multiple_valid_hits_had_one_in_expected_window = 0;
    int n_invalid_hits = 0;
    int n_events_out_of_bounds = 0;

    vector<int> hit_intervals;
    for (int time_low = -65; time_low < 36; time_low += 10)
        hit_intervals.push_back(time_low);

    int n_T5_scintillators = 8;

    vector<event_T5_detection> all_T5_hits;
    vector<int> n_hits_in_scints(n_T5_scintillators);
    for (auto &hit : n_hits_in_scints) {
        hit = 0;
    }
    if (run_directory.is_processed) {

        arr_pmt_times = new std::vector<double>();
        arr_pmt_charges = new std::vector<double>();
        arr_pmt_ids = new std::vector<int>();
        arr_mpmt_ids = new std::vector<int>();
    }

    // Save SiPM times to a temporary CSV file, to be later parsed by a python
    // script
    auto file_out = open_out_file(run_directory);

    cout << "Starting event loop over " << n_events << " events..." << endl;

    for (long long i = 0; i < n_events; i++) {
        if (run_directory.is_processed) {
            auto entry = tree->LoadTree(i);
            if (entry < 0)
                break;
            b_nhits->GetEntry(i);
            if (nhits > MAX_HITS || nhits < 0) {
                continue;
            }
            arr_pmt_times->assign(hit_time_carray, hit_time_carray + nhits);
            arr_pmt_charges->assign(hit_charge_carray,
                                    hit_charge_carray + nhits);
            arr_pmt_ids->assign(hit_pmt_carray, hit_pmt_carray + nhits);
            arr_mpmt_ids->assign(hit_card_carray, hit_card_carray + nhits);
        }

        tree->GetEntry(i);
        // Print progress
        event_T5_detection detections;

        if (i % verb == 0)
            cout << "\rAnalyzed " << i << " of " << n_events
                 << std::setprecision(2) << std::fixed << " events ("
                 << static_cast<float>(i) / n_events * 100 << " %)"
                 << std::flush << endl;

        // RVecI bm_time_ids(arr_bm_time_ids->data(), arr_bm_time_ids->size());
        // RVecI bm_charge_ids(arr_bm_charge_ids->data(),
        //                     arr_bm_charge_ids->size());
        // RVecF bm_times(arr_bm_times->data(), arr_bm_times->size());
        // RVecF bm_charges(arr_bm_charges->data(), arr_bm_charges->size());

        RVecD pmt_times(arr_pmt_times->data(), arr_pmt_times->size());
        RVecI pmt_ids(arr_pmt_ids->data(), arr_pmt_ids->size());
        RVecI mpmt_ids(arr_mpmt_ids->data(), arr_mpmt_ids->size());

        if (!cut.hit_T5(mpmt_ids, pmt_ids)) {
            detections.event_nr = i;
            all_T5_hits.push_back(detections);
            continue;
        }

        n_pass_cut++;

        auto mask_T5_board = (mpmt_ids == cut.get_T5_board());
        auto T5_board_ids = pmt_ids[mask_T5_board];
        auto T5_board_times = pmt_times[mask_T5_board];

        // RVecI wf_mpmt_ids(pmt_waveform_card_ids->data(),
        //                   pmt_waveform_card_ids->size());
        // RVecI wf_pmt_ids(pmt_waveform_pmt_ids->data(),
        //                  pmt_waveform_pmt_ids->size());
        // RVecD wf_start_times(pmt_waveform_times->data(),
        //                      pmt_waveform_times->size());
        // RVec<std::vector<double>> wf_pmt_waveforms(pmt_waveforms->data(),
        //                                            pmt_waveforms->size());
        //
        // auto wf_T5_card_positions = VecOps::Nonzero(wf_mpmt_ids == 132);
        // auto wf_pmts_filtered = VecOps::Take(wf_pmt_ids,
        // wf_T5_card_positions); wf_start_times = VecOps::Take(wf_start_times,
        // wf_T5_card_positions); wf_pmt_waveforms =
        // VecOps::Take(wf_pmt_waveforms, wf_T5_card_positions);
        //
        // auto lmbd_is_T5 = [wf_pmts_filtered, &cut]() {
        //     RVecI mask(wf_pmts_filtered.size(), 0);
        //     for (const auto &id : cut.Get_T5_ids()) {
        //         mask = mask || (wf_pmts_filtered == id);
        //     }
        //     return mask;
        // };
        // auto t5_mask = lmbd_is_T5();
        // wf_pmts_filtered = wf_pmts_filtered[t5_mask];
        // wf_start_times = wf_start_times[t5_mask];
        // wf_pmt_waveforms = wf_pmt_waveforms[t5_mask];
        detections = recon.Return_position(i, *arr_mpmt_ids, *arr_pmt_ids,
                                           *arr_pmt_times, *arr_pmt_charges);

        if (detections.HasValidHit) {
            n_T5_valid_events++;
            for (const auto &hit : detections.T5_hits) {
                if (!hit.is_valid_hit)
                    continue;
                hists.fill("valid_hit_times", hit.hit_time);
                hists.fill(Form("hit_raw_times_%i", hit.scintillator_id),
                           hit.raw_time);
                hists.fill("trigger_times", hit.trigger_time);
                file_out << run_directory.run_number << "," << hit.hit_time
                         << "," << hit.raw_time << "\n";
            }
        }
        // cout << "Event " << i << ": " << detections.T5_hits.size() << " hits"
        //      << endl;
        // for (const auto &hit : detections.T5_hits) {
        //     cout << "\tSiPM: " << hit.scintillator_id
        //          << " Time: " << hit.hit_time << " X: " << hit.position_x
        //          << " Y: " << hit.position_y;
        //     if (hit.is_valid_hit)
        //         cout << " Valid ";
        //     cout << endl;
        // }

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
        hists.fill("n_event_hits", n_hits_in_T5_in_single_event);

        for (const auto &hit : detections.T5_hits) {
            bool is_in_time_window = true;

            if (!hit.is_valid_hit)
                //|| hit.quality != HitQuality::Perfect)
                continue;
            // if (hit.sipm_a_charge < 1000 || hit.sipm_b_charge < 1000)
            //     continue;
            hists.fill(Form("hit_charges_2D_%i", hit.scintillator_id),
                       hit.sipm_a_charge, hit.sipm_b_charge);
            hists.fill(Form("hit_charges_%i", hit.scintillator_id),
                       hit.total_hit_charge);
            hists.fill("positions", hit.position_x, hit.position_y);
            n_hits_in_scints.at(hit.scintillator_id)++;
        }

        for (auto &hit : detections.T5_hits) {
            auto sipm_a_Q_threshold = 1000;
            auto sipm_b_Q_threshold = 1000;
            if (hit.scintillator_id == 2)
                sipm_b_Q_threshold = 400;
            if (hit.sipm_a_charge < sipm_a_Q_threshold ||
                hit.sipm_b_charge < sipm_b_Q_threshold || !hit.is_valid_hit)
                continue;
            hists.fill("positions_chargecut", hit.position_x, hit.position_y);
        }

        std::unordered_set<double> bad_timestamps;

        for (auto &hit : detections.T5_hits) {

            if (hit.sipm_time_a == hit.sipm_time_b) {
                // cout << "Suspicious event! the detection times are exactly
                // the "
                //         "same!"
                // << endl;
                bad_timestamps.insert(hit.sipm_time_a);
            }
        }
        if (!bad_timestamps.empty()) {
            detections.T5_hits.erase(
                std::remove_if(
                    detections.T5_hits.begin(), detections.T5_hits.end(),
                    [&bad_timestamps](const auto &hit) {
                        bool has_bad_time_a =
                            (bad_timestamps.count(hit.sipm_time_a) > 0);
                        bool has_bad_time_b =
                            (bad_timestamps.count(hit.sipm_time_b) > 0);
                        return has_bad_time_a || has_bad_time_b;
                    }),
                detections.T5_hits.end());
        }
        for (const auto &hit : detections.T5_hits) {
            if (!hit.is_valid_hit)
                continue;
            hists.fill("positions_timecut", hit.position_x, hit.position_y);

            auto sipm_a_Q_threshold = 1000;
            auto sipm_b_Q_threshold = 1000;
            if (hit.scintillator_id == 2)
                sipm_b_Q_threshold = 400;
            if (hit.total_hit_charge < 1000)
                continue;
            hists.fill("positions_timechargecut", hit.position_x,
                       hit.position_y);
        }

        all_T5_hits.push_back(detections);
    }

    fit_positions_2D(fit_results, hists, recon);

    fit_positions_1D(fit_results, hists, recon);

    TString fit_result_1D_filename = run_directory.plots_directory;
    if (!fit_result_1D_filename.EndsWith("/"))
        fit_result_1D_filename.Append("/");
    fit_result_1D_filename.Append("fit_1D.csv");
    save_1D_fit_result(fit_results, fit_result_1D_filename);

    if (run_directory.plots_directory == "") {
        run_directory.plots_directory =
            "plots/Run_" + std::to_string(run_directory.run_number);
    }

    gSystem->Exec("mkdir -p " + run_directory.plots_directory);
    gSystem->cd(run_directory.plots_directory);
    hists.print_exclusive("positions", 1000, 900);
    hists.print_all();
    hists.save_all("hists");
    hists.print_exclusive_log("positions", 1000, 900);
    hists.print_exclusive_log("valid_hit_times", 1800, 900);
    hists.print_exclusive("positions_chargecut", 1000, 900);
    hists.print_exclusive("positions_timecut", 1000, 900);
    hists.print_exclusive("positions_timechargecut", 1000, 900);

    cout << endl
         << n_pass_cut << " events out of " << n_events << " passed cuts"
         << endl;
    cout
        << n_T5_valid_events << " events got a valid reconstruction -- "
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
        << " events had hits in multiple scintillators -- in the expected time "
           "window"
        << endl
        << endl;
    cout << endl;
    gSystem->cd("/eos/user/f/fhruby/projects/T5_analysis/");
}

void analyze_raw(TFile *file, run_config run_directory,
                 Fit_results &fit_results) {

    cout << "Loading Tree..." << endl;

    auto tree = file->Get<TTree>("WCTEReadoutWindows");

    cout << "tree loaded, printing structure  ..." << endl;

    tree->Print();

    cout << "Finished printing tree structure, setting branch addresses ..."
         << endl;

    // vector<float> *arr_bm_times = nullptr;
    // vector<float> *arr_bm_charges = nullptr;
    // vector<int> *arr_bm_time_ids = nullptr;
    // vector<int> *arr_bm_charge_ids = nullptr;

    vector<double> *arr_pmt_times = nullptr;
    vector<double> *arr_pmt_charges = nullptr;
    vector<int> *arr_pmt_ids = nullptr;
    vector<int> *arr_mpmt_ids = nullptr;

    // vector<vector<double>> *pmt_waveforms = nullptr;
    // vector<double> *pmt_waveform_times = nullptr;
    // vector<int> *pmt_waveform_card_ids = nullptr;
    // vector<int> *pmt_waveform_pmt_ids = nullptr;

    tree->SetBranchStatus("*", 0);
    // tree->SetBranchStatus("beamline_pmt_qdc_ids", 1);
    // tree->SetBranchStatus("beamline_pmt_tdc_ids", 1);
    // tree->SetBranchStatus("beamline_pmt_tdc_times", 1);
    // tree->SetBranchStatus("beamline_pmt_qdc_charges", 1);
    tree->SetBranchStatus("hit_pmt_times", 1);
    tree->SetBranchStatus("hit_pmt_charges", 1);
    tree->SetBranchStatus("hit_mpmt_card_ids", 1);
    tree->SetBranchStatus("hit_pmt_channel_ids", 1);

    // tree->SetBranchStatus("pmt_waveforms", 1);
    // tree->SetBranchStatus("pmt_waveform_times", 1);
    // tree->SetBranchStatus("pmt_waveform_mpmt_card_ids", 1);
    // tree->SetBranchStatus("pmt_waveform_pmt_channel_ids", 1);

    // tree->SetBranchAddress("beamline_pmt_qdc_ids", &arr_bm_charge_ids);
    // tree->SetBranchAddress("beamline_pmt_tdc_ids", &arr_bm_time_ids);
    // tree->SetBranchAddress("beamline_pmt_tdc_times", &arr_bm_times);
    // tree->SetBranchAddress("beamline_pmt_qdc_charges", &arr_bm_charges);
    tree->SetBranchAddress("hit_pmt_times", &arr_pmt_times);
    tree->SetBranchAddress("hit_pmt_charges", &arr_pmt_charges);
    tree->SetBranchAddress("hit_mpmt_card_ids", &arr_mpmt_ids);
    tree->SetBranchAddress("hit_pmt_channel_ids", &arr_pmt_ids);

    // tree->SetBranchAddress("pmt_waveforms", &pmt_waveforms);
    // tree->SetBranchAddress("pmt_waveform_times", &pmt_waveform_times);
    // tree->SetBranchAddress("pmt_waveform_mpmt_card_ids",
    //                        &pmt_waveform_card_ids);
    // tree->SetBranchAddress("pmt_waveform_pmt_channel_ids",
    //                        &pmt_waveform_pmt_ids);

    Cuts cut;
    TOF_reconstructor recon;
    Histograms hists;
    setup_histograms(hists, recon);
    int n_pass_cut = 0;
    int n_T5_valid_events = 0;
    auto n_events = tree->GetEntries();
    int verb = 1000;
    int n_events_with_multiple_valid_hits = 0;
    int n_events_with_valid_hits_in_expected_window = 0;
    int n_events_with_multiple_scint_hits = 0;
    int n_events_with_multiple_valid_hits_had_one_in_expected_window = 0;
    int n_invalid_hits = 0;
    int n_events_out_of_bounds = 0;

    vector<int> hit_intervals;
    for (int time_low = -65; time_low < 36; time_low += 10)
        hit_intervals.push_back(time_low);

    int n_T5_scintillators = 8;

    vector<event_T5_detection> all_T5_hits;
    vector<int> n_hits_in_scints(n_T5_scintillators);
    for (auto &hit : n_hits_in_scints) {
        hit = 0;
    }

    auto file_out = open_out_file(run_directory);

    cout << "Starting event loop over " << n_events << " events..." << endl;

    for (long long i = 0; i < n_events; i++) {
        tree->GetEntry(i);
        // Print progress
        event_T5_detection detections;

        if (i % verb == 0)
            cout << "\rAnalyzed " << i << " of " << n_events
                 << std::setprecision(2) << std::fixed << " events ("
                 << static_cast<float>(i) / n_events * 100 << " %)"
                 << std::flush << endl;

        // RVecI bm_time_ids(arr_bm_time_ids->data(), arr_bm_time_ids->size());
        // RVecI bm_charge_ids(arr_bm_charge_ids->data(),
        //                     arr_bm_charge_ids->size());
        // RVecF bm_times(arr_bm_times->data(), arr_bm_times->size());
        // RVecF bm_charges(arr_bm_charges->data(), arr_bm_charges->size());

        RVecD pmt_times(arr_pmt_times->data(), arr_pmt_times->size());
        RVecI pmt_ids(arr_pmt_ids->data(), arr_pmt_ids->size());
        RVecI mpmt_ids(arr_mpmt_ids->data(), arr_mpmt_ids->size());

        if (!cut.hit_T5(mpmt_ids, pmt_ids)) {
            detections.event_nr = i;
            all_T5_hits.push_back(detections);
            continue;
        }

        n_pass_cut++;

        auto mask_T5_board = (mpmt_ids == cut.get_T5_board());
        auto T5_board_ids = pmt_ids[mask_T5_board];
        auto T5_board_times = pmt_times[mask_T5_board];

        // RVecI wf_mpmt_ids(pmt_waveform_card_ids->data(),
        //                   pmt_waveform_card_ids->size());
        // RVecI wf_pmt_ids(pmt_waveform_pmt_ids->data(),
        //                  pmt_waveform_pmt_ids->size());
        // RVecD wf_start_times(pmt_waveform_times->data(),
        //                      pmt_waveform_times->size());
        // RVec<std::vector<double>> wf_pmt_waveforms(pmt_waveforms->data(),
        //                                            pmt_waveforms->size());
        //
        // auto wf_T5_card_positions = VecOps::Nonzero(wf_mpmt_ids == 132);
        // auto wf_pmts_filtered = VecOps::Take(wf_pmt_ids,
        // wf_T5_card_positions); wf_start_times = VecOps::Take(wf_start_times,
        // wf_T5_card_positions); wf_pmt_waveforms =
        // VecOps::Take(wf_pmt_waveforms, wf_T5_card_positions);
        //
        // auto lmbd_is_T5 = [wf_pmts_filtered, &cut]() {
        //     RVecI mask(wf_pmts_filtered.size(), 0);
        //     for (const auto &id : cut.Get_T5_ids()) {
        //         mask = mask || (wf_pmts_filtered == id);
        //     }
        //     return mask;
        // };
        // auto t5_mask = lmbd_is_T5();
        // wf_pmts_filtered = wf_pmts_filtered[t5_mask];
        // wf_start_times = wf_start_times[t5_mask];
        // wf_pmt_waveforms = wf_pmt_waveforms[t5_mask];
        detections = recon.Return_position(i, *arr_mpmt_ids, *arr_pmt_ids,
                                           *arr_pmt_times, *arr_pmt_charges);

        if (detections.HasValidHit) {
            n_T5_valid_events++;
            for (const auto &hit : detections.T5_hits) {
                if (!hit.is_valid_hit)
                    continue;
                hists.fill("valid_hit_times", hit.hit_time);
                hists.fill(Form("hit_raw_times_%i", hit.scintillator_id),
                           hit.raw_time);
                hists.fill("trigger_times", hit.trigger_time);
                file_out << run_directory.run_number << "," << hit.hit_time
                         << "," << hit.raw_time << "\n";
            }
        }
        // cout << "Event " << i << ": " << detections.T5_hits.size() << " hits"
        //      << endl;
        // for (const auto &hit : detections.T5_hits) {
        //     cout << "\tSiPM: " << hit.scintillator_id
        //          << " Time: " << hit.hit_time << " X: " << hit.position_x
        //          << " Y: " << hit.position_y;
        //     if (hit.is_valid_hit)
        //         cout << " Valid ";
        //     cout << endl;
        // }

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
        hists.fill("n_event_hits", n_hits_in_T5_in_single_event);

        for (const auto &hit : detections.T5_hits) {
            bool diff_trigger = false;
            if (run_directory.run_number > 2067 &&
                run_directory.run_number < 2220) {
                diff_trigger = true;
            }
            bool is_in_time_window = true;
            if (diff_trigger) {
                double MAX_TIME = 20;
                double MIN_TIME = -10;
                is_in_time_window =
                    (hit.hit_time > MIN_TIME && hit.hit_time < MAX_TIME);
            }

            if (!hit.is_valid_hit)
                //|| hit.quality != HitQuality::Perfect)
                continue;
            // if (hit.sipm_a_charge < 1000 || hit.sipm_b_charge < 1000)
            //     continue;
            hists.fill(Form("hit_charges_2D_%i", hit.scintillator_id),
                       hit.sipm_a_charge, hit.sipm_b_charge);
            hists.fill(Form("hit_charges_%i", hit.scintillator_id),
                       hit.total_hit_charge);
            hists.fill("positions", hit.position_x, hit.position_y);
            n_hits_in_scints.at(hit.scintillator_id)++;
        }

        for (auto &hit : detections.T5_hits) {
            auto sipm_a_Q_threshold = 1000;
            auto sipm_b_Q_threshold = 1000;
            if (hit.scintillator_id == 2)
                sipm_b_Q_threshold = 400;
            if (hit.sipm_a_charge < sipm_a_Q_threshold ||
                hit.sipm_b_charge < sipm_b_Q_threshold || !hit.is_valid_hit)
                continue;
            hists.fill("positions_chargecut", hit.position_x, hit.position_y);
        }

        std::unordered_set<double> bad_timestamps;

        for (auto &hit : detections.T5_hits) {

            if (hit.sipm_time_a == hit.sipm_time_b) {
                // cout << "Suspicious event! the detection times are exactly
                // the "
                //         "same!"
                // << endl;
                bad_timestamps.insert(hit.sipm_time_a);
            }
        }
        if (!bad_timestamps.empty()) {
            detections.T5_hits.erase(
                std::remove_if(
                    detections.T5_hits.begin(), detections.T5_hits.end(),
                    [&bad_timestamps](const auto &hit) {
                        bool has_bad_time_a =
                            (bad_timestamps.count(hit.sipm_time_a) > 0);
                        bool has_bad_time_b =
                            (bad_timestamps.count(hit.sipm_time_b) > 0);
                        return has_bad_time_a || has_bad_time_b;
                    }),
                detections.T5_hits.end());
        }
        for (const auto &hit : detections.T5_hits) {
            if (!hit.is_valid_hit)
                continue;
            hists.fill("positions_timecut", hit.position_x, hit.position_y);

            auto sipm_a_Q_threshold = 1000;
            auto sipm_b_Q_threshold = 1000;
            if (hit.scintillator_id == 2)
                sipm_b_Q_threshold = 400;
            if (hit.total_hit_charge < 1000)
                continue;
            hists.fill("positions_timechargecut", hit.position_x,
                       hit.position_y);
        }

        all_T5_hits.push_back(detections);
    }
    fit_positions_2D(fit_results, hists, recon);
    fit_positions_1D(fit_results, hists, recon);

    TString fit_result_1D_filename = run_directory.plots_directory;
    if (!fit_result_1D_filename.EndsWith("/"))
        fit_result_1D_filename.Append("/");
    fit_result_1D_filename.Append("fit_1D.csv");
    save_1D_fit_result(fit_results, fit_result_1D_filename);

    if (run_directory.plots_directory == "") {
        run_directory.plots_directory =
            "plots/Run_" + std::to_string(run_directory.run_number);
    }

    gSystem->Exec("mkdir -p " + run_directory.plots_directory);
    gSystem->cd(run_directory.plots_directory);
    hists.print_exclusive("positions", 1000, 900);
    hists.print_all();
    hists.save_all("hists");
    hists.print_exclusive_log("positions", 1000, 900);
    hists.print_exclusive_log("valid_hit_times", 1800, 900);
    hists.print_exclusive("positions_chargecut", 1000, 900);
    hists.print_exclusive("positions_timecut", 1000, 900);
    hists.print_exclusive("positions_timechargecut", 1000, 900);

    cout << endl
         << n_pass_cut << " events out of " << n_events << " passed cuts"
         << endl;
    cout
        << n_T5_valid_events << " events got a valid reconstruction -- "
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
        << " events had hits in multiple scintillators -- in the expected time "
           "window"
        << endl
        << endl;
    cout << endl;
    gSystem->cd("/eos/user/f/fhruby/projects/T5_analysis/");
}

int main(int argc, char **argv) {

    run_config run_directory;

    int opt;
    while ((opt = getopt(argc, argv, "r:i:o:p:c:")) != -1) {
        switch (opt) {
        case 'r':
            run_directory.run_number = std::stoi(optarg);
            break;
        case 'i':
            run_directory.input_path = optarg;
            break;
        case 'o':
            run_directory.output_path = optarg;
            break;
        case 'p':
            run_directory.plots_directory = optarg;
            break;
        case 'c':
            run_directory.config_file_path = optarg;
            break;
        default:
            cerr << "Usage: " << argv[0]
                 << " -r <run_number> [-i <input_path>] [-o <output_file>]"
                 << endl;
            return -1;
        }
    }

    RUN_NUMBER = run_directory.run_number;

    if (run_directory.config_file_path == "") {
        string config_file_path = "../configs/LEMB_runs.json";
    }
    cout << "Opening a config file " << run_directory.config_file_path << endl;
    std::ifstream config_file(run_directory.config_file_path);
    if (!config_file || !config_file.is_open()) {
        cerr << "ERROR: could not open config file" << endl;
        return 1;
    }

    nlohmann::json config;
    config_file >> config;

    cout << "Loading run configuration from json...";
    nlohmann::json *run_config;
    for (auto &run : config) {
        if (run["run_number"] == "")
            continue;
        cout << "Run " << run["run_number"] << "...";
        string run_nr_str = run["run_number"];

        if (std::stoi(run_nr_str) == run_directory.run_number) {
            run_config = &run;
            string run_momentum_str = run.value("beam_momentum", "0");
            string run_configuration = run.value("run_config", "");
            string beam_config = run.value("beam_config", "");
            string trigger_config = run.value("trigger_config", "");
            if (run_configuration.find("mpmt_beam") != std::string::npos)
                run_directory.is_self_trigger = true;
            else if (run_configuration.find("hardware_trigger") !=
                     std::string::npos)
                run_directory.is_hardware_trigger = true;
            if (beam_config.find("hadron") != std::string::npos)
                run_directory.is_charged_hadron = true;
            else if (beam_config.find("tagged gamma") != std::string::npos)
                run_directory.is_tagged_gamma = true;
            if (trigger_config.find("LEMB") != std::string::npos)
                run_directory.is_minimum_bias = true;
            BEAM_MOMENTUM = std::stoi(run_momentum_str);
            break;
        }
    }
    cout << "Found\t Done, closing config file" << endl;
    config_file.close();
    cout << "Config file loaded, beam momentum is " << BEAM_MOMENTUM << "MeV/c"
         << endl;

    TFile *file;
    Fit_results fit_results;

    if (run_directory.is_self_trigger) {
        cout << "Run was diagnosed as self-trigger, running self trigger "
                "analysis"
             << endl;
        run_directory.input_path = "/eos/experiment/wcte/data/"
                                   "2025_commissioning/offline_data_vme_match/";

        file = open_file(run_directory);
        analyze_raw(file, run_directory, fit_results);
    } else if (run_directory.is_hardware_trigger) {
        cout << "Run was diagnosed as hardware_trigger, running processed "
                "analysis"
             << endl;
        analyze_processed(file, run_directory, fit_results);
    }

    file->Close();

    // Calculate the effective increase in width of the measured sigma_x due to
    // SiPM resolution

    // Load fit values from T5 model fit, extract SiPM resolutions

    // cout << " Calculating correction to x beam width: " << endl;
    //
    // cout << " Loading resolutions from a file" << endl;
    // FitParameters T5_parameters;
    // cout << " Loaded effective speed, it is " << T5_parameters.Get_veff()
    //      << " cm/ns" << endl;
    // auto resolutions = T5_parameters.GetResolutions();
    // cout << "Done, " << resolutions.size()
    //      << " resolutions extracted, they are: " << endl;
    // for (const auto &res : resolutions) {
    //   cout << "SiPM " << res.sipm_nr << ": " << res.width_cm << " cm ("
    //        << res.time_ns << " ns)" << endl;
    // }
    // cout << "Accumulating the total number of events" << endl;
    // auto total_events =
    //     std::accumulate(n_hits_in_scints.begin(), n_hits_in_scints.end(), 0);
    // cout << " Done, the total number of events is " << total_events << endl;
    // cout << "Calculating the variance of effective sigma..." << endl;
    // double effective_sigma_sipms_mm_square = 0;
    // for (int i = 0; i < n_T5_scintillators; i++) {
    //   auto weight = static_cast<double>(n_hits_in_scints.at(i)) /
    //   total_events; effective_sigma_sipms_mm_square +=
    //       weight * pow(resolutions.at(i).width_cm * 10, 2);
    // }
    // auto effective_sigma_sipms_mm = sqrt(effective_sigma_sipms_mm_square);
    // cout << "Done, it is " << effective_sigma_sipms_mm_square << endl;
    //
    // double variance_error_square = 0;
    // for (int i = 0; i < n_T5_scintillators; i++) {
    //   double r_i_sq = pow(resolutions.at(i).width_cm, 2);
    //   double n_hits = static_cast<double>(n_hits_in_scints.at(i));
    //
    //   // Propagate the Poisson error of the hits (delta N_i = sqrt(N_i))
    //   variance_error_square +=
    //       n_hits *
    //       pow((r_i_sq - effective_sigma_sipms_mm_square) / total_events, 2);
    // }
    // auto delta_V = sqrt(variance_error_square);
    //
    // auto sig_x_corr = sqrt(pow(sig_x, 2) - effective_sigma_sipms_mm_square);
    // auto sig_x_corr_error = delta_V / (2.0 * effective_sigma_sipms_mm);
    // pow(effective_sigma_sipms * ))

    // (*run_config)["T5_beam_sigma_x"] = fit_results.sigma_x;
    // (*run_config)["T5_beam_sigma_y"] = fit_results.sigma_y;
    // (*run_config)["T5_beam_sigma_x_error"] = fit_results.sigma_x_error;
    // (*run_config)["T5_beam_sigma_y_error"] = fit_results.sigma_y_error;
    // (*run_config)["T5_beam_mean_x"] = fit_results.mean_x;
    // (*run_config)["T5_beam_mean_y"] = fit_results.mean_y;
    // (*run_config)["T5_beam_mean_x_error"] = fit_results.mean_x_error;
    // (*run_config)["T5_beam_mean_y_error"] = fit_results.mean_y_error;
    //
    // string sigmas_output;
    // if (run_directory.output_path == "") {
    //     sigmas_output = "config_out.json";
    // } else {
    //     sigmas_output = run_directory.output_path;
    // }

    // cout << "Opening output file: " << sigmas_output
    //      << " to dump fit results in" << endl;
    // std::ofstream config_file_out(sigmas_output);
    // if (config_file_out.is_open()) {
    //     // The '.dump(4)' method adds a 4-space indentation for pretty
    //     // formatting
    //     config_file_out << config.dump(4) << std::endl;
    //     config_file_out.close();
    //     cout << "Successfully updated config.json with sigma_x ("
    //          << fit_results.sigma_x << ") and sigma_y (" <<
    //          fit_results.sigma_y
    //          << ")." << endl;
    // } else {
    //     cerr << "ERROR: Could not open config.json for writing!" << endl;
    //     cerr << "System_error: " << std::strerror(errno) << endl;
    // }

    // TFile* output_file = TFile::Open(output_path, "RECREATE");
    // if (!output_file || output_file->IsZombie()){
    // 	cerr << "ERROR: Did not open output file" << endl;
    // 	return -1;
    // }
    //
    // output_file->cd();
    // TTree* output_tree = new TTree("T5_Events", "Reconstructed T5 events");
    // // --- Single-value branches (per event) ---
    // int b_n_particles = 0;
    // int b_event_nr = 0;
    // bool b_HasValidHit;
    // bool b_HasMultipleScintillatorsHit;
    // bool b_HasOutOfTimeWindow;
    // bool b_HasInTimeWindow;
    //
    // output_tree->Branch("event_nr", &b_event_nr, "event_nr/I");
    // output_tree->Branch("T5_particle_nr", &b_n_particles,
    // "T5_particle_nr/I"); output_tree->Branch("T5_HasValidHit",
    // &b_HasValidHit, "T5_HasValidHit/O");
    // output_tree->Branch("T5_HasMultipleScintillatorsHit",
    // &b_HasMultipleScintillatorsHit, "T5_HasMultipleScintillatorsHit/O");
    // output_tree->Branch("T5_HasOutOfTimeWindow", &b_HasOutOfTimeWindow,
    // "T5_HasOutOfTimeWindow/O"); output_tree->Branch("T5_HasInTimeWindow",
    // &b_HasInTimeWindow, "T5_HasInTimeWindow/O");
    //
    // // --- Vector branches (multiple hits per event) ---
    // // Primary hits -- hits in the expected timeframe
    // std::vector<int>* b_hit_is_in_bounds = new std::vector<int>();
    // std::vector<double>* b_hit_pos_x = new std::vector<double>();
    // std::vector<double>* b_hit_pos_y = new std::vector<double>();
    // std::vector<double>* b_hit_time = new std::vector<double>();
    //
    // output_tree->Branch("T5_hit_is_in_bounds", &b_hit_is_in_bounds);
    // output_tree->Branch("T5_hit_pos_x", &b_hit_pos_x);
    // output_tree->Branch("T5_hit_pos_y", &b_hit_pos_y);
    // output_tree->Branch("T5_hit_time", &b_hit_time);
    //
    // // Secondary hits -- hits outside of the main bunch
    // std::vector<bool>* b_secondary_hit_is_in_bounds = new
    // std::vector<bool>(); std::vector<double>* b_secondary_hit_pos_x = new
    // std::vector<double>(); std::vector<double>* b_secondary_hit_pos_y = new
    // std::vector<double>(); std::vector<double>* b_secondary_hit_time = new
    // std::vector<double>();
    //
    // output_tree->Branch("T5_secondary_hit_is_in_bounds",
    // &b_secondary_hit_is_in_bounds);
    // output_tree->Branch("T5_secondary_hit_pos_x", &b_secondary_hit_pos_x);
    // output_tree->Branch("T5_secondary_hit_pos_y", &b_secondary_hit_pos_y);
    // output_tree->Branch("T5_secondary_hit_time", &b_secondary_hit_time);
    //
    //
    // for (const auto& event : all_T5_hits){
    // 	b_n_particles = 0;
    // 	b_event_nr = event.event_nr;
    // 	b_HasValidHit = false;
    // 	b_HasMultipleScintillatorsHit = false;
    // 	b_HasOutOfTimeWindow = false;
    // 	b_HasInTimeWindow = false;
    //
    // 	b_hit_is_in_bounds->clear();
    // 	b_hit_pos_x->clear();
    // 	b_hit_pos_y->clear();
    // 	b_hit_time->clear();
    //
    // 	b_secondary_hit_is_in_bounds->clear();
    // 	b_secondary_hit_pos_x->clear();
    // 	b_secondary_hit_pos_y->clear();
    // 	b_secondary_hit_time->clear();
    //
    // 	if (!event.HasValidHit){
    // 		output_tree->Fill();
    // 		continue;
    // 	}
    // 	b_HasValidHit = event.HasValidHit;
    // 	b_HasMultipleScintillatorsHit = event.HasMultipleScintillatorsHit;
    // 	b_HasInTimeWindow = event.HasInTimeWindow;
    // 	b_HasOutOfTimeWindow = event.HasOutOfTimeWindow;
    //
    // 	for (const auto& hit : event.T5_hits){
    // 		if (hit.quality == HitQuality::AccidentalCoincidence) continue;
    // 		if (hit.is_in_time_window){
    // 			b_hit_time->push_back(hit.hit_time);
    // 			b_hit_pos_x->push_back(hit.position_x);
    // 			b_hit_pos_y->push_back(hit.position_y);
    // 			if (hit.quality == HitQuality::Perfect)
    // b_hit_is_in_bounds->push_back(true); 			else
    // b_hit_is_in_bounds->push_back(false);
    // 		}
    // 		else{
    // 			b_secondary_hit_time->push_back(hit.hit_time);
    // 			b_secondary_hit_pos_x->push_back(hit.position_x);
    // 			b_secondary_hit_pos_y->push_back(hit.position_y);
    // 			if (hit.quality == HitQuality::Perfect)
    // b_secondary_hit_is_in_bounds->push_back(true); 			else
    // b_secondary_hit_is_in_bounds->push_back(false);
    //
    // 		}
    // 		b_n_particles++;
    // 	}
    //
    // 	output_tree->Fill();
    // }
    // output_tree->Write();
    //
    // delete b_hit_time;
    // delete b_hit_pos_x;
    // delete b_hit_pos_y;
    // delete b_hit_is_in_bounds;
    // delete b_secondary_hit_time;
    // delete b_secondary_hit_pos_x;
    // delete b_secondary_hit_pos_y;
    // delete b_secondary_hit_is_in_bounds;
    //
    // output_file->Close();

    // std::ofstream file_out;
    // file_out.open("Beam_profile_widths.dat", std::ofstream::app);
    // file_out << run_number << "\t"
    // 	<< BEAM_MOMENTUM << "\t"
    // 	<< sigma_x << "\t"
    // 	<< sigma_y << "\t" << endl;

    //	app.Run();

    return 0;
}
