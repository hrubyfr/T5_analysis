#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
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

namespace {

void write_mask_summary_csv(const std::string &summary_path, int run_number,
                            const std::string &input_file,
                            long long total_events, long long good_events,
                            long long suspicious_events,
                            long long unclean_events, long long empty_events,
                            long long OutOfBoundEvents) {
    const bool file_exists = std::ifstream(summary_path).good();
    std::ofstream summary_file(summary_path);
    if (!summary_file.is_open()) {
        std::cerr << "ERROR: Could not open summary CSV: " << summary_path
                  << std::endl;
        return;
    }

    if (!file_exists) {
        summary_file
            << "run_number,input_file,total_events,good_events,good_fraction_"
               "pct,"
               "events_with_multiple_main_window_hits_on_same_scintillator,"
               "unclean_events(>1 hit in main time window),empty_events(no hit "
               "in main time window),main_events_OutOfBounds\n";
    }

    const double good_fraction_pct =
        total_events > 0 ? (100.0 * good_events) / total_events : 0.0;

    summary_file << run_number << ',' << input_file << ',' << total_events
                 << ',' << good_events << ',' << std::fixed
                 << std::setprecision(3) << good_fraction_pct << ','
                 << suspicious_events << ',' << unclean_events << ','
                 << empty_events << ',' << OutOfBoundEvents << '\n';
    summary_file.close();
}

} // namespace

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
                 << " -r <run_number> [-i <input_file>] [-o "
                    "<output_file>] [-d]"
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

        if (!current_output_path.EndsWith("/")) {
            current_output_path += "/";
        }
        current_output_path = output_path + base_name;

        auto file = TFile::Open(filename, "READ");

        if (!file || file->IsZombie()) {
            cerr << "ERROR, file did not open: " << filename << endl;
            continue;
        }

        // std::ifstream
        // config_file("/home/frantisek/scripts/config.json");
        // nlohmann::json config = nlohmann::json:: parse(config_file);
        // auto& run_config = config[std::to_string(run_number)];
        // BEAM_MOMENTUM = run_config.value("Beam momentum (MeV/c)", 0);
        // cout << "Config file loaded, beam momentum is " <<
        // BEAM_MOMENTUM << "MeV/c" << endl ;

        bool hardware_processed_data = false;
        string tree_name = "";
        if (file->GetListOfKeys()->Contains("ProcessedWaveforms")) {
            tree_name = "ProcessedWaveforms";
            hardware_processed_data = true;
        } else if (file->GetListOfKeys()->Contains("WCTEReadoutWindows")) {
            tree_name = "WCTEReadoutWindows";
        } else {
            cerr << "ERROR: Input file does not contain recognized "
                    "tree "
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
        std::unique_ptr<TTreeReaderValue<vector<float>>> arr_pmt_charges;
        std::unique_ptr<TTreeReaderValue<vector<int>>> arr_pmt_ids;
        std::unique_ptr<TTreeReaderValue<vector<int>>> arr_mpmt_ids;

        // needed for array input of the processed waveforms file

        const int MAX_HITS = 4000;

        std::unique_ptr<TTreeReaderValue<int>> n_pmt_times;
        std::unique_ptr<TTreeReaderValue<int>> n_pmt_charges;
        std::unique_ptr<TTreeReaderValue<int>> n_pmt_chans;
        std::unique_ptr<TTreeReaderValue<int>> n_pmt_cards;

        std::unique_ptr<TTreeReaderArray<double>> processed_hit_times;
        std::unique_ptr<TTreeReaderArray<float>> processed_hit_charges;
        std::unique_ptr<TTreeReaderArray<int>> processed_hit_cards;
        std::unique_ptr<TTreeReaderArray<int>> processed_hit_chans;

        if (hardware_processed_data) {
            n_pmt_times =
                std::make_unique<TTreeReaderValue<int>>(tree, "nhit_time");
            n_pmt_charges =
                std::make_unique<TTreeReaderValue<int>>(tree, "nhit_charge");
            n_pmt_chans =
                std::make_unique<TTreeReaderValue<int>>(tree, "nhit_chan");
            n_pmt_cards =
                std::make_unique<TTreeReaderValue<int>>(tree, "nhit_card");

            processed_hit_times =
                std::make_unique<TTreeReaderArray<double>>(tree, "hit_time");
            processed_hit_charges =
                std::make_unique<TTreeReaderArray<float>>(tree, "hit_charge");
            processed_hit_cards =
                std::make_unique<TTreeReaderArray<int>>(tree, "hit_card");
            processed_hit_chans =
                std::make_unique<TTreeReaderArray<int>>(tree, "hit_chan");

        } else {
            arr_pmt_times = std::make_unique<TTreeReaderValue<vector<double>>>(
                tree, "hit_pmt_times");
            arr_pmt_charges = std::make_unique<TTreeReaderValue<vector<float>>>(
                tree, "hit_pmt_charges");
            arr_pmt_ids = std::make_unique<TTreeReaderValue<vector<int>>>(
                tree, "hit_pmt_channel_ids");
            arr_mpmt_ids = std::make_unique<TTreeReaderValue<vector<int>>>(
                tree, "hit_mpmt_card_ids");
        }

        // tree-> SetBranchAddress("beamline_pmt_qdc_ids",
        // &arr_bm_charge_ids); tree->
        // SetBranchAddress("beamline_pmt_tdc_ids", &arr_bm_time_ids);
        // tree-> SetBranchAddress("beamline_pmt_tdc_times",
        // &arr_bm_times); tree->
        // SetBranchAddress("beamline_pmt_qdc_charges",
        // &arr_bm_charges);

        Cuts cut;
        TOF_reconstructor recon;
        // Histograms hists;
        // setup_histograms(hists, recon);
        int n_pass_cut = 0;
        int n_T5_valid_events = 0;
        long long total_events = tree.GetEntries();
        long long good_events = 0;
        long long suspicious_main_window_events = 0;
        auto n_events = total_events;
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

        int n_clean_events = 0;
        int n_empty_events = 0;
        int n_unclean_events = 0;

        vector<event_T5_detection> all_T5_hits;

        while (tree.Next()) {
            auto i = tree.GetCurrentEntry();

            if (i >= n_events)
                break;

            vector<double> vec_hit_time;
            vector<double> vec_hit_charge;
            vector<int> vec_hit_chan;
            vector<int> vec_hit_card;

            event_T5_detection detections;

            if (hardware_processed_data) {
                if (**n_pmt_times >= MAX_HITS) {
                    // continue and push back dummy
                    // detections
                    std::cerr << "ERROR: Too many hits in event "
                              << tree.GetCurrentEntry() << " nhits "
                              << **n_pmt_times << " (max " << MAX_HITS << ")"
                              << std::endl;
                    detections.event_nr = tree.GetCurrentEntry();
                    all_T5_hits.push_back(detections);
                    continue;
                }
                // change array output into std vectors
                for (int j = 0; j < **n_pmt_times; j++) {
                    if ((*processed_hit_cards)[j] == T5_CONFIG::T5_MPMT_ID) {
                        vec_hit_time.push_back((*processed_hit_times)[j]);
                        vec_hit_charge.push_back((*processed_hit_charges)[j]);
                        vec_hit_card.push_back((*processed_hit_cards)[j]);
                        vec_hit_chan.push_back((*processed_hit_chans)[j]);
                    }
                }
            } else {
                auto &cards = (**arr_mpmt_ids);
                for (int j = 0; j < cards.size(); j++) {
                    if (cards[j] == T5_CONFIG::T5_MPMT_ID) {
                        vec_hit_time.push_back((**arr_pmt_times)[j]);
                        vec_hit_charge.push_back((**arr_pmt_charges)[j]);
                        vec_hit_card.push_back((**arr_mpmt_ids)[j]);
                        vec_hit_chan.push_back((**arr_pmt_ids)[j]);
                    }
                }
            }

            if (i % verb == 0)
                cout << "\rAnalyzed " << i << " of " << n_events
                     << std::setprecision(2) << std::fixed << " events ("
                     << static_cast<float>(i) / n_events * 100 << " %)"
                     << std::flush << endl;

            // RVecI bm_time_ids(arr_bm_time_ids->data(),
            // arr_bm_time_ids->size()); RVecI
            // bm_charge_ids(arr_bm_charge_ids->data(),
            // arr_bm_charge_ids->size()); RVecF
            // bm_times(arr_bm_times->data(), arr_bm_times->size());
            // RVecF bm_charges(arr_bm_charges->data(),
            // arr_bm_charges->size());
            if (vec_hit_card.empty()) {
                detections.event_nr = i;
                all_T5_hits.push_back(detections);
                continue;
            }
            RVecI pmt_ids(vec_hit_chan.data(), vec_hit_chan.size());
            RVecI mpmt_ids(vec_hit_card.data(), vec_hit_card.size());

            if (!cut.hit_T5(mpmt_ids, pmt_ids)) {
                detections.event_nr = i;
                all_T5_hits.push_back(detections);
                continue;
            }

            // n_pass_cut++;

            // auto mask_T5_board = (mpmt_ids == cut.get_T5_board());
            // auto T5_board_ids = pmt_ids[mask_T5_board];
            // auto T5_board_times = pmt_times[mask_T5_board];

            detections = recon.Return_position(i, vec_hit_card, vec_hit_chan,
                                               vec_hit_time, vec_hit_charge);

            std::array<int, 8> main_window_hits_per_scintillator = {};
            int n_main_window_hits = 0;
            for (const auto &hit : detections.T5_hits) {
                if (hit.is_in_time_window)
                    n_main_window_hits++;
                if (hit.is_in_time_window && hit.is_valid_hit &&
                    hit.scintillator_id >= 0 &&
                    hit.scintillator_id <
                        static_cast<int>(
                            main_window_hits_per_scintillator.size())) {
                    ++main_window_hits_per_scintillator[hit.scintillator_id];
                }
            }
            const bool has_multiple_hits_on_same_scintillator =
                std::any_of(main_window_hits_per_scintillator.begin(),
                            main_window_hits_per_scintillator.end(),
                            [](int count) { return count > 1; });
            if (has_multiple_hits_on_same_scintillator) {
                ++suspicious_main_window_events;
                std::cerr << "WARNING: event " << i
                          << " has multiple main-window hits on the same "
                             "scintillator; review the reconstruction output."
                          << std::endl;
            }
            // cout << "Event " << i << " Diagnostics: " << endl;
            // cout << "\tnumber of main window hits: " << n_main_window_hits
            //      << endl;
            if (detections.IsClean) {
                n_clean_events++;
                // cout << "\tevent is clean" << endl;
            }
            // if (detections.HasInTimeWindow)
            //     cout << "\tevent has a main time window hit" << endl;
            // else if (!detections.HasInTimeWindow)
            //     cout << "\tevent does not have a main time window hit" <<
            //     endl;

            if (detections.IsClean && detections.HasInTimeWindow) {
                ++good_events;
            }
            if (n_main_window_hits == 0)
                n_empty_events++;
            if (n_main_window_hits > 1)
                n_unclean_events++;

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

            all_T5_hits.push_back(detections);
        }

        cout << "Run " << RUN_NUMBER << " Diagnostics: " << endl;
        cout << "\t" << n_events << " total events" << endl;
        cout << "\t" << n_clean_events << " events were clean" << endl;
        cout << "\t" << n_unclean_events
             << " events were unclean (more than 1 hit in main time window)"
             << endl;
        cout << "\t" << n_empty_events
             << " events were empty (no hit in main time window)" << endl;

        file->Close();

        // ROOT output-file writing is disabled for file-by-file Python-driven
        // analysis runs that only need the summary statistics.
        TFile *output_file = TFile::Open(current_output_path, "RECREATE");
        if (!output_file || output_file->IsZombie()) {
            cerr << "ERROR: Did not open output file: " << current_output_path
                 << endl;
            continue;
        }

        // output_file->cd();
        // TH1D *h_main_window_hits_per_scintillator =
        //     new TH1D("T5_main_window_hits_per_scintillator",
        //              "Main-window hits per scintillator;Scintillator
        //              ID;Count", 8, -0.5, 7.5);
        // TH1D *h_suspicious_events =
        //     new TH1D("T5_suspicious_multi_hits_same_scintillator",
        //              "Events with >1 main-window hit on the same
        //              scintillator;" "Event count;Entries", 2, -0.5, 1.5);
        // h_suspicious_events->SetBinContent(2, suspicious_main_window_events);
        // h_suspicious_events->SetBinContent(
        //     1, total_events - suspicious_main_window_events);

        // for (const auto &event : all_T5_hits) {
        //     for (const auto &hit : event.T5_hits) {
        //         if (hit.is_in_time_window && hit.is_valid_hit &&
        //             hit.scintillator_id >= 0 && hit.scintillator_id < 8) {
        //             h_main_window_hits_per_scintillator->Fill(
        //                 hit.scintillator_id);
        //         }
        //     }
        // }

        // h_main_window_hits_per_scintillator->Write();
        // h_suspicious_events->Write();

        TTree *output_tree = new TTree("T5_Events", "Reconstructed T5 events");
        // --- Single-value branches (per event) ---
        int b_n_main_particles = 0;
        int b_event_nr = 0;
        double b_main_hit_time = -9999;
        double b_main_hit_charge = -9999;
        double b_main_hit_sipm_time_left = -9999;
        double b_main_hit_sipm_time_right = -9999;
        double b_main_hit_sipm_charge_left = -9999;
        double b_main_hit_sipm_charge_right = -9999;
        double b_main_position_x = -9999;
        double b_main_position_y = -9999;
        double b_main_position_x_error = -9999;
        double b_main_position_y_error = -9999;
        // int b_main_hit_scint_ID = -9999;
        int b_T5_hit_bitmask;

        // Diagnostics

        int n_main_OutOfBounds = 0;

        output_tree->Branch("event_nr", &b_event_nr, "event_nr/I");
        output_tree->Branch("t5_hit_bitmask", &b_T5_hit_bitmask,
                            "t5_hit_bitmask/I");
        output_tree->Branch("t5_n_main_bunch_particles", &b_n_main_particles,
                            "t5_n_main_bunch_particles/I");
        output_tree->Branch("t5_main_hit_time", &b_main_hit_time,
                            "t5_main_hit_time/D");
        output_tree->Branch("t5_main_hit_charge", &b_main_hit_charge,
                            "t5_main_hit_charge/D");
        output_tree->Branch("t5_main_hit_sipm_time_left",
                    &b_main_hit_sipm_time_left,
                    "t5_main_hit_sipm_time_left/D");
        output_tree->Branch("t5_main_hit_sipm_time_right",
                    &b_main_hit_sipm_time_right,
                    "t5_main_hit_sipm_time_right/D");
        output_tree->Branch("t5_main_hit_sipm_charge_left",
                    &b_main_hit_sipm_charge_left,
                    "t5_main_hit_sipm_charge_left/D");
        output_tree->Branch("t5_main_hit_sipm_charge_right",
                    &b_main_hit_sipm_charge_right,
                    "t5_main_hit_sipm_charge_right/D");
        output_tree->Branch("t5_main_hit_pos_x", &b_main_position_x,
                            "t5_main_hit_pos_x/D");
        output_tree->Branch("t5_main_hit_pos_y", &b_main_position_y,
                            "t5_main_hit_pos_y/D");
        output_tree->Branch("t5_main_hit_pos_x_error", &b_main_position_x_error,
                            "t5_main_hit_pos_x_error/D");
        output_tree->Branch("t5_main_hit_pos_y_error", &b_main_position_y_error,
                            "t5_main_hit_pos_y_error/D");
        // output_tree->Branch("t5_main_hit_scint_id", &b_main_hit_scint_ID,
        //                     "t5_main_hit_scint_id/I");

        // --- Vector branches (multiple hits per event) ---
        // All hits -- Any hit other and including the main hit
        std::vector<double> *b_all_hits_pos_x = new std::vector<double>();
        std::vector<double> *b_all_hits_pos_y = new std::vector<double>();
        std::vector<double> *b_all_hits_pos_x_error = new std::vector<double>();
        std::vector<double> *b_all_hits_pos_y_error = new std::vector<double>();
        std::vector<double> *b_all_hits_time = new std::vector<double>();
        std::vector<double> *b_all_hits_charge = new std::vector<double>();
        std::vector<double> *b_all_hits_sipm_time_left =
            new std::vector<double>();
        std::vector<double> *b_all_hits_sipm_time_right =
            new std::vector<double>();
        std::vector<double> *b_all_hits_sipm_charge_left =
            new std::vector<double>();
        std::vector<double> *b_all_hits_sipm_charge_right =
            new std::vector<double>();
        // std::vector<int> *b_all_hits_scintillator_id = new
        // std::vector<int>();
        vector<int> *b_all_hits_is_in_time_window = new std::vector<int>();
        vector<int> *b_all_hits_is_in_bounds = new std::vector<int>();

        output_tree->Branch("t5_all_hits_pos_x", &b_all_hits_pos_x);
        output_tree->Branch("t5_all_hits_pos_y", &b_all_hits_pos_y);
        output_tree->Branch("t5_all_hits_pos_x_error", &b_all_hits_pos_x_error);
        output_tree->Branch("t5_all_hits_pos_y_error", &b_all_hits_pos_y_error);
        output_tree->Branch("t5_all_hits_time", &b_all_hits_time);
        output_tree->Branch("t5_all_hits_charge", &b_all_hits_charge);
        output_tree->Branch("t5_all_hits_sipm_time_left",
                    &b_all_hits_sipm_time_left);
        output_tree->Branch("t5_all_hits_sipm_time_right",
                    &b_all_hits_sipm_time_right);
        output_tree->Branch("t5_all_hits_sipm_charge_left",
                    &b_all_hits_sipm_charge_left);
        output_tree->Branch("t5_all_hits_sipm_charge_right",
                    &b_all_hits_sipm_charge_right);
        // output_tree->Branch("t5_all_hits_scintillator_id",
        //                     &b_all_hits_scintillator_id);
        output_tree->Branch("t5_all_hits_is_in_time_window",
                            &b_all_hits_is_in_time_window);
        output_tree->Branch("t5_all_hits_is_in_bounds",
                            &b_all_hits_is_in_bounds);

        // y error calculation
        auto sigma_y_error = T5_CONFIG::SCINT_BLOCK_HEIGHT / sqrt(12);

        int n_hits_good = 0, n_hits_accidental = 0;
        int n_multiple_main_hit_good_main_accidental_secondary_events = 0;

        int n_events_bitmask_0 = 0, n_events_bitmask_1 = 0,
            n_events_bitmask_2 = 0, n_events_bitmask_4 = 0, n_events_others = 0;

        for (const auto &event : all_T5_hits) {
            b_event_nr = event.event_nr;
            b_T5_hit_bitmask = 0;

            b_n_main_particles = 0;

            b_main_hit_charge = -9999;
            b_main_hit_time = -9999;
            b_main_hit_sipm_time_left = -9999;
            b_main_hit_sipm_time_right = -9999;
            b_main_hit_sipm_charge_left = -9999;
            b_main_hit_sipm_charge_right = -9999;
            b_main_position_x = -9999;
            b_main_position_y = -9999;
            b_main_position_x_error = -9999;
            b_main_position_y_error = -9999;
            // b_main_hit_scint_ID = -1;

            b_all_hits_time->clear();
            b_all_hits_charge->clear();
            b_all_hits_sipm_time_left->clear();
            b_all_hits_sipm_time_right->clear();
            b_all_hits_sipm_charge_left->clear();
            b_all_hits_sipm_charge_right->clear();
            b_all_hits_pos_x->clear();
            b_all_hits_pos_y->clear();
            b_all_hits_pos_x_error->clear();
            b_all_hits_pos_y_error->clear();
            // b_all_hits_scintillator_id->clear();
            b_all_hits_is_in_time_window->clear();
            b_all_hits_is_in_bounds->clear();

            bool had_accidental_main_hit_candidate = false;
            const T5_hit *main_hit_candidate = nullptr;
            for (const auto &hit : event.T5_hits) {
                b_all_hits_time->push_back(hit.raw_time);
                b_all_hits_charge->push_back(hit.hit_charge);
                b_all_hits_sipm_time_left->push_back(hit.sipm_time_b);
                b_all_hits_sipm_time_right->push_back(hit.sipm_time_a);
                b_all_hits_sipm_charge_left->push_back(hit.sipm_charge_b);
                b_all_hits_sipm_charge_right->push_back(hit.sipm_charge_a);
                b_all_hits_pos_x->push_back(hit.position_x);
                b_all_hits_pos_y->push_back(hit.position_y);
                b_all_hits_pos_x_error->push_back(hit.uncertainty);
                b_all_hits_pos_y_error->push_back(sigma_y_error);
                if (hit.quality == HitQuality::AccidentalCoincidence) {
                    n_hits_accidental++;
                    b_all_hits_is_in_bounds->push_back(0);
                } else {
                    b_all_hits_is_in_bounds->push_back(1);
                    n_hits_good++;
                }
                if (hit.is_in_time_window)
                    b_all_hits_is_in_time_window->push_back(1);
                else
                    b_all_hits_is_in_time_window->push_back(0);

                if (!hit.is_in_time_window)
                    continue;

                b_n_main_particles++;
                bool is_accidental =
                    hit.quality == HitQuality::AccidentalCoincidence;
                if (is_accidental)
                    had_accidental_main_hit_candidate = true;
                bool candidate_is_accidental =
                    main_hit_candidate && (main_hit_candidate->quality ==
                                           HitQuality::AccidentalCoincidence);
                if (!main_hit_candidate ||
                    (!is_accidental && candidate_is_accidental) ||
                    (!is_accidental && !candidate_is_accidental &&
                     hit.raw_time < main_hit_candidate->raw_time))
                    main_hit_candidate = &hit;
            }
            if (b_n_main_particles == 0)
                b_T5_hit_bitmask |= (1 << 0);
            if (main_hit_candidate) {
                bool is_accidental = main_hit_candidate->quality ==
                                     HitQuality::AccidentalCoincidence;

                if (b_n_main_particles > 1)
                    b_T5_hit_bitmask |= (1 << 1);
                else if (is_accidental) {
                    b_T5_hit_bitmask |= (1 << 2);
                } else if (b_T5_hit_bitmask == 0) {

                    b_main_hit_time = main_hit_candidate->raw_time;
                    b_main_hit_charge = main_hit_candidate->hit_charge;
                    b_main_hit_sipm_time_left =
                        main_hit_candidate->sipm_time_b;
                    b_main_hit_sipm_time_right =
                        main_hit_candidate->sipm_time_a;
                    b_main_hit_sipm_charge_left =
                        main_hit_candidate->sipm_charge_b;
                    b_main_hit_sipm_charge_right =
                        main_hit_candidate->sipm_charge_a;
                    b_main_position_x = main_hit_candidate->position_x;
                    b_main_position_y = main_hit_candidate->position_y;
                    b_main_position_x_error = main_hit_candidate->uncertainty;
                    b_main_position_y_error = sigma_y_error;
                }
            }

            if (b_T5_hit_bitmask == 2 && had_accidental_main_hit_candidate)
                n_multiple_main_hit_good_main_accidental_secondary_events++;
            if (b_T5_hit_bitmask == 0)
                n_events_bitmask_0++;
            else if (b_T5_hit_bitmask == 1)
                n_events_bitmask_1++;
            else if (b_T5_hit_bitmask == 2)
                n_events_bitmask_2++;
            else if (b_T5_hit_bitmask == 4)
                n_events_bitmask_4++;
            else
                n_events_others++;
            output_tree->Fill();
        }
        cout << "Output diagnostics: " << endl;
        cout << "\tnumber of good hits: " << n_hits_good << endl;
        cout << "\tnumber of accidental hits: " << n_hits_accidental << endl;
        cout << "\tnumber of events with a good main hit candidate, but an "
                "accidental hit in the main time window: "
             << n_multiple_main_hit_good_main_accidental_secondary_events
             << endl
             << "\tnumber of events with bitmask 0: " << n_events_bitmask_0
             << endl
             << "\tnumber of events with bitmask 1: " << n_events_bitmask_1
             << endl
             << "\tnumber of events with bitmask 2: " << n_events_bitmask_2
             << endl
             << "\tnumber of events with bitmask 4: " << n_events_bitmask_4
             << endl
             << "\tnumber of events with any other bitmask: " << n_events_others
             << endl;

        output_tree->Write();

        delete b_all_hits_charge;
        delete b_all_hits_time;
        delete b_all_hits_sipm_time_left;
        delete b_all_hits_sipm_time_right;
        delete b_all_hits_sipm_charge_left;
        delete b_all_hits_sipm_charge_right;
        delete b_all_hits_pos_x;
        delete b_all_hits_pos_y;
        delete b_all_hits_pos_x_error;
        delete b_all_hits_pos_y_error;
        // delete b_all_hits_scintillator_id;
        delete b_all_hits_is_in_bounds;
        delete b_all_hits_is_in_time_window;

        // write_mask_summary_csv(
        //     (std::string(current_output_path.Data()) + "run.summary.csv"),
        //     run_number, filename.Data(), total_events, good_events,
        //     suspicious_main_window_events, n_unclean_events, n_empty_events,
        //     n_main_OutOfBounds);

        output_file->Close();
    } // end loop over input files

    return 0;
}
