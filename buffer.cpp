#include "buffer.h"
#include "TString.h"
#include "return_TOF_position.h"
#include <string>

std::vector<std::string> special_hists = {};

void setup_histograms(Histograms &hists, TOF_reconstructor &recon) {
    int n_scints = 8;
    int n_SiPMs = 16;

    hists.book1D("T5_to_PMT_times", "Time of T5 to tank PMTs", 500, -5000,
                 5000);

    hists.book1D("trigger_times", "Time of the trigger (channel 19)", 250,
                 -2500, 2500);

    hists.book1D("valid_hit_times",
                 "Average time of valid hits;time - trigger[ns];count", 300,
                 -250, 250);
    hists.book1D("invalid_hit_times",
                 "Average times in invalid hits;time - trigger[ns];count", 500,
                 -5000, 5000);
    hists.book1D("multiple_hit_times",
                 "Average times of multiple hits;time - trigger[ns];count", 200,
                 -250, 250);
    hists.book1D("n_event_hits", "Number of hits in a single event", 9, 1, 10);
    for (int i = 0; i < n_scints; i++) {
        hists.book2D(Form("hit_charges_2D_%i", i),
                     Form("Charges in scintillator %i;Q_{A};Q_{B};count", i),
                     250, 0, 16000, 250, 0, 16000);
        hists.book1D(Form("hit_charges_%i", i),
                     Form("Total charge in scintillator %i;Q_{total};count", i),
                     500, 0, 16000);
        hists.book1D(Form("hit_raw_times_%i", i),
                     Form("Time in scintillator %i -- uncorrected by trigger "
                          "time;raw_time [ns];count",
                          i),
                     500, -5000, 5000);
    }
    for (int i = 0; i < n_SiPMs; i++) {
        hists.book1D(Form("T5_number_of_hits_%i", i),
                     Form("Number of hits in SiPM %i;count", i), 9, 1, 10);
    }
    for (int time_low = -65; time_low < 36; time_low += 10) {
        int iter = (time_low + 65) / 10;
        hists.book2D(
            Form("cutup_positions_%i", iter),
            Form("Reconstructed T5 positions for time - trigger (%i to %i)",
                 time_low, time_low + 10),
            20, -recon.GetScintDimensionX(3), recon.GetScintDimensionX(3), 8,
            recon.GetScintPositionY(7) - 16.25 / 2.0,
            recon.GetScintPositionY(0) + 16.25 / 2.0);
    }
    hists.book2D("positions", "Reconstructed T5 positions;x[mm];y[mm];count",
                 30, -recon.GetScintDimensionX(3), recon.GetScintDimensionX(3),
                 8, recon.GetScintPositionY(7) - 16.25 / 2.0,
                 recon.GetScintPositionY(0) + 16.25 / 2.0);
    hists.book2D("positions_timecut",
                 "Reconstructed T5 positions with time veto;x[mm];y[mm];count",
                 30, -recon.GetScintDimensionX(3), recon.GetScintDimensionX(3),
                 8, recon.GetScintPositionY(7) - 16.25 / 2.0,
                 recon.GetScintPositionY(0) + 16.25 / 2.0);
    hists.book2D(
        "positions_chargecut",
        "Reconstructed T5 positions with charge veto;x[mm];y[mm];count", 30,
        -recon.GetScintDimensionX(3), recon.GetScintDimensionX(3), 8,
        recon.GetScintPositionY(7) - 16.25 / 2.0,
        recon.GetScintPositionY(0) + 16.25 / 2.0);
    hists.book2D("positions_timechargecut",
                 "Reconstructed T5 positions with charge and time "
                 "veto;x[mm];y[mm];count",
                 30, -recon.GetScintDimensionX(3), recon.GetScintDimensionX(3),
                 8, recon.GetScintPositionY(7) - 16.25 / 2.0,
                 recon.GetScintPositionY(0) + 16.25 / 2.0);
    for (int i = 0; i < n_scints; i++) {
        hists.book1D(
            Form("positions_%i", i),
            Form("Reconstructed positions in scintillator %i;x[mm];count", i),
            300, -recon.GetScintDimensionX(i), recon.GetScintDimensionX(i));
    }
}
