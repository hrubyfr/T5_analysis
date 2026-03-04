#include "buffer.h"
#include "return_TOF_position.h"
#include <string>

std::vector<std::string> special_hists = {};

void setup_histograms(Histograms &hists, TOF_reconstructor& recon){
	int n_scints = 8;
	int n_SiPMs = 16;
	hists.book1D("invalid_hit_times", "Average times in invalid hits;time - trigger[ns];count");
	hists.book1D("multiple_hit_times", "Average times of multiple hits;time - trigger[ns];count");
	hists.book1D("n_event_hits", "Number of hits in a single event", 9, 1, 10);
	for (int i = 0; i < n_SiPMs; i++){
		hists.book1D(Form("T5_number_of_hits_%i", i), 
			     Form("Number of hits in SiPM %i;count", i),
			     9, 1, 10);
	}
	hists.book2D("positions", "Reconstructed T5 positions;x[mm];y[mm];count",
			50,
			-recon. GetScintDimensionX(3)/2,
			recon. GetScintDimensionX(3)/2,
			8,
			recon. GetScintPositionY(7) - 16.25/2.0,
			recon. GetScintPositionY(0) + 16.25/2.0);
	for(int i = 0; i < n_scints; i++){
		hists. book1D(Form("positions_%i", i), Form("Reconstructed positions in scintillator %i;x[mm];count", i),
			      50, -recon. GetScintDimensionX(i)/2, recon. GetScintDimensionX(i)/2);
	}
}
