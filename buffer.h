#pragma once
#include "return_TOF_position.h"
#include "utils.h"

extern std::vector<std::string> special_hists;
void setup_histograms(Histograms& hists, TOF_reconstructor& recon);
