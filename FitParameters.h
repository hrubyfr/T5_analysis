#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct SiPM_resolution {
    int sipm_nr;
    double time_ns;
    double uncertainty_ns;
    double width_cm;
    double uncertainty_cm;
};

class FitParameters {
  public:
    explicit FitParameters(
        std::string fname =
            "/home/frantisek/Analysis/configs/T5_reconstruction_values.json");

    std::vector<SiPM_resolution> GetResolutions() const;
    double Get_veff() const { return f_v_eff_cmns; };

  private:
    void Load_config();
    void Convert_to_width(SiPM_resolution &);

    std::string f_config_file_name;
    std::vector<SiPM_resolution> f_resolutions;
    double f_v_eff_cmns;
    double f_v_eff_uncertainty_cmns;
};
