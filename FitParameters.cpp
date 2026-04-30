#include "FitParameters.h"
#include <cmath>
#include <iostream>
#include <nlohmann/json.hpp>
#include <ostream>

FitParameters::FitParameters(std::string fname) : f_config_file_name(fname) {
    Load_config();
}

void FitParameters::Convert_to_width(SiPM_resolution &res) {
    res.width_cm = f_v_eff_cmns * res.time_ns / 2;
    res.uncertainty_cm =
        res.width_cm * sqrt(pow(res.uncertainty_ns / res.time_ns, 2) +
                            pow(f_v_eff_uncertainty_cmns / f_v_eff_cmns, 2));
}

void FitParameters::Load_config() {
    std::ifstream file(f_config_file_name);

    if (!file || !file.is_open()) {
        std::cerr << "ERROR: config file loading SiPM resolutions did not open"
                  << std::endl;
        return;
    }

    nlohmann::json j = nlohmann::json::parse(file);

    f_v_eff_cmns = j.at("v_eff").get<double>() / 10;
    f_v_eff_uncertainty_cmns = j.at("v_eff_error").get<double>() / 10;

    for (const auto &item : j.items()) {
        auto key = item.key();
        if (key.find("sigma_sipm_") == 0 &&
            key.find("_uncertainty") == std::string::npos) {
            SiPM_resolution res;
            res.sipm_nr = std::stoi(item.key().substr(11));
            res.time_ns = item.value().get<double>();
            res.uncertainty_ns = j.at(key + "_uncertainty").get<double>();
            Convert_to_width(res);
            std::cout << "Resolution converted, pushing back sigma_sipm_cm = "
                      << res.width_cm << " (calculated with effective speed "
                      << f_v_eff_cmns << ")" << std::endl;
            f_resolutions.push_back(res);
        }
    }
    file.close();
};

std::vector<SiPM_resolution> FitParameters::GetResolutions() const {
    return f_resolutions;
}
