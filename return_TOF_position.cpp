#include "./return_TOF_position.h"
#include <cstdlib>
#include <ranges>
#include <set>
#include <vector>

using namespace std;
using namespace T5_CONFIG;

// TOF_reconstructor class constructor -- private values contain only fit values
// -- can change depending on
TOF_reconstructor::TOF_reconstructor(double v_eff)
    : _v_eff(v_eff), v_eff_uncertainty(33.7644),
      sigma_sipm_i({0.3623299596435143, 0.41314637562035883, 0.4214532281907809,
                    0.34389155240346503, 0.2769688564632169, 0.2944165197343617,
                    0.3388218044201263, 0.28871976899865825}),
      sigma_sipm_i_uncertainties({0.007049698438571145, 0.03205329644174973,
                                  0.04459054671484605, 0.06588460376135835,
                                  0.08178483906009885, 0.06380836194128096,
                                  0.039075506816850014, 0.008781947884922536}),

      _verbose(0) {}

void TOF_reconstructor::SetVeff(double v) { _v_eff = v; }
void TOF_reconstructor::SetVeffUncertainty(double uncertainty) {
    v_eff_uncertainty = uncertainty;
}
void TOF_reconstructor::SetVerbosity(int i) { _verbose = i; }

double TOF_reconstructor::GetVeff() const { return _v_eff; }
double TOF_reconstructor::GetVeff_uncertainty() const {
    return v_eff_uncertainty;
}
bool TOF_reconstructor::GetVerbosity() const { return _verbose; }

double TOF_reconstructor::GetScintDimensionX(int i) const {
    return SCINT_DIMENSIONS[i];
}
double TOF_reconstructor::GetScintPositionY(int i) const {
    return SCINT_Y_POSITIONS[i];
}
double TOF_reconstructor::Get_scint_xmax(int i) const {
    return SCINT_DIMENSIONS[i] / 2;
}
double TOF_reconstructor::Get_scint_xmin(int i) const {
    return -SCINT_DIMENSIONS[i] / 2;
}
double TOF_reconstructor::Get_ymax() const {
    double scintillator_block_halfheight = SCINT_BLOCK_HEIGHT / 2.0;
    return SCINT_Y_POSITIONS[0] + scintillator_block_halfheight;
}
event_T5_detection TOF_reconstructor::Return_position(
    const int event_nr, const vector<int> &hit_mpmt_ids,
    const vector<int> &hit_pmt_ids, const vector<double> &hit_pmt_times,
    const vector<double> &hit_pmt_charges) {

    event_T5_detection detection;
    detection.event_nr = event_nr;
    vector<T5_hit> all_hits;

    // Create 16 vectors where to store all the times detected by T5 SiPMs,
    // some will be empty, some will have multiple hits, ideally paired --
    // later will check for hits in paired detectors
    vector<vector<double>> T5_times(N_T5_SIPMS);
    vector<vector<double>> T5_charges(N_T5_SIPMS);
    vector<vector<double>> T5_raw_times(N_T5_SIPMS);
    bool is_valid_hit = false;
    double trigger_time = 0;
    for (size_t i = 0; i < hit_pmt_ids.size(); i++) {
        if (hit_mpmt_ids.at(i) == T5_MPMT_ID &&
            hit_pmt_ids.at(i) == T5_TRIGGER_ID) {
            // Take the first trigger event (should be only one,
            // hopefully), and set it as trigger time, then break
            // the for loop
            trigger_time = hit_pmt_times.at(i);
            break;
        }
    }
    for (size_t i = 0; i < hit_mpmt_ids.size(); i++) {
        auto mPMT_id = hit_mpmt_ids.at(i);
        if (mPMT_id != T5_MPMT_ID)
            continue;
        auto PMT_id = hit_pmt_ids.at(i);
        auto SiPM_index = GetSiPMIndex(PMT_id);
        if (SiPM_index == -1)
            continue;
        // correct the measured T5 sipm time by the trigger time
        T5_times.at(SiPM_index).push_back(hit_pmt_times.at(i) - trigger_time);
        T5_raw_times.at(SiPM_index).push_back(hit_pmt_times.at(i));
        T5_charges.at(SiPM_index).push_back(hit_pmt_charges.at(i));
    }
    // Loop over all saved times in the corresponding vectors, compare all
    // the times
    for (int i = 0; i < N_T5_SCINTS; i++) {
        if (T5_times[i].empty() || T5_times[i + 8].empty())
            continue;
        for (int isipm_a = 0; isipm_a < T5_times[i].size(); isipm_a++) {
            auto sipm_time_a = T5_times[i].at(isipm_a);
            auto sipm_raw_time_a = T5_raw_times[i].at(isipm_a);
            auto sipm_charge_a = T5_charges[i].at(isipm_a);
            bool time_a_valid = (sipm_time_a < EXPECTED_DETECTION_TIME_MAX &&
                                 sipm_time_a > EXPECTED_DETECTION_TIME_MIN);
            for (int isipm_b = 0; isipm_b < T5_times[i + 8].size(); isipm_b++) {
                auto sipm_time_b = T5_times[i + 8].at(isipm_b);
                auto sipm_charge_b = T5_charges[i + 8].at(isipm_b);
                auto sipm_raw_time_b = T5_raw_times[i + 8].at(isipm_b);
                bool time_b_valid =
                    (sipm_time_b < EXPECTED_DETECTION_TIME_MAX &&
                     sipm_time_b > EXPECTED_DETECTION_TIME_MIN);
                T5_hit hit;
                hit.sipm_time_a = sipm_time_a;
                hit.sipm_time_b = sipm_time_b;
                hit.sipm_charge_a = sipm_charge_a;
                hit.sipm_charge_b = sipm_charge_b;
                hit.raw_time = (sipm_raw_time_a + sipm_raw_time_b) / 2.0;
                hit.hit_charge = sipm_charge_a + sipm_charge_b;
                hit.is_in_time_window = time_a_valid && time_b_valid;
                if (!hit.is_in_time_window && _verbose) {
                    cout << "WARNING: SiPM time is out of "
                            "expected hit times"
                         << endl;
                }
                auto time_diff = sipm_time_a - sipm_time_b;
                double avg_time = (sipm_time_a + sipm_time_b) / 2;
                hit.hit_time = avg_time;
                is_valid_hit = true;
                double position_x =
                    (time_diff - SCINT_BIAS.at(i)) * _v_eff / 2.0;
                double position_y = SCINT_Y_POSITIONS[i];
                double uncertainty =
                    sqrt(pow(v_eff_uncertainty, 2) * pow(time_diff / 2, 2) +
                         pow(sigma_sipm_i[i], 2) * pow(_v_eff / 2, 2));
                double margin_of_error = 3 * uncertainty;

                if (abs(position_x) <= SCINT_DIMENSIONS[i] / 2) {
                    is_valid_hit = true;
                    hit.quality = HitQuality::Perfect;
                } else if (abs(position_x) <=
                           (SCINT_DIMENSIONS[i] / 2 + margin_of_error)) {
                    is_valid_hit = true;
                    hit.quality = HitQuality::OutOfBounds;
                    if (_verbose)
                        cout << "Hit was reconstructed "
                                "out of bounds -- due "
                                "to detector "
                                "smearing?"
                             << endl;
                } else {
                    is_valid_hit = false;
                    hit.quality = HitQuality::AccidentalCoincidence;
                    if (_verbose)
                        cout << "Hit was extremely out "
                                "of bounds, is not valid"
                             << endl;
                }
                // _verbose checks
                if (_verbose)
                    cout << "X coordinate is:\t" << position_x;

                // calculate the position uncertainty
                if (_verbose)
                    cout << "\tX uncertainty is:\t" << uncertainty << endl;
                hit.uncertainty = uncertainty;
                hit.position_x = position_x;
                hit.position_y = position_y;
                hit.is_valid_hit = is_valid_hit;
                hit.scintillator_id = i;
                // save the hit to the vector of all hits in the
                // event
                all_hits.push_back(hit);
            }
        }
    }

    std::set<int> unique_scints;
    int n_valid_hits = 0;

    if (!(all_hits.size() < 1))
        detection.HasHit = true;
    if (all_hits.size() > 1)
        detection.HasMultipleHits = true;
    auto n_hits_in_time_window = 0;
    for (const auto &hit : all_hits) {
        if (hit.is_valid_hit) {
            detection.HasValidHit = true;
            n_valid_hits++;
        }
        if (hit.is_valid_hit && hit.is_in_time_window) {
            unique_scints.insert(hit.scintillator_id);
            n_hits_in_time_window++;
        }
        if (hit.quality == HitQuality::OutOfBounds)
            detection.HasOutOfBounds = true;
        if (!hit.is_in_time_window &&
            hit.quality != HitQuality::AccidentalCoincidence)
            detection.HasOutOfTimeWindow = true;
        if (hit.is_in_time_window && hit.is_valid_hit)
            detection.HasInTimeWindow = true;
    }
    if (n_hits_in_time_window == 1)
        detection.IsClean = true;
    if (n_valid_hits > 1)
        detection.HasMultipleValidHits = true;
    detection.n_main_window_events = n_hits_in_time_window;
    detection.HasMultipleScintillatorsHit = (unique_scints.size() > 1);
    detection.T5_hits = all_hits;

    return detection;
}
