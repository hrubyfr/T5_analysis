#include "./return_TOF_position.h"
#include <vector>

using namespace std;
using namespace T5_CONFIG;

// TOF_reconstructor class constructor -- private values contain only fit values -- can change depending on 
TOF_reconstructor::TOF_reconstructor(double v_eff) : _v_eff(v_eff),
				v_eff_uncertainty(33.7644),
				sigma_sipm_i({0.333076, 0.360282, 0.356124, 0.304019, 0.226525, 0.243068, 0.306572, 0.273243}),
				sigma_sipm_i_uncertainties({0.0101467, 0.0490062, 0.0703601, 0.0993561, 0.133235, 0.103015, 0.0575791, 0.0123247}),

				_verbose(0)
{}


void TOF_reconstructor::SetVeff(double v){ _v_eff = v; }
void TOF_reconstructor::SetVerbosity(int i){ _verbose = i; }

double TOF_reconstructor::GetVeff() const { return _v_eff; }
bool TOF_reconstructor::GetVerbosity(){ return _verbose;}

double TOF_reconstructor::GetScintDimensionX(int i){
	return SCINT_DIMENSIONS[i];
}
double TOF_reconstructor::GetScintPositionY(int i){
	return SCINT_Y_POSITIONS[i];
}
double TOF_reconstructor::Get_scint_xmax(int i){
	return SCINT_DIMENSIONS[i]/2;
}
double TOF_reconstructor::Get_scint_xmin(int i){
	return -SCINT_DIMENSIONS[i]/2;
}
double TOF_reconstructor::Get_ymax(){
	double scintillator_block_halfheight = SCINT_BLOCK_HEIGHT / 2.0;
	return SCINT_Y_POSITIONS[0] + scintillator_block_halfheight;
}
double TOF_reconstructor::Get_ymin(){
	double scintillator_block_halfheight = SCINT_BLOCK_HEIGHT / 2.0;
	return SCINT_Y_POSITIONS[7] - scintillator_block_halfheight;
}
event_T5_detection TOF_reconstructor::Return_position(const int event_nr,
						      const vector<int>& hit_mpmt_ids,
						      const vector<int>& hit_pmt_ids, 
						      const vector<double>& hit_pmt_times) {

	event_T5_detection detection;
	detection.event_nr = event_nr;
	vector<T5_hit> all_hits;

	// Create 16 vectors where to store all the times detected by T5 SiPMs, some will be empty, some will have multiple hits, ideally paired -- later will check for hits in paired detectors
	vector<vector<double>> T5_times(N_T5_SIPMS);	
	bool valid_hit = false;
	double trigger_time = 0;
	for (size_t i = 0; i < hit_pmt_ids.size(); i++){
		if (hit_mpmt_ids.at(i) == T5_MPMT_ID && hit_pmt_ids.at(i) == T5_TRIGGER_ID){
			// Take the first trigger event (should be only one, hopefully), and set it as trigger time, then break the for loop
			trigger_time = hit_pmt_times.at(i);	
			break;
		}
	}
	for (size_t i = 0; i < hit_mpmt_ids.size(); i++){
		auto mPMT_id = hit_mpmt_ids.at(i);	
		if (mPMT_id != T5_MPMT_ID) continue;
		auto PMT_id = hit_pmt_ids.at(i);
		auto SiPM_index = GetSiPMIndex(PMT_id);
		if (SiPM_index == -1) continue;
		// correct the measured T5 sipm time by the trigger time
		T5_times.at(SiPM_index).push_back(hit_pmt_times.at(i) - trigger_time);
	}
	// Loop over all saved times in the corresponding vectors, compare all the times
	for (int i = 0; i < N_T5_SCINTS; i++){	
		if (T5_times[i].empty() || T5_times[i+8].empty()) continue;
		for (const auto& sipm_time_a : T5_times[i]){
			if (sipm_time_a > EXPECTED_DETECTION_TIME_MAX || sipm_time_a < EXPECTED_DETECTION_TIME_MIN) {
				if (_verbose) cout << "WARNING: SiPM time is out of expected hit times" << endl; 
				// warning that the detected time is outside the expected event peak, as detected in low energy events -- not tested in high energy events
			}
			for (const auto& sipm_time_b : T5_times[i+8]){
				T5_hit hit;
				hit.sipm_time_a = sipm_time_a;
				hit.sipm_time_b = sipm_time_b;
				if (sipm_time_b > EXPECTED_DETECTION_TIME_MAX || sipm_time_b < EXPECTED_DETECTION_TIME_MIN) {
					if (_verbose) cout << "WARNING: SiPM time is out of expected hit times" << endl;
				}
				auto time_diff = sipm_time_a - sipm_time_b;
				double avg_time = (sipm_time_a + sipm_time_b) / 2;
				hit.hit_time = avg_time;
				valid_hit = true;
				pair<double,double> position;
				position.first = (time_diff - SCINT_BIAS.at(i)) * _v_eff / 2.0;
				position.second = SCINT_Y_POSITIONS[i];

				if (abs(position.first) > SCINT_DIMENSIONS[i]/2){
					if (_verbose)cout << "Reconstructed position is out of bounds" << endl;
					// If the reconstruction is outside of the corresponding scintillator bounds, set hit to false, return error value
					position.first = -999;
					valid_hit = false;
				}
				// _verbose checks
				if (_verbose) cout << "X coordinate is:\t" << position.first;

				//calculate the position uncertainty
				double uncertainty = sqrt(pow(v_eff_uncertainty, 2) * pow(time_diff/2, 2) + pow(sigma_sipm_i[i], 2) * pow(_v_eff/2, 2));
				if (_verbose) cout << "\tX uncertainty is:\t" << uncertainty << endl;
				hit.uncertainty = uncertainty;
				hit.position = position;
				hit.valid_hit = valid_hit;
				hit.scintillator_id = i;
				// save the hit to the vector of all hits in the event
				all_hits.push_back(hit);
			}
		}
	}

	if (all_hits.size() != 1 && _verbose)cout << "WARNING: More than 1 scintillator was hit" << endl;

	detection.T5_hits = all_hits;

	return detection;
}
bool TOF_reconstructor::HasMultiValidHits(const std::vector<T5_hit>& hits) {
	// If the vector is too small, we can return false immediately
	if (hits.size() < 2) return false;

	int valid_count = 0;
	for (const auto& hit : hits) {
		if (hit.valid_hit) {
			valid_count++;
		}
		// "Smart" Optimization: Stop as soon as we confirm > 1
		if (valid_count > 1) {
			return true;
		}
	}
	return false;
}

int TOF_reconstructor::HasWeirdHits(const std::vector<T5_hit>& hits){
	// return 0 if event has hits reconstructed outside of the scintillator scope, 1 if the event has multiple hits, -1 if event is empty (no pair hit)
	for (const auto& hit : hits){
		if (!hit.valid_hit && hit.position) return 0;
	}
	if (hits.size() > 1) return 1;
	return -1;
}
bool TOF_reconstructor::HasValidHits(const std::vector<T5_hit> &hits){
	// checks if the event has at least one valid hit
	for (const auto& hit : hits){
		if (hit.valid_hit) return true;
	}
	return false;
}

