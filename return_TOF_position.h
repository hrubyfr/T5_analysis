#ifndef RETURN_TOF_POSITION_H
#define RETURN_TOF_POSITION_H

#include <optional>
#include <utility>
#include <array>
#include <vector>
#include <set>
#include <iostream>
#include <math.h>

namespace T5_CONFIG { // a namespace storing the parameters of the T5 detector
		inline constexpr int T5_MPMT_ID = 132; //Beam monitors readout board number -- board 132 had T5 and T4 detectors
		inline constexpr int T5_TRIGGER_ID = 19; //Trigger was connected to board 132 to channel 19
		inline constexpr int N_T5_SCINTS = 8; //8 scintillator blocks in T5
		inline constexpr int N_T5_SIPMS = 16; //Each block had 2 SiPMs connected
		inline constexpr double SCINT_BLOCK_HEIGHT = 16.25; //Height of the scintillator blocks
		inline constexpr double EXPECTED_DETECTION_TIME_MIN = -170; 
		inline constexpr double EXPECTED_DETECTION_TIME_MAX = -140;

		//Shift in measured time in different scintillators due to detector bias 
		inline constexpr std::array<double, 8> SCINT_BIAS = {
			-0.0253795, 0.307779, -0.0568409, 0.437488, -0.137346, -0.350588, -0.153047, -0.374798
		};
		//y-coordinate of the middle of a scintillator block 
		inline constexpr std::array<double, 8> SCINT_Y_POSITIONS = {
			60.025, 42.875, 25.725, 8.575, -8.575, -25.725, -42.875, -60.025
		};
		//length of the 8 scintillator blocks	
		inline constexpr std::array<double, 8> SCINT_DIMENSIONS = {
			40.92, 94.0, 112.0, 123.0, 123.0, 112.0, 94.0, 40.92	
		};
		//a map of corresponding channels, connected to a single scintillator
		inline constexpr std::array<std::pair<int, int>, 16> T5_SIPM_INDICES = {{
			{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}, {7, 7},
			{8, 8}, {10, 9}, {11, 10}, {12, 11}, {13, 12}, {14, 13}, {15, 14}, {16, 15}
		}};
		//Channel numbers of T5 SiPMs
		inline constexpr std::array<int, 16> T5_IDS = {
			0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16	
		};	
		constexpr int GetSiPMIndex(int pmt_id) {
			for (const auto& pair : T5_SIPM_INDICES) {
				if (pair.first == pmt_id) return pair.second;
			}
			return -1; // Indicates invalid/not found
		}
}

enum class HitQuality{
	Perfect, // Hit is inside the scintillator bars
	OutOfBounds, // Hit is outside the bounds, but within a 3-sigma error margin
	AccidentalCoincidence // Way out of bounds -- completely different particles, or comparing with dark noise
};
struct T5_hit{
	bool is_valid_hit = false; // Tells if the hit was out of bounds or was reconstructed inside
	bool is_in_time_window = false;
	HitQuality quality = HitQuality::AccidentalCoincidence;
	std::optional<std::pair<double, double>> position; //X and Y coordinates of the T5 hit
	std::optional<double> uncertainty; //Position uncertainty, assumed that v_eff and sigma_sipm are uncorrelated
	std::optional<double> hit_time; // Time of the detection, calculated as average time of the two SiPMs
	std::optional<double> sipm_time_a;
	std::optional<double> sipm_time_b;
	std::optional<int> scintillator_id;
};

struct event_T5_detection{
	int event_nr = -1;
	bool HasHit = false;
	bool HasMultipleHits = false;
	bool HasValidHit = false;
	bool HasMultipleValidHits = false;
	bool HasMultipleScintillatorsHit = false;
	bool HasOutOfBounds = false;
	bool HasOutOfTimeWindow = false;
	bool HasInTimeWindow = false;
	std::vector<T5_hit> T5_hits;
};

// The class TOF reconstructor has several useful functions. The main one is the Return_position, which takes the data, analyzes it and returns a vcector of the T5 hit structure that contains the validity of the hit, the position, uncertainty of the hit and the hit time of all hits. It also has the method HasMultiHits, which checks if the event has multiple valid hits (i.e. several hits in the expected timeframe of time - trigger between -140 and -170 ns). This method looks at the vector of T5 hits returned by Return_position. Then there is IsEventValid method, which checks, if there are any valid hits inside of the event. 

class TOF_reconstructor{
	public:
		explicit TOF_reconstructor(double v_eff = 181.974);
		void SetVeff(double v);
		void SetVerbosity(int i);

		double GetVeff() const;
		bool GetVerbosity();
		double GetScintDimensionX(int i);
		double GetScintPositionY(int i);
		double Get_scint_xmin(int i);
		double Get_scint_xmax(int i);
		double Get_ymax();
		double Get_ymin();

		event_T5_detection Return_position(const int event_nr,
				const std::vector<int>& hit_mpmt_ids,
				const std::vector<int>& hit_pmt_ids, 
				const std::vector<double>& hit_pmt_times);
		static bool HasMultiValidHits(const std::vector<T5_hit>& hits);
		static int HasWeirdHits(const std::vector<T5_hit>& hits);
		static bool HasValidHits(const std::vector<T5_hit>& hits);
	private:
		double _v_eff;
		const double v_eff_uncertainty;
		//Scintillator time uncertainty, acquired through minimizing chi2 function
		std::array<double, 8> sigma_sipm_i;
		// uncertainty of the scintillator resolution
		std::array<double, 8> sigma_sipm_i_uncertainties;

		bool _verbose;
};



#endif 
