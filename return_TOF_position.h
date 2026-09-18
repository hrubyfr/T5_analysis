#ifndef RETURN_TOF_POSITION_H
#define RETURN_TOF_POSITION_H

#include <array>
#include <iostream>
#include <math.h>
#include <utility>
#include <vector>

namespace T5_CONFIG { // a namespace storing the parameters of the T5 detector
inline constexpr int T5_MPMT_ID = 132; // Beam monitors readout board number --
                                       // board 132 had T5 and T4 detectors
inline constexpr int T5_TRIGGER_ID =
    19; // Trigger was connected to board 132 to channel 19
inline constexpr int N_T5_SCINTS = 8; // 8 scintillator blocks in T5
inline constexpr int N_T5_SIPMS = 16; // Each block had 2 SiPMs connected
inline constexpr double SCINT_BLOCK_HEIGHT =
    16.25; // Height of the scintillator blocks
inline constexpr double EXPECTED_DETECTION_TIME_MIN = -190;
inline constexpr double EXPECTED_DETECTION_TIME_MAX = -120;

// Shift in measured time in different scintillators due to detector bias
inline constexpr std::array<double, 8> SCINT_BIAS = {
    -0.0253795, 0.307779,  -0.0568409, 0.437488,
    -0.137346,  -0.350588, -0.153047,  -0.374798};
// y-coordinate of the middle of a scintillator block
inline constexpr std::array<double, 8> SCINT_Y_POSITIONS = {
    60.025, 42.875, 25.725, 8.575, -8.575, -25.725, -42.875, -60.025};
// length of the 8 scintillator blocks
inline constexpr std::array<double, 8> SCINT_DIMENSIONS = {
    40.92, 94.0, 112.0, 123.0, 123.0, 112.0, 94.0, 40.92};
// a map of corresponding channels, connected to a single scintillator
inline constexpr std::array<std::pair<int, int>, 16> T5_SIPM_INDICES = {
    {{0, 0},
     {1, 1},
     {2, 2},
     {3, 3},
     {4, 4},
     {5, 5},
     {6, 6},
     {7, 7},
     {8, 8},
     {10, 9},
     {11, 10},
     {12, 11},
     {13, 12},
     {14, 13},
     {15, 14},
     {16, 15}}};
// Channel numbers of T5 SiPMs
inline constexpr std::array<int, 16> T5_IDS = {0, 1,  2,  3,  4,  5,  6,  7,
                                               8, 10, 11, 12, 13, 14, 15, 16};
constexpr int GetSiPMIndex(int pmt_id) {
    for (const auto &pair : T5_SIPM_INDICES) {
        if (pair.first == pmt_id)
            return pair.second;
    }
    return -1; // Indicates invalid/not found
}
} // namespace T5_CONFIG

enum class HitQuality {
    Perfect,     // Hit is inside the scintillator bars
    OutOfBounds, // Hit is outside the bounds, but within a 3-sigma error margin
    AccidentalCoincidence // Way out of bounds -- completely different
                          // particles, or comparing with dark noise
};
struct T5_hit {
    // Tells if the hit was out of bounds or was reconstructed inside
    bool is_valid_hit = false;
    bool is_in_time_window = false;
    // Says where the event was reconstructed:
    // 	Perfect = inside of a scintillator;
    //	OutOfBounds = outside of scintillator dimensions;
    //	Accidental coincidence = Way out of the scintillator dimensions --
    // comparison of two completely different times/different particles
    HitQuality quality = HitQuality::AccidentalCoincidence;
    double position_x = -999;
    double position_y = -999;  // X and Y coordinates of the T5 hit
    double uncertainty = -999; // Position uncertainty, assumed that v_eff and
                               // sigma_sipm are uncorrelated
    double hit_time = -999; // Time of the detection, calculated as average time
                            // of the two SiPMs
    double sipm_time_a = -999;
    double sipm_time_b = -999;
    double raw_time = -999;
    double hit_charge = -999;
    int scintillator_id = -1;
};

struct event_T5_detection {
    // Event number
    int event_nr = -1;
    // The number of hits in the main time window (i.e. time minus trigger
    // between expected time window max and expected time window min)
    int n_main_window_events = -1;
    // Says, if the event is clean (i.e. has only a single hit in the main time
    // window)
    bool IsClean = false;
    // Says, if the event has at least one hit -- can be an accidental
    // coincidence
    bool HasHit = false;
    // Does event have more than one hit (vector of T5 hits is larger than 1)
    bool HasMultipleHits = false;
    // Does the event have one valid hit -- either in the scintillator, or out
    // of OutOfBounds
    bool HasValidHit = false;
    // Are there more than one hits, that are valid?
    bool HasMultipleValidHits = false;
    // Were more than one scintillators hit in a single event? -- Hit was in the
    // expected timeframe, and was valid
    bool HasMultipleScintillatorsHit = false;
    // Does the event have a reconstruction out of bounds?
    bool HasOutOfBounds = false;
    // Does the event have at least one event out of the expected time frame?
    bool HasOutOfTimeWindow = false;
    // Does the event have at least one hit in the expected time window? (main
    // time peak -- between -140 and -170 ns)
    bool HasInTimeWindow = false;
    // A vector of the final T5 hits
    std::vector<T5_hit> T5_hits;
};

// The class TOF reconstructor has several useful functions. The main one is the
// Return_position, which takes the data, analyzes it and returns a vcector of
// the T5 hit structure that contains the validity of the hit, the position,
// uncertainty of the hit and the hit time of all hits. It also has the method
// HasMultiHits, which checks if the event has multiple valid hits (i.e. several
// hits in the expected timeframe of time - trigger between -140 and -170 ns).
// This method looks at the vector of T5 hits returned by Return_position. Then
// there is IsEventValid method, which checks, if there are any valid hits
// inside of the event.

class TOF_reconstructor {
  public:
    explicit TOF_reconstructor(double v_eff = 181.974);
    void SetVeff(double v);
    void SetVeffUncertainty(double uncertainty);
    void SetVerbosity(int i);

    double GetVeff() const;
    double GetVeff_uncertainty() const;
    bool GetVerbosity() const;
    double GetScintDimensionX(int i) const;
    double GetScintPositionY(int i) const;
    double Get_scint_xmin(int i) const;
    double Get_scint_xmax(int i) const;
    double Get_ymax() const;
    double Get_ymin() const;

    event_T5_detection
    Return_position(const int event_nr, const std::vector<int> &hit_mpmt_ids,
                    const std::vector<int> &hit_pmt_ids,
                    const std::vector<double> &hit_pmt_times,
                    const std::vector<double> &hit_pmt_charges);

  private:
    double _v_eff;
    double v_eff_uncertainty;
    // Scintillator time uncertainty, acquired through minimizing chi2 function
    std::array<double, 8> sigma_sipm_i;
    // uncertainty of the scintillator resolution
    std::array<double, 8> sigma_sipm_i_uncertainties;

    bool _verbose;
};

#endif
