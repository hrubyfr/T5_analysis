# T5 detector analysis

This project reconstructs T5 detector hits from WCTE input ROOT files and writes a simplified per-event ROOT tree for downstream analysis.

## Dependencies

The code is built against ROOT 6 (the Makefile uses `root-config`).

Typical setup:

- ROOT 6.x
- C++ compiler with C++17 support
- GNU Make

## Build

From the project root:

```sh
make
```

This produces the executable `analyze_T5`.

## Usage

The current executable accepts the following command-line options:

```sh
./analyze_T5 -r <run_number> -i <input_file.root> [-o <output_dir>] [-d]
```

Options:

- `-r <run_number>`: required run number used for the analysis metadata.
- `-i <input_file.root>`: input ROOT file to analyze. You can pass this option multiple times.
- `-o <output_dir>`: output directory prefix used to write one reconstructed file per input file.
- `-d`: enable debug mode (currently limits the event loop to the first 5000 entries).

Example:

```sh
./analyze_T5 -r 1361 -i /path/to/run_1361.root -o output/
```

The program recognizes the following input trees:

- `ProcessedWaveforms`
- `WCTEReadoutWindows`

If neither tree is found, the program stops with an error.

## What the program does

For each input file, the analysis:

1. Opens the input ROOT file.
2. Detects the appropriate TTree.
3. Reconstructs T5 hits using the T5 SiPM timing information.
4. Writes a new ROOT file containing a simplified event tree.

## Output tree structure

The output ROOT file contains one tree named `T5_Events` with the following branches:

| Branch name | Type | Description |
| :--- | :---: | :--- |
| `event_nr` | `Int_t` | Event number from the input tree |
| `t5_hit_bitmask` | `Int_t` | A bitmask  |
| `t5_n_main_bunch_particles` | `Int_t` | Number of reconstructed main-bunch hits written for this event |
| `t5_main_hit_time` | `Double_t` | Time of the primary reconstructed hit |
| `t5_main_hit_charge` | `Double_t` | Per-bar charge of the primary reconstructed hit in raw ADC, calculated as the geometric mean of the left and right SiPM charges; not calibrated to p.e. |
| `t5_main_hit_sipm_time_left` | `Double_t` | Time measured by the left SiPM of the primary hit (SiPM indices 8-15) |
| `t5_main_hit_sipm_time_right` | `Double_t` | Time measured by the right SiPM of the primary hit (SiPM indices 0-7) |
| `t5_main_hit_sipm_charge_left` | `Double_t` | Raw ADC charge measured by the left SiPM of the primary hit (SiPM indices 8-15), not p.e. |
| `t5_main_hit_sipm_charge_right` | `Double_t` | Raw ADC charge measured by the right SiPM of the primary hit (SiPM indices 0-7), not p.e. |
| `t5_main_hit_pos_x` | `Double_t` | Reconstructed x position of the primary hit [mm] |
| `t5_main_hit_pos_y` | `Double_t` | Reconstructed y position of the primary hit [mm] |
| `t5_main_hit_pos_x_error` | `Double_t` | Uncertainty of the reconstructed x position of the primary hit [mm] |
| `t5_main_hit_pos_y_error` | `Double_t` | Uncertainty of the reconstructed y position of the primary hit [mm] |
| `t5_all_hits_pos_x` | `vector<double>` | x positions of all reconstructed hits in an event |
| `t5_all_hits_pos_y` | `vector<double>` | y positions of all reconstructed hits in an event |
| `t5_all_hits_pos_x_error` | `vector<double>` | x positions of all reconstructed hits in an event |
| `t5_all_hits_pos_y_error` | `vector<double>` | y positions of all reconstructed hits in an event |
| `t5_all_hits_time` | `vector<double>` | Times of all reconstructed hits in an event |
| `t5_all_hits_charge` | `vector<double>` | Per-bar charges of all reconstructed hits in raw ADC, calculated as the geometric mean of the two SiPM charges; not calibrated to p.e. |
| `t5_all_hits_sipm_time_left` | `vector<double>` | Left SiPM times for all reconstructed hits (SiPM indices 8-15) |
| `t5_all_hits_sipm_time_right` | `vector<double>` | Right SiPM times for all reconstructed hits (SiPM indices 0-7) |
| `t5_all_hits_sipm_charge_left` | `vector<double>` | Raw ADC left SiPM charges for all reconstructed hits (SiPM indices 8-15), not p.e. |
| `t5_all_hits_sipm_charge_right` | `vector<double>` | Raw ADC right SiPM charges for all reconstructed hits (SiPM indices 0-7), not p.e. |

## T5 hit bitmask

The output file has a T5 hit bitmask variable, that describes the error states that occured in the event. It can have the following values:

- `0` -- The event is clean, has one main hit reconstructed in the T5 detector bounds (or is within uncertainty), there are no additional hits in the main bunch of particles
- `1` -- The event is empty -- there was no hit in the main bunch of particles
- `2` -- The event has multiple hits in the main time window -- a 'main' hit cannot be reliably chosen
- `4` -- The event had a single hit in the main time window, but it was not successfully reconstructed within T5 bounds -- could indicate an SiPM time measurement error, or just an invalid hit.

## Notes

- The program currently writes one output file per input file, using the input file basename with `_T5.root` appended.
- The output directory path passed to `-o` is used as the base location for those generated files.
- The reconstruction logic is implemented mainly in `return_TOF_position.cpp` and uses the detector geometry constants defined in `return_TOF_position.h`.
