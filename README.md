# T5 detector analysis

This is the code to analyze the T5 detector using the matched VME + BRB data files of the WCTE experiment.

## How to run

To run this code, compile it using the 
```
make
```
command, and run the analysis using 
```
./analyze_T5 -r <run_number>
```
command

The program can take additional arguments, `-o` will modify the name of the output root file and output path, `-i` will modify, where the program searches for the files

## Output file data structure

The root output file has the following structure:

| Variable Name | Type | Description |
| :--- | :---: | :--- |
| **event_nr** | `Int` | Unique identification number for the physics event. |
| **T5_particle_nr** | `Int` | Total number of particles detected by T5 in this event -- both in and out of the main time window |
| **T5_HasValidHit** | `Bool` | Flag indicating if a valid signal was recorded. |
| **T5_HasMultipleScintillatorsHit** | `Bool` | True if the event had particles in the defined main time window going through more than one scintillator -- multiple particles in one bunch, or secondary particles created by interactions. |
| **T5_HasOutOfTimeWindow** | `Bool` | Flag saying, if the event has hits falling outside the defined trigger window. |
| **T5_HasInTimeWindow** | `Bool` | Flag saying, if the event has hits falling within the defined trigger window. |
| **T5_hit_is_in_bounds** | `vector<int>` | Validation that the hit is within physical detector limits (1 = true, 0 = false). |
| **T5_hit_pos_x** | `vector<float>` | X-coordinate of the hit in the T5 detector [mm]. |
| **T5_hit_pos_y** | `vector<float>` | Y-coordinate of the primary hit in the T5 detector [mm]. |
| **T5_hit_time** | `vector<float>` | Time measured by T5 for the primary hit. |
| **T5_secondary_hit_is_in_bounds** | `vector<int>` | Says if the hit out of the main time window is in the scintillator T5_secondary_hit_is_in_bounds|
| **T5_secondary_hit_pos_x** | `vector<float>` | X-coordinate of the secondary hit -- hit outside of the main time window |
| **T5_secondary_hit_pos_y** | `vector<float>` | Y-coordinate of the secondary hit -- hit outside of the main time window |
| **T5_secondary_hit_time** | `vector<float>` | Time measured by T5 for the secondary hit -- hit outside of the main time window |
