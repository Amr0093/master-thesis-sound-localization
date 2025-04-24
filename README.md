# master-thesis-sound-localization
This folder will contain the full documentation and implementation of my Master’s thesis, focused on acoustic source localization using vehicle-mounted microphone arrays in complex urban environments. The work addresses challenges such as Doppler shift, reverberation, and low signal-to-noise ratio through advanced beamforming algorithms. Full content will be uploaded upon official submission.

🧠 SC-DAMAS Simulation: Acoustic Imaging in Reverberant Fields
This MATLAB-based simulation implements the beamforming and deconvolution techniques described in the thesis paper:


It supports:

Room geometry setup

Room impulse response (RIR) modeling using the mirror source method

Signal simulation with added noise

Acoustic source imaging via CBF, OMP-DAMAS, and SC-DAMAS

🔁 Workflow Overview

A[1. room_geometry_setup.m] --> B[2. rir_generator.m (custom)]
B --> C[3. simulate_microphone_signals.m]
C --> D[4. sc_damas.m (beamforming + deconvolution)]

📂 File Descriptions
room_geometry_setup.m ✅
Purpose:
Defines the dimensions of the room, source location, and sensor array. Calculates frequency-dependent reflection coefficients for each wall.

Output:

room_geometry_setup.mat (contains: roomDim, sourcePos, micPositions, reflection_coef, etc.)

rir_generator.m (custom mirror source implementation) ✅
Purpose:
Implements the mirror image source method to compute the room impulse response (RIR) and its frequency-domain version.

Output:

room_impulse_response.mat (contains: impulse_responses, H, freq_vector, t, etc.)

simulate_microphone_signals.m ✅
Purpose:
Convolves a test signal (chirp) with the RIRs to simulate what each microphone "hears", and adds noise to simulate low SNR conditions.

Output:

simulated_microphone_signals.mat (contains: clean_signals, noisy_signals, SNR_dB, fs, t)

sc_damas.m ✅
Purpose:
Performs acoustic source localization using:

CBF (Conventional Beamforming)

OMP-DAMAS

SC-DAMAS (Sparse CoSaMP-based variant)

Visualizes results as beamformed power maps and compares accuracy (vs. ground truth).

Output:

sc_damas_results.mat (beamforming results and localization metrics)

📊 Visualization Matches Paper
The generated figures match:

Fig. 3: Sound pressure distributions

Figs. 4–6: Acoustic maps of beamforming methods

Table 3: Localization performance

⚙️ Requirements
MATLAB R2020+ (earlier versions should work too)

Signal Processing Toolbox (for chirp, fft, etc.)

Plots and saving handled natively (no additional dependencies)

✅ Status

Step	Task	Status
1	Room geometry setup	✅ Complete
2	Generate room impulse responses	✅ Complete
3	Simulate received signals (with noise)	✅ Complete
4	Implement SC-DAMAS	✅ Complete
5	Compare with CBF and OMP-DAMAS	✅ Complete
