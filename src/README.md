Sound Source Localization – MATLAB/Simulink (TUHH × Valeo)

End-to-end toolbox for simulating, training, and evaluating acoustic source localization using:

Synthetic + real impulse responses (IRs)

Classical beamforming (CBF / SC-DAMAS / OMP-DAMAS)

Deep neural networks (DNNs) for azimuth, elevation, and distance

GUI front-ends for both simulation and analysis

Tip: Start with AcousticAnalysisGUI.m (main launcher).
The simulator GUI (FastHybridDopplerReverbSimulator_GUI.m) can be run standalone or from the main GUI.
____________________________________________________________________________________________________________________________________________________________________________________________________________________
File map
GUIs (entry points)

AcousticAnalysisGUI.m
Main control panel. Orchestrates dataset generation, DNN training/inference, beamforming pipelines, and IR benchmarking. Provides access to all modules below.

FastHybridDopplerReverbSimulator_GUI.m
Acoustic scene simulator. Generates signals/IRs under configurable conditions:
mic-array geometry (1D → 3D), source motion (Doppler), reverberation, noise type/SNR, time/frequency/hybrid processing, and IR convolution/deconvolution.
____________________________________________________________________________________________________________________________________________________________________________________________________________________
Dataset & IR tools

generateAcousticDataset.m
Builds balanced datasets (synthetic + real) by stacking/combining IRs across scenarios.

IRs_comparison_analysis_tool3.m
Aligns and compares synthetic vs. real recordings (time/frequency-domain metrics) to validate simulation fidelity.

reflection_understanding.m
Helpers for modeling/inspecting reflections & reverberation characteristics.
____________________________________________________________________________________________________________________________________________________________________________________________________________________
Beamforming & classical DSP

room_beamforming_and_comparison.m
Runs CBF, SC-DAMAS, and OMP-DAMAS for AoA (azimuth/elevation) and distance estimation; benchmarks versus DNN outputs.
____________________________________________________________________________________________________________________________________________________________________________________________________________________
Deep learning (training/inference)

DNN_azimuth_est_clean_noisy1.m
Train/infer DNN for azimuth localization on clean/noisy data.

DNN_elevation_est_clean_noisy1.m
Train/infer DNN for elevation localization on clean/noisy data.

DNN_distance_est_clean_noisy1.m
Train/infer DNN for distance estimation.

DNN_combined_3D_localization.m
Joint model for azimuth + elevation + distance (3D localization).
____________________________________________________________________________________________________________________________________________________________________________________________________________________
Typical workflows
1) Generate synthetic data / IRs

Run FastHybridDopplerReverbSimulator_GUI.m.

Choose array geometry (1D/2D/3D), source motion, SNR/noise, reverberation model, and processing mode (time/freq/hybrid).

Export audio/IRs or use generateAcousticDataset.m to build a balanced dataset.

2) Validate sim vs. real recordings

Use IRs_comparison_analysis_tool3.m to align and benchmark synthetic vs. real recordings (FFT/DSP metrics, IR similarity).

3) Train & evaluate models

Beamforming: room_beamforming_and_comparison.m

DNNs: one of the DNN_* scripts (or the combined model).
Compare accuracy/latency/robustness; hybridize if helpful (e.g., OMP-DAMAS for AoA + DNN for distance).

4) End-to-end control

Launch AcousticAnalysisGUI.m to run the full pipeline via GUI (dataset → training/inference → comparisons → plots/exports).
____________________________________________________________________________________________________________________________________________________________________________________________________________________
Key features

Mic arrays with customizable geometry (linear, planar, volumetric)

Doppler and reverberation modeling

SNR & noise type controls; signal conditioning and filtering

IR convolution/deconvolution

Beamforming baselines vs. DNN models

Balanced synthetic + real datasets

Plots, metrics, and reproducible exports
____________________________________________________________________________________________________________________________________________________________________________________________________________________
Requirements

MATLAB (R2021b+ recommended)
Toolboxes: Signal Processing, DSP System, Phased Array (recommended), Deep Learning.

Optional: Simulink (for some blocks/workflows).

Real recordings (optional) for benchmarking.
____________________________________________________________________________________________________________________________________________________________________________________________________________________
Quick start

Open MATLAB in this folder.

Run:

AcousticAnalysisGUI
% or, to build synthetic data first
FastHybridDopplerReverbSimulator_GUI


Follow on-screen steps to generate data, train, and evaluate.
____________________________________________________________________________________________________________________________________________________________________________________________________________________
Notes

The simulator supports both time-domain, frequency-domain, and hybrid pipelines.

DNN scripts accept either raw audio or GAF/feature representations, depending on configuration inside each file.

For reproducibility across environments, packaging via Docker was used in the broader project (outside pure MATLAB scope).

Citation / Thesis

This code accompanies my TUHH × Valeo master’s thesis on acoustic source localization with synthetic/real IRs, beamforming baselines, and DNNs.
The full thesis document is linked in my LinkedIn profile under Education → TUHH.
____________________________________________________________________________________________________________________________________________________________________________________________________________________
License
For academic evaluation and non-commercial use. Contact the author for other licensing terms.
