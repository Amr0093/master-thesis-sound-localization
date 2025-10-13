🧠 Source Code Overview (/src)

This folder contains all MATLAB source files used in the Sound Source Localization – TUHH × Valeo project.
It includes the simulation environment, data generation tools, classical DSP baselines, and deep learning (DNN) models.
____________________________________________________________________________________________________________________________________________________________________
🎛️ GUI Entry Points
File	Description
AcousticAnalysisGUI.m	Main control panel. Launch this to access the full pipeline — dataset generation, DNN training/inference, beamforming, benchmarking, and visualization. It serves as the project’s front-end and integrates all other modules.
FastHybridDopplerReverbSimulator_GUI.m	Acoustic simulation GUI. Generates synthetic impulse responses (IRs) and datasets under configurable conditions — mic array geometry (1D–3D), source motion, SNR, noise type, and reverberation. Can run standalone or be called from the main GUI.
____________________________________________________________________________________________________________________________________________________________________
🧩 Dataset & IR Utilities
File	Description
generateAcousticDataset.m	Builds balanced datasets by combining synthetic and real IRs across multiple environments. Supports automatic labeling for DNN training.
IRs_comparison_analysis_tool3.m	Aligns synthetic vs. real recordings in time/frequency domains. Computes similarity metrics (e.g., correlation, spectral distance) to assess IR realism.
reflection_understanding.m	Utility for modeling and visualizing reflections and reverberation behavior to better understand acoustic propagation effects.
🎧 Beamforming & Classical DSP
File	Description
room_beamforming_and_comparison.m	Implements and compares beamforming algorithms (CBF, SC-DAMAS, OMP-DAMAS) for azimuth, elevation, and distance estimation. Produces performance metrics and plots for benchmarking against DNN models.
____________________________________________________________________________________________________________________________________________________________________
🤖 Deep Learning Models (Training & Inference)
File	Description
DNN_azimuth_est_clean_noisy1.m	Trains and evaluates deep neural networks for azimuth estimation using clean and noisy datasets.
DNN_elevation_est_clean_noisy1.m	DNN for elevation angle estimation under varying noise levels and SNR conditions.
DNN_distance_est_clean_noisy1.m	DNN for source distance prediction based on IRs or convolved signals.
DNN_combined_3D_localization.m	Joint 3D model predicting azimuth, elevation, and distance simultaneously. Integrates all previous models for end-to-end localization.
____________________________________________________________________________________________________________________________________________________________________
🔬 Typical Workflow (Inside /src)

Run FastHybridDopplerReverbSimulator_GUI.m
→ Generate synthetic IRs and configure microphone geometry, noise, and environment parameters.

Optionally run generateAcousticDataset.m
→ Build structured datasets for DNN training or hybrid testing.

Train and evaluate models:

DNNs → Run one of the DNN_* scripts.

Classical → Run room_beamforming_and_comparison.m.

Compare results in the main GUI (AcousticAnalysisGUI.m)
→ Plot, analyze, and export metrics.
____________________________________________________________________________________________________________________________________________________________________
⚙️ Notes

All scripts are modular — you can run them independently or through the main GUI.

The simulator supports time-domain, frequency-domain, and hybrid signal processing pipelines.

DNNs can be trained on synthetic data, real data, or a mix (transfer learning ready).
