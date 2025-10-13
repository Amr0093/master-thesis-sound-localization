🎧 Sound Source Localization – MATLAB/Simulink (TUHH × Valeo)

An end-to-end MATLAB toolbox for acoustic source localization, combining classical beamforming and deep learning (DNN) approaches.
Developed as part of my Master’s thesis in Mechatronics at TUHH, in collaboration with Valeo (ADAS Team).
____________________________________________________________________________________________________________________________________________________________________
🚀 Overview

This project enables simulation, training, and evaluation of 3D sound localization systems under realistic acoustic conditions.

It provides:

Synthetic & real impulse response (IR) generation

Beamforming (CBF / SC-DAMAS / OMP-DAMAS) and DNN-based localization

Configurable microphone array geometries (1D → 3D)

Noise, Doppler, and reverberation modeling

GUI-based end-to-end control

Reproducible results & benchmarks

💡 Start here: AcousticAnalysisGUI.m — this launches the full interactive GUI.
____________________________________________________________________________________________________________________________________________________________________
📂 Repository Structure
Folder	Description
src/	Contains all MATLAB source code, including GUIs, DNN modules, simulation tools, and DSP utilities.
results/	Stores simulation and evaluation outputs — plots, metrics, and benchmark data comparing DNN vs. beamforming.
docs/	Contains the master’s thesis PDF and additional documentation describing the algorithms and experiments.
____________________________________________________________________________________________________________________________________________________________________
🧩 Key Features

Hybrid simulation environment: time-domain, frequency-domain, or mixed processing

DNN models for azimuth, elevation, and distance estimation

Beamforming baselines for comparative benchmarking

Acoustic simulator GUI for controlled dataset generation

Balanced datasets (synthetic + real recordings)

Dockerized environment for reproducibility (outside MATLAB scope)
____________________________________________________________________________________________________________________________________________________________________
🧠 Example Workflow

Generate synthetic IRs / datasets
Run FastHybridDopplerReverbSimulator_GUI.m to define mic geometry, source motion, and acoustic conditions.

Train DNN models
Use one of the DNN_* scripts or the combined 3D model.

Evaluate & compare
Run room_beamforming_and_comparison.m to benchmark beamforming vs. DNN performance.

Full pipeline GUI
Use AcousticAnalysisGUI.m to manage all stages from dataset generation to analysis.
____________________________________________________________________________________________________________________________________________________________________
⚙️ Requirements

MATLAB R2021b+

Toolboxes: Signal Processing, DSP System, Phased Array (recommended), Deep Learning

Optional: Simulink

Optional: Real recordings for IR benchmarking
____________________________________________________________________________________________________________________________________________________________________
🧾 Thesis Reference

This repository accompanies my Master’s thesis:

“Deep Learning-Based Acoustic Source Localization Using Synthetic and Real Impulse Responses.”
The full thesis PDF is located under /docs
____________________________________________________________________________________________________________________________________________________________________

📜 License

This repository is released for academic and non-commercial use.
For other use cases, please contact the author directly.
