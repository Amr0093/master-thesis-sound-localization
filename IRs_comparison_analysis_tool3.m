%% Combined Template-Aligned Impulse Response Comparison Tool
% Performs template matching alignment and comprehensive acoustic analysis
% Combines templateMatchTrim() and IRs_comparison_analysis_tool3() functionality

function IRs_comparison_analysis_tool3()
    close all;
    clear;
    clc;
    
    fprintf('========================================================\n');
    fprintf('    COMBINED TEMPLATE-ALIGNED IR COMPARISON TOOL       \n');
    fprintf('========================================================\n\n');
    
    %% =================================================================
    %% PHASE 1: TEMPLATE MATCHING AND ALIGNMENT
    %% =================================================================
    fprintf('=== PHASE 1: TEMPLATE MATCHING AUDIO TRIMMER ===\n\n');
    
    % === GENERATE REFERENCE SIGNAL ===
    fprintf('Generating reference signal...\n');
    fs_ref = 48000;
    duration_ref = 3.0;
    freqStart = 20;
    freqEnd = 20000;
    reference_signal = generateExpSweep(fs_ref, duration_ref, freqStart, freqEnd);
    
    % === EXTRACT DISTINCTIVE TEMPLATE ===
    fprintf('Extracting distinctive template...\n');
    
    % Find the most energetic region (around 500-2000 Hz sweep)
    template_duration = 0.9; % 900ms template
    template_length = round(template_duration * fs_ref);
    
    % Extract template from middle region of sweep (most distinctive)
    start_template = round(0.3 * length(reference_signal)); % Skip fade-in
    end_template = start_template + template_length - 1;
    template = reference_signal(start_template:end_template);
    
    fprintf('  Template: %.2f s (samples %d to %d)\n', template_duration, start_template, end_template);
    fprintf('  Template energy: %.4f\n', rms(template));
    
    % === LOAD SIGNALS ===
    fprintf('\nLoading signals...\n');
    
    % Load recorded signal
    [recorded_signal, fs_rec, recorded_found] = loadRecordedSignal();
    
    % Load simulated signal  
    [simulated_signal, fs_sim, simulated_found] = loadSimulatedSignal();
    
    if ~recorded_found || ~simulated_found
        error('Could not load both recorded and simulated signals');
    end
    
    % === UNIFY SAMPLING RATES ===
    target_fs = max([fs_rec, fs_sim, fs_ref]);
    fprintf('Unifying to %.0f Hz sampling rate...\n', target_fs);
    
    if fs_rec ~= target_fs
        recorded_signal = resample(recorded_signal, target_fs, fs_rec);
    end
    if fs_sim ~= target_fs
        simulated_signal = resample(simulated_signal, target_fs, fs_sim);
    end
    if fs_ref ~= target_fs
        template = resample(template, target_fs, fs_ref);
        template_length = length(template);
    end
    
    % === TEMPLATE MATCHING ===
    fprintf('\nPerforming template matching...\n');
    
    % Match template in recorded signal
    fprintf('  Matching in recorded signal...\n');
    [rec_match_pos, rec_correlation] = findTemplateMatch(recorded_signal, template, target_fs);
    
    % Match template in simulated signal
    fprintf('  Matching in simulated signal...\n');
    [sim_match_pos, sim_correlation] = findTemplateMatch(simulated_signal, template, target_fs);
    
    fprintf('  Recorded: Best match at %.2f s (correlation: %.3f)\n', rec_match_pos/target_fs, rec_correlation);
    fprintf('  Simulated: Best match at %.2f s (correlation: %.3f)\n', sim_match_pos/target_fs, sim_correlation);
    
    % === EXTRACT ALIGNED SEGMENTS ===
    fprintf('\nExtracting aligned segments...\n');
    
    % Define extraction parameters
    extraction_samples = template_length;
    
    % Calculate start positions (align template matches)
    rec_start = rec_match_pos;
    sim_start = sim_match_pos;
    
    % Calculate end positions
    rec_end = min(rec_start + extraction_samples - 1, length(recorded_signal));
    sim_end = min(sim_start + extraction_samples - 1, length(simulated_signal));
    
    % Extract segments
    if rec_start > 0 && rec_start <= length(recorded_signal)
        trimmed_recorded = recorded_signal(rec_start:rec_end);
    else
        error('Invalid recorded signal start position: %d', rec_start);
    end
    
    if sim_start > 0 && sim_start <= length(simulated_signal)
        trimmed_simulated = simulated_signal(sim_start:sim_end);
    else
        error('Invalid simulated signal start position: %d', sim_start);
    end
    
    % Ensure both signals have same length
    min_length = min(length(trimmed_recorded), length(trimmed_simulated));
    trimmed_recorded = trimmed_recorded(1:min_length);
    trimmed_simulated = -trimmed_simulated(1:min_length);
    
    % === VERIFY ALIGNMENT ===
    fprintf('\nVerifying final alignment...\n');
    
    if length(trimmed_recorded) > template_length && length(trimmed_simulated) > template_length
        % Compare first part of both trimmed signals
        verify_length = min(template_length, min_length);
        rec_segment = trimmed_recorded(1:verify_length);
        sim_segment = trimmed_simulated(1:verify_length);
        if std(rec_segment) > eps && std(sim_segment) > eps
            verification_corr = corrcoef(rec_segment, sim_segment);
            fprintf('  Verification correlation: %.3f\n', verification_corr(1,2));
            
            if verification_corr(1,2) < 0.3
                fprintf('  ⚠ Warning: Low correlation between aligned signals\n');
            else
                fprintf('  ✓ Good alignment achieved\n');
            end
        else
            fprintf('  Verification correlation: NaN (insufficient variance)\n');
            fprintf('  ⚠ Warning: Cannot verify alignment due to insufficient signal variance\n');
        end
    end
    
    % === SAVE ALIGNED SEGMENTS ===
    fprintf('\nSaving aligned segments...\n');
    
    audiowrite('template_trimmed_recorded.wav', trimmed_recorded, target_fs);
    audiowrite('template_trimmed_simulated.wav', trimmed_simulated, target_fs);
    audiowrite('reference_template.wav', template, target_fs);
    
    fprintf('  ✓ template_trimmed_recorded.wav\n');
    fprintf('  ✓ template_trimmed_simulated.wav\n');
    fprintf('  ✓ reference_template.wav (for reference)\n');
    
    % === TEMPLATE MATCHING SUMMARY ===
    fprintf('\n=== TEMPLATE MATCHING SUMMARY ===\n');
    fprintf('Template matching results:\n');
    fprintf('  Template duration: %.2f s\n', template_duration);
    fprintf('  Recorded match: sample %d (%.2f s), correlation %.3f\n', rec_match_pos, rec_match_pos/target_fs, rec_correlation);
    fprintf('  Simulated match: sample %d (%.2f s), correlation %.3f\n', sim_match_pos, sim_match_pos/target_fs, sim_correlation);
    
    fprintf('\nExtracted segments:\n');
    fprintf('  Duration: %.2f s (%d samples)\n', min_length/target_fs, min_length);
    fprintf('  Sampling rate: %.0f Hz\n', target_fs);
    fprintf('\n');
    
    %% =================================================================
    %% PHASE 2: IMPULSE RESPONSE ANALYSIS
    %% =================================================================
    fprintf('=== PHASE 2: IMPULSE RESPONSE ANALYSIS ===\n\n');
    
    % Use the aligned signals from Phase 1
    measured_signal = trimmed_recorded;
    simulated_signal = trimmed_simulated;
    
    unified_fs = target_fs;
    
    fprintf('Step 1: Using pre-aligned signals from template matching\n');
    fprintf('  Measured: %.2f s, %.0f Hz\n', length(measured_signal)/unified_fs, unified_fs);
    fprintf('  Simulated: %.2f s, %.0f Hz\n', length(simulated_signal)/unified_fs, unified_fs);
    fprintf('\n');
    
    %% =================================================================
    %% STEP 2: GENERATE REFERENCE SIGNAL FOR DECONVOLUTION
    %% =================================================================
    fprintf('Step 2: Preparing reference signal for deconvolution...\n');
    
    % Use the same reference signal parameters
    reference_signal_deconv = generateExpSweep(unified_fs, duration_ref, freqStart, freqEnd);
    
    fprintf('  Reference signal: %.1f s, %d-%d Hz, %.0f Hz sampling\n', ...
        duration_ref, freqStart, freqEnd, unified_fs);
    fprintf('✓ Reference signal prepared\n\n');
    
    %% =================================================================
    %% STEP 3: DECONVOLVE TO EXTRACT IMPULSE RESPONSES
    %% =================================================================
    fprintf('Step 3: Deconvolving aligned signals to extract impulse responses...\n');
    
    % Deconvolution parameters
    reg_value = 0.001;  % Regularization parameter
    signal_choice = 'measurement_sweep';
    use_custom_config = false;
    
    % Deconvolve measured signal
    fprintf('  Deconvolving measured signal...\n');
    measured_IR = deconvolveImpulseResponse(measured_signal, reference_signal_deconv, ...
        unified_fs, signal_choice, use_custom_config, reg_value);
    
    % Deconvolve simulated signal
    fprintf('  Deconvolving simulated signal...\n');
    simulated_IR = deconvolveImpulseResponse(simulated_signal, reference_signal_deconv, ...
        unified_fs, signal_choice, use_custom_config, reg_value);
    
    % Ensure both IRs have same length
    min_ir_length = min(length(measured_IR), length(simulated_IR));
    measured_IR = measured_IR(1:min_ir_length);
    simulated_IR = simulated_IR(1:min_ir_length);
    
    fprintf('✓ Deconvolution completed\n');
    fprintf('  IR length: %.2f s (%d samples)\n', min_ir_length/unified_fs, min_ir_length);
    fprintf('\n');
    
   


    %% =================================================================
    %% STEP 4: COMPREHENSIVE ACOUSTIC ANALYSIS
    %% =================================================================
    fprintf('Step 4: Performing comprehensive acoustic analysis...\n');
    
    % Initialize results structure
    comparison_results = struct();
    comparison_results.metadata = struct();
    comparison_results.metadata.sampling_rate = unified_fs;
    comparison_results.metadata.ir_length_samples = min_ir_length;
    comparison_results.metadata.ir_length_seconds = min_ir_length / unified_fs;
    comparison_results.metadata.analysis_timestamp = datetime('now');
    comparison_results.metadata.reference_signal = reference_signal_deconv;
    comparison_results.metadata.deconvolution_params = struct('reg_value', reg_value, 'signal_choice', signal_choice);
    comparison_results.metadata.template_matching = struct('template_duration', template_duration, ...
        'rec_correlation', rec_correlation, 'sim_correlation', sim_correlation, ...
        'rec_match_time', rec_match_pos/target_fs, 'sim_match_time', sim_match_pos/target_fs);
    
    % Frequency bands for analysis
    freq_bands = [125, 250, 500, 1000, 2000, 4000]; % Standard octave bands
    comparison_results.metadata.frequency_bands = freq_bands;
    
    % === TIME DOMAIN ANALYSIS ===
    fprintf('  Time domain analysis...\n');
    comparison_results.time_domain = analyzeTimeDomain(measured_IR, simulated_IR, unified_fs);
    
    % === FREQUENCY DOMAIN ANALYSIS ===
    fprintf('  Frequency domain analysis...\n');
    comparison_results.frequency_domain = analyzeFrequencyDomain(measured_IR, simulated_IR, unified_fs);
    
    % === ACOUSTIC METRICS ===
    fprintf('  Acoustic metrics analysis...\n');
    comparison_results.acoustic_metrics = analyzeAcousticMetrics(measured_IR, simulated_IR, unified_fs, freq_bands);
    
    % === STATISTICAL ANALYSIS ===
    fprintf('  Statistical analysis...\n');
    comparison_results.statistical_analysis = analyzeStatisticalMetrics(measured_IR, simulated_IR);
    
    fprintf('✓ Analysis completed\n\n');
    
    %% =================================================================
    %% STEP 5: GENERATE COMPREHENSIVE VISUALIZATIONS
    %% =================================================================
    fprintf('Step 5: Generating comparison visualizations...\n');
    
    generateComparisonPlots(measured_IR, simulated_IR, unified_fs, comparison_results);
    
    fprintf('✓ All comparison plots generated\n\n');
    
    %% =================================================================
    %% STEP 6: SAVE RESULTS AND GENERATE REPORT
    %% =================================================================
    fprintf('Step 6: Saving results and generating report...\n');
    
    % Create output directory
    output_dir = 'combined_template_ir_comparison_results';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Save comprehensive results
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    results_filename = fullfile(output_dir, sprintf('Combined_Template_IR_comparison_%s.mat', timestamp_str));
    
    save(results_filename, 'comparison_results', 'measured_IR', 'simulated_IR', ...
         'unified_fs', 'reference_signal_deconv', 'trimmed_recorded', 'trimmed_simulated', ...
         'template', 'rec_match_pos', 'sim_match_pos');
    
    fprintf('✓ Results saved to: %s\n', results_filename);
    
    % Generate console report
    generateConsoleReport(comparison_results);
    
    % Save figures
    saveFigures(output_dir, timestamp_str);
    
    fprintf('\n========================================================\n');
    fprintf('    COMBINED TEMPLATE-ALIGNED ANALYSIS COMPLETE         \n');
    fprintf('========================================================\n');
end

%% =================================================================
%% TEMPLATE MATCHING FUNCTIONS
%% =================================================================

function [match_position, max_correlation] = findTemplateMatch(signal, template, fs)
    % Find best template match using normalized cross-correlation
    
    if length(signal) < length(template)
        error('Signal shorter than template');
    end
    
    % Method 1: FFT-based correlation (fast)
    fprintf('    Using FFT-based template matching...\n');
    
    % Pad template to signal length for FFT
    template_padded = [template; zeros(length(signal) - length(template), 1)];
    
    % FFT-based normalized cross-correlation
    signal_fft = fft(signal);
    template_fft = fft(flipud(template_padded)); % Flip for correlation
    
    % Cross-correlation in frequency domain
    xcorr_fft = ifft(signal_fft .* template_fft);
    xcorr_fft = real(xcorr_fft);
    
    % Normalize correlation
    template_energy = sqrt(sum(template.^2));
    
    % Calculate local energy for normalization
    signal_energy = zeros(length(signal) - length(template) + 1, 1);
    for i = 1:length(signal_energy)
        window = signal(i:i+length(template)-1);
        signal_energy(i) = sqrt(sum(window.^2));
    end
    
    % Normalize correlation values
    valid_length = length(signal_energy);
    % Proper element-wise normalization
    normalized_corr = xcorr_fft(1:valid_length) ./ (signal_energy .* template_energy + eps);
    
    % Find maximum correlation
    [max_correlation, match_position] = max(normalized_corr);
    
    fprintf('    Best correlation: %.3f at sample %d (%.2f s)\n', max_correlation, match_position, match_position/fs);
    
    % Additional verification with direct correlation at peak
    if match_position > 0 && match_position + length(template) - 1 <= length(signal)
        direct_window = signal(match_position:match_position + length(template) - 1);
        if std(template) > eps && std(direct_window) > eps
            direct_corr = corrcoef(template, direct_window);
            fprintf('    Direct verification: %.3f\n', direct_corr(1,2));
        else
            fprintf('    Direct verification: NaN (zero variance)\n');
        end
    end
end

function [signal, fs, found] = loadRecordedSignal()
    % Load recorded signal
    wav_files = dir('acoustic_measurement_results_*/raw_data/*_mic1.wav');
    
    if isempty(wav_files)
        fprintf('  ✗ No recorded mic1 WAV file found\n');
        signal = []; fs = []; found = false;
        return;
    end
    
    wav_path = fullfile(wav_files(end).folder, wav_files(end).name);
    [signal, fs] = audioread(wav_path);
    
    % Convert to mono
    if size(signal, 2) > 1
        signal = mean(signal, 2);
    end
    
    fprintf('  ✓ Loaded recorded: %s (%.2f s, %.0f Hz)\n', wav_files(end).name, length(signal)/fs, fs);
    found = true;
end

function [signal, fs, found] = loadSimulatedSignal()
    % Load simulated signal
    sim_patterns = {'output_data/*_mic_1_signal.wav', 'output_data/*mic1*.wav', '*mic1*signal*.wav'};
    sim_files = [];
    
    for pattern = sim_patterns
        sim_files = dir(pattern{1});
        if ~isempty(sim_files)
            break;
        end
    end
    
    if isempty(sim_files)
        fprintf('  ✗ No simulated mic1 WAV file found\n');
        signal = []; fs = []; found = false;
        return;
    end
    
    sim_path = fullfile(sim_files(end).folder, sim_files(end).name);
    [signal, fs] = audioread(sim_path);
    
    % Convert to mono
    if size(signal, 2) > 1
        signal = mean(signal, 2);
    end
    
    fprintf('  ✓ Loaded simulated: %s (%.2f s, %.0f Hz)\n', sim_files(end).name, length(signal)/fs, fs);
    found = true;
end

%% =================================================================
%% REFERENCE SIGNAL GENERATION
%% =================================================================

function sweepSignal = generateExpSweep(fs, duration, fStart, fEnd)
    % Generate exponential sine sweep (matching simulator)
    N = round(duration * fs);
    t = (0:N-1)' / fs;
    
    % Calculate sweep parameters
    K = (duration * fStart) / log(fEnd/fStart);
    L = duration / log(fEnd/fStart);
    
    % Generate sweep
    sweepSignal = sin(2 * pi * K * (exp(t/L) - 1));
    
    % Apply fade in/out
    fadeLen = round(0.01 * fs); % 10ms fade
    fadeIn = (0:fadeLen-1)' / fadeLen;
    fadeOut = flipud(fadeIn);
    
    sweepSignal(1:fadeLen) = sweepSignal(1:fadeLen) .* fadeIn;
    sweepSignal(end-fadeLen+1:end) = sweepSignal(end-fadeLen+1:end) .* fadeOut;
    
    % Normalize
    sweepSignal = sweepSignal / max(abs(sweepSignal));
end

%% =================================================================
%% DECONVOLUTION FUNCTION
%% =================================================================

function ir = deconvolveImpulseResponse(recordedSignal, referenceSignal, fs, signal_choice, use_custom_config, reg_value)
    % Perform FFT-based deconvolution to extract impulse response
    
    % Ensure signals are column vectors
    recordedSignal = recordedSignal(:);
    referenceSignal = referenceSignal(:);
    
    % Zero-pad to avoid circular convolution artifacts
    N = length(recordedSignal) + length(referenceSignal) - 1;
    N = 2^nextpow2(N); % Use power of 2 for efficient FFT
    
    % Calculate FFTs
    X = fft(referenceSignal, N);
    Y = fft(recordedSignal, N);
    
    reg = reg_value;
    
    H = Y .* conj(X) ./ (X .* conj(X) + reg * max(X .* conj(X)));
    
    % Convert back to time domain
    ir = real(ifft(H));
    
    % Extract useful portion (twice the room's reverberation time or 2 seconds)
    irLength = min(length(ir), round(fs * 2));
    ir = ir(1:irLength);
    
    % Apply fade out to remove artifacts
    fadeLen = round(0.05 * fs); % 50ms fade
    if fadeLen < length(ir)
        fadeOut = ones(size(ir));
        fadeOut(end-fadeLen+1:end) = (fadeLen-1:-1:0)' / fadeLen;
        ir = ir .* fadeOut;
    end
    
    % Normalize
    ir = ir / max(abs(ir));
end

%% =================================================================
%% ANALYSIS FUNCTIONS
%% =================================================================

function results = analyzeTimeDomain(measured, simulated, fs)
    % Comprehensive time domain analysis
    
    results = struct();
    
    % Cross-correlation analysis
    [xcorr_vals, lags] = xcorr(measured, simulated, 'normalized');
    [max_corr, max_idx] = max(abs(xcorr_vals));
    optimal_lag = lags(max_idx);
    time_delay = optimal_lag / fs;
    
    results.cross_correlation.max_correlation = max_corr;
    results.cross_correlation.optimal_lag_samples = optimal_lag;
    results.cross_correlation.time_delay = time_delay;
    results.cross_correlation.values = xcorr_vals;
    results.cross_correlation.lags = lags;
    
    % RMS comparison
    rms_measured = rms(measured);
    rms_simulated = rms(simulated);
    
    results.rms.measured = rms_measured;
    results.rms.simulated = rms_simulated;
    
    if rms_measured > eps
        results.rms.ratio = rms_simulated / rms_measured;
        results.rms.difference_db = 20 * log10(results.rms.ratio);
    else
        results.rms.ratio = Inf;
        results.rms.difference_db = Inf;
    end
    
    % Peak analysis
    results.peak.measured = max(abs(measured));
    results.peak.simulated = max(abs(simulated));
    
    if results.peak.measured > eps
        results.peak.ratio = results.peak.simulated / results.peak.measured;
    else
        results.peak.ratio = Inf;
    end
    
    % Energy decay analysis
    results.energy_decay = analyzeEnergyDecay(measured, simulated, fs);
end

function results = analyzeFrequencyDomain(measured, simulated, fs)
    % Comprehensive frequency domain analysis
    
    results = struct();
    
    % Compute frequency responses
    N = 2^nextpow2(max(length(measured), length(simulated)));
    freqs = (0:N/2) * fs / N;
    
    H_measured = fft(measured, N);
    H_simulated = fft(simulated, N);
    
    % Keep only positive frequencies
    H_measured = H_measured(1:N/2+1);
    H_simulated = H_simulated(1:N/2+1);
    
    results.frequency = freqs;
    results.magnitude.measured_db = 20 * log10(abs(H_measured) + eps);
    results.magnitude.simulated_db = 20 * log10(abs(H_simulated) + eps);
    results.magnitude.difference_db = results.magnitude.simulated_db - results.magnitude.measured_db;
    
    % Phase analysis
    results.phase.measured = angle(H_measured);
    results.phase.simulated = angle(H_simulated);
    results.phase.difference = angle(H_simulated .* conj(H_measured));
    
    % Coherence
    results.coherence = abs(H_measured .* conj(H_simulated)) ./ (abs(H_measured) .* abs(H_simulated) + eps);
    
    % Frequency band analysis
    freq_bands = [20 200; 200 2000; 2000 8000]; % Low, Mid, High
    band_names = {'low', 'mid', 'high'};
    
    for i = 1:size(freq_bands, 1)
        band_indices = freqs >= freq_bands(i,1) & freqs <= freq_bands(i,2);
        
        measured_band = mean(results.magnitude.measured_db(band_indices));
        simulated_band = mean(results.magnitude.simulated_db(band_indices));
        
        results.frequency_bands.(band_names{i}).measured_db = measured_band;
        results.frequency_bands.(band_names{i}).simulated_db = simulated_band;
        results.frequency_bands.(band_names{i}).difference_db = simulated_band - measured_band;
    end
end

function results = analyzeAcousticMetrics(measured, simulated, fs, freq_bands)
    % Advanced acoustic metrics
    
    results = struct();
    
    % RT60 Analysis
    results.rt60 = struct();
    for i = 1:length(freq_bands)
        rt60_measured = calculateRT60Simple(measured, fs, freq_bands(i));
        rt60_simulated = calculateRT60Simple(simulated, fs, freq_bands(i));
        
        band_name = sprintf('f_%dHz', freq_bands(i));
        results.rt60.(band_name).measured = rt60_measured;
        results.rt60.(band_name).simulated = rt60_simulated;
        results.rt60.(band_name).difference = rt60_simulated - rt60_measured;
        if rt60_measured > 0
            results.rt60.(band_name).relative_error = (rt60_simulated - rt60_measured) / rt60_measured * 100;
        else
            results.rt60.(band_name).relative_error = 0;
        end
    end
    
    % Clarity metrics (C50, C80)
    results.clarity = calculateClarityMetrics(measured, simulated, fs);
    
    % Early-to-late energy ratio
    results.early_late_ratio = calculateEarlyLateRatio(measured, simulated, fs);
end

function results = analyzeStatisticalMetrics(measured, simulated)
    % Statistical comparison metrics
    
    results = struct();
    
    % Basic statistics
    corr_matrix = corrcoef(measured, simulated);
    results.correlation_coefficient = corr_matrix(1,2);
    
    % Error metrics
    error_signal = simulated - measured;
    results.mse = mean(error_signal.^2);
    results.rmse = sqrt(results.mse);
    results.mae = mean(abs(error_signal));
    
    % Normalized metrics
    signal_power = mean(measured.^2);
    if signal_power > 0
        results.nmse = results.mse / signal_power;
        results.nrmse = sqrt(results.nmse);
        
        % Signal-to-noise ratio (treating difference as noise)
        noise_power = mean(error_signal.^2);
        results.snr_db = 10 * log10(signal_power / (noise_power + eps));
    else
        results.nmse = Inf;
        results.nrmse = Inf;
        results.snr_db = -Inf;
    end
end

%% =================================================================
%% HELPER FUNCTIONS FOR ACOUSTIC ANALYSIS
%% =================================================================

function rt60 = calculateRT60Simple(ir, fs, center_freq)
    % Simple RT60 calculation for a specific frequency band
    
    % Apply bandpass filter around center frequency
    if center_freq > 0
        % Design simple bandpass filter
        low_freq = center_freq / sqrt(2);
        high_freq = center_freq * sqrt(2);
        
        % Ensure frequencies are within Nyquist limit
        nyquist = fs / 2;
        if high_freq > nyquist
            high_freq = nyquist * 0.9;
        end
        if low_freq < 20
            low_freq = 20;
        end
        
        try
            [b, a] = butter(4, [low_freq high_freq] / nyquist, 'bandpass');
            filtered_ir = filter(b, a, ir);
        catch
            filtered_ir = ir; % Fallback to unfiltered
        end
    else
        filtered_ir = ir;
    end
    
    % Calculate energy decay curve
    energy = filtered_ir.^2;
    
    % Reverse integrate (Schroeder integration)
    edc = cumsum(energy(end:-1:1));
    edc = edc(end:-1:1);
    
    % Convert to dB
    edc_db = 10 * log10(edc / max(edc) + eps);
    
    % Find -5dB and -35dB points
    idx_5db = find(edc_db <= -5, 1, 'first');
    idx_35db = find(edc_db <= -35, 1, 'first');
    
    if isempty(idx_5db) || isempty(idx_35db) || idx_35db <= idx_5db
        rt60 = 0.1; % Default value
    else
        % Calculate T30 and extrapolate to RT60
        t30 = (idx_35db - idx_5db) / fs;
        rt60 = 2 * t30; % RT60 = 2 * T30
    end
    
    % Clamp to reasonable values
    rt60 = max(0.01, min(rt60, 10));
end

function results = calculateClarityMetrics(measured, simulated, fs)
    % Calculate C50 and C80 metrics
    
    results = struct();
    
    % Find direct sound peak
    [~, peak_idx_m] = max(abs(measured));
    [~, peak_idx_s] = max(abs(simulated));
    
    % Calculate C50 (50ms clarity)
    ms50_samples = round(0.05 * fs);
    
    % Measured
    early_end_m = min(length(measured), peak_idx_m + ms50_samples);
    late_start_m = min(length(measured), peak_idx_m + ms50_samples + 1);
    
    early_energy_m = sum(measured(peak_idx_m:early_end_m).^2);
    late_energy_m = sum(measured(late_start_m:end).^2);
    c50_measured = 10 * log10(early_energy_m / (late_energy_m + eps));
    
    % Simulated  
    early_end_s = min(length(simulated), peak_idx_s + ms50_samples);
    late_start_s = min(length(simulated), peak_idx_s + ms50_samples + 1);
    
    early_energy_s = sum(simulated(peak_idx_s:early_end_s).^2);
    late_energy_s = sum(simulated(late_start_s:end).^2);
    c50_simulated = 10 * log10(early_energy_s / (late_energy_s + eps));
    
    results.c50.measured = c50_measured;
    results.c50.simulated = c50_simulated;
    results.c50.difference = c50_simulated - c50_measured;
    
    % Calculate C80 (80ms clarity)
    ms80_samples = round(0.08 * fs);
    
    % Measured
    early_end_m = min(length(measured), peak_idx_m + ms80_samples);
    late_start_m = min(length(measured), peak_idx_m + ms80_samples + 1);
    
    early_energy_m = sum(measured(peak_idx_m:early_end_m).^2);
    late_energy_m = sum(measured(late_start_m:end).^2);
    c80_measured = 10 * log10(early_energy_m / (late_energy_m + eps));
    
    % Simulated
    early_end_s = min(length(simulated), peak_idx_s + ms80_samples);
    late_start_s = min(length(simulated), peak_idx_s + ms80_samples + 1);
    
    early_energy_s = sum(simulated(peak_idx_s:early_end_s).^2);
    late_energy_s = sum(simulated(late_start_s:end).^2);
    c80_simulated = 10 * log10(early_energy_s / (late_energy_s + eps));
    
    results.c80.measured = c80_measured;
    results.c80.simulated = c80_simulated;
    results.c80.difference = c80_simulated - c80_measured;
end

function results = calculateEarlyLateRatio(measured, simulated, fs)
    % Calculate early-to-late energy ratio
    
    results = struct();
    
    % Define early period (first 50ms after peak)
    [~, peak_idx_m] = max(abs(measured));
    [~, peak_idx_s] = max(abs(simulated));
    
    early_samples = round(0.05 * fs); % 50ms
    
    % Measured
    early_end_m = min(length(measured), peak_idx_m + early_samples);
    early_energy_m = sum(measured(peak_idx_m:early_end_m).^2);
    total_energy_m = sum(measured.^2);
    ratio_measured = early_energy_m / total_energy_m;
    
    % Simulated
    early_end_s = min(length(simulated), peak_idx_s + early_samples);
    early_energy_s = sum(simulated(peak_idx_s:early_end_s).^2);
    total_energy_s = sum(simulated.^2);
    ratio_simulated = early_energy_s / total_energy_s;
    
    results.measured = ratio_measured;
    results.simulated = ratio_simulated;
    results.difference = ratio_simulated - ratio_measured;
end

function results = analyzeEnergyDecay(measured, simulated, fs)
    % Analyze energy decay characteristics
    
    results = struct();
    
    % Calculate energy decay curves
    energy_m = measured.^2;
    energy_s = simulated.^2;
    
    % Reverse cumulative sum (Schroeder integration)
    edc_m = cumsum(energy_m(end:-1:1));
    edc_m = edc_m(end:-1:1);
    edc_s = cumsum(energy_s(end:-1:1));
    edc_s = edc_s(end:-1:1);
    
    % Normalize and convert to dB
    edc_m_db = 10 * log10(edc_m / max(edc_m) + eps);
    edc_s_db = 10 * log10(edc_s / max(edc_s) + eps);
    
    % Time vector
    time_vec = (0:length(edc_m_db)-1) / fs;
    
    results.time = time_vec;
    results.measured_db = edc_m_db;
    results.simulated_db = edc_s_db;
    results.difference_db = edc_s_db - edc_m_db;
    
    % Calculate decay rates in different ranges
    % Early decay (0-10dB)
    idx_0db = 1;
    idx_10db_m = find(edc_m_db <= -10, 1, 'first');
    idx_10db_s = find(edc_s_db <= -10, 1, 'first');
    
    if ~isempty(idx_10db_m) && idx_10db_m > idx_0db
        early_decay_rate_m = -10 / (time_vec(idx_10db_m) - time_vec(idx_0db));
    else
        early_decay_rate_m = NaN;
    end
    
    if ~isempty(idx_10db_s) && idx_10db_s > idx_0db
        early_decay_rate_s = -10 / (time_vec(idx_10db_s) - time_vec(idx_0db));
    else
        early_decay_rate_s = NaN;
    end
    
    results.early_decay_rate.measured = early_decay_rate_m;
    results.early_decay_rate.simulated = early_decay_rate_s;
end

%% =================================================================
%% VISUALIZATION FUNCTIONS
%% =================================================================

function generateComparisonPlots(measured_IR, simulated_IR, fs, results)
    % Generate comprehensive comparison plots
    
    % Time vector
    time_vec = (0:length(measured_IR)-1) / fs;
    
    % Figure 1: Time Domain Comparison
    figure('Position', [100, 100, 1200, 600], 'Name', 'Time Domain Comparison');
    
    subplot(2, 2, 1);
    plot(time_vec, measured_IR, 'b-', 'LineWidth', 1.5); hold on;
    plot(time_vec, simulated_IR, 'r--', 'LineWidth', 1.5);
    title('Impulse Response Comparison');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Measured', 'Simulated', 'Location', 'best');
    grid on;
    
    % Cross-correlation subplot
    subplot(2, 2, 2);
    plot(results.time_domain.cross_correlation.lags / fs, results.time_domain.cross_correlation.values, 'g-');
    title(sprintf('Cross Correlation (Max: %.3f)', results.time_domain.cross_correlation.max_correlation));
    xlabel('Lag (s)');
    ylabel('Correlation');
    grid on;
    
    % Energy decay curves
    subplot(2, 2, 3);
    if isfield(results.time_domain, 'energy_decay')
        plot(results.time_domain.energy_decay.time, results.time_domain.energy_decay.measured_db, 'b-', 'LineWidth', 2); hold on;
        plot(results.time_domain.energy_decay.time, results.time_domain.energy_decay.simulated_db, 'r--', 'LineWidth', 2);
        title('Energy Decay Curves');
        xlabel('Time (s)');
        ylabel('Level (dB)');
        legend('Measured', 'Simulated', 'Location', 'best');
        grid on;
        ylim([-60 0]);
    end
    
    % RMS and peak comparison
    subplot(2, 2, 4);
    metrics = [results.time_domain.rms.measured, results.time_domain.rms.simulated; ...
               results.time_domain.peak.measured, results.time_domain.peak.simulated];
    bar(metrics);
    title('RMS and Peak Comparison');
    set(gca, 'XTickLabel', {'RMS', 'Peak'});
    ylabel('Amplitude');
    legend('Measured', 'Simulated', 'Location', 'best');
    grid on;
    
    sgtitle('Time Domain Analysis', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Figure 2: Frequency Domain Comparison
    figure('Position', [150, 150, 1200, 600], 'Name', 'Frequency Domain Comparison');
    
    subplot(2, 2, 1);
    semilogx(results.frequency_domain.frequency, results.frequency_domain.magnitude.measured_db, 'b-', 'LineWidth', 1.5); hold on;
    semilogx(results.frequency_domain.frequency, results.frequency_domain.magnitude.simulated_db, 'r--', 'LineWidth', 1.5);
    title('Frequency Response');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    legend('Measured', 'Simulated', 'Location', 'best');
    grid on;
    xlim([20 fs/2]);
    
    subplot(2, 2, 2);
    semilogx(results.frequency_domain.frequency, results.frequency_domain.magnitude.difference_db, 'k-', 'LineWidth', 1.5);
    title('Magnitude Difference');
    xlabel('Frequency (Hz)');
    ylabel('Difference (dB)');
    grid on;
    xlim([20 fs/2]);
    
    subplot(2, 2, 3);
    semilogx(results.frequency_domain.frequency, results.frequency_domain.phase.measured, 'b-', 'LineWidth', 1.5); hold on;
    semilogx(results.frequency_domain.frequency, results.frequency_domain.phase.simulated, 'r--', 'LineWidth', 1.5);
    title('Phase Response');
    xlabel('Frequency (Hz)');
    ylabel('Phase (rad)');
    legend('Measured', 'Simulated', 'Location', 'best');
    grid on;
    xlim([20 fs/2]);
    
    subplot(2, 2, 4);
    semilogx(results.frequency_domain.frequency, results.frequency_domain.coherence, 'g-', 'LineWidth', 1.5);
    title('Coherence');
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    grid on;
    xlim([20 fs/2]);
    ylim([0 1]);
    
    sgtitle('Frequency Domain Analysis', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Figure 3: Acoustic Metrics Summary
    figure('Position', [200, 200, 1000, 600], 'Name', 'Acoustic Metrics Comparison');
    
    % RT60 comparison
    subplot(2, 2, 1);
    freq_bands = results.metadata.frequency_bands;
    rt60_measured = zeros(1, length(freq_bands));
    rt60_simulated = zeros(1, length(freq_bands));
    
    for i = 1:length(freq_bands)
        band_name = sprintf('f_%dHz', freq_bands(i));
        if isfield(results.acoustic_metrics.rt60, band_name)
            rt60_measured(i) = results.acoustic_metrics.rt60.(band_name).measured;
            rt60_simulated(i) = results.acoustic_metrics.rt60.(band_name).simulated;
        end
    end
    
    semilogx(freq_bands, rt60_measured, 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
    semilogx(freq_bands, rt60_simulated, 'rs--', 'LineWidth', 2, 'MarkerSize', 8);
    title('RT60 Comparison');
    xlabel('Frequency (Hz)');
    ylabel('RT60 (s)');
    legend('Measured', 'Simulated', 'Location', 'best');
    grid on;
    
    % Correlation and error metrics
    subplot(2, 2, 2);
    corr_coeff = results.statistical_analysis.correlation_coefficient;
    nrmse = results.statistical_analysis.nrmse;
    snr_db = results.statistical_analysis.snr_db;
    
    metrics_data = [corr_coeff; nrmse; snr_db/10]; % Scale SNR for visualization
    bar(metrics_data, 'FaceColor', [0.3 0.7 0.9]);
    title('Statistical Metrics');
    set(gca, 'XTickLabel', {'Correlation', 'NRMSE', 'SNR/10'});
    ylabel('Value');
    grid on;
    
    % Clarity metrics
    subplot(2, 2, 3);
    if isfield(results.acoustic_metrics, 'clarity')
        c50_data = [results.acoustic_metrics.clarity.c50.measured, results.acoustic_metrics.clarity.c50.simulated];
        c80_data = [results.acoustic_metrics.clarity.c80.measured, results.acoustic_metrics.clarity.c80.simulated];
        
        clarity_data = [c50_data; c80_data];
        bar(clarity_data);
        title('Clarity Metrics');
        set(gca, 'XTickLabel', {'C50', 'C80'});
        ylabel('dB');
        legend('Measured', 'Simulated', 'Location', 'best');
        grid on;
    end
    
    % Frequency band comparison
    subplot(2, 2, 4);
    band_names = {'Low', 'Mid', 'High'};
    measured_bands = [results.frequency_domain.frequency_bands.low.measured_db, ...
                      results.frequency_domain.frequency_bands.mid.measured_db, ...
                      results.frequency_domain.frequency_bands.high.measured_db];
    simulated_bands = [results.frequency_domain.frequency_bands.low.simulated_db, ...
                       results.frequency_domain.frequency_bands.mid.simulated_db, ...
                       results.frequency_domain.frequency_bands.high.simulated_db];
    
    band_data = [measured_bands; simulated_bands];
    bar(band_data);
    title('Frequency Band Levels');
    set(gca, 'XTickLabel', band_names);
    ylabel('Level (dB)');
    legend('Measured', 'Simulated', 'Location', 'best');
    grid on;
    
    sgtitle('Acoustic Metrics Summary', 'FontSize', 16, 'FontWeight', 'bold');
end

function saveFigures(output_dir, timestamp_str)
    % Save all open figures
    fig_handles = findall(0, 'Type', 'figure');
    
    for i = 1:length(fig_handles)
        fig = fig_handles(i);
        if ~isempty(fig.Name)
            filename = sprintf('%s_%s.png', strrep(fig.Name, ' ', '_'), timestamp_str);
            filepath = fullfile(output_dir, filename);
            saveas(fig, filepath);
            fprintf('  ✓ Saved figure: %s\n', filename);
        end
    end
end

%% =================================================================
%% REPORTING FUNCTIONS
%% =================================================================

function generateConsoleReport(results)
    % Generate comprehensive console report
    
    fprintf('\n========================================================\n');
    fprintf('         COMBINED TEMPLATE-ALIGNED IR COMPARISON         \n');
    fprintf('========================================================\n\n');
    
    % Template matching summary
    if isfield(results.metadata, 'template_matching')
        tm = results.metadata.template_matching;
        fprintf('TEMPLATE MATCHING RESULTS:\n');
        fprintf('  Template duration: %.2f s\n', tm.template_duration);
        fprintf('  Recorded correlation: %.3f at %.2f s\n', tm.rec_correlation, tm.rec_match_time);
        fprintf('  Simulated correlation: %.3f at %.2f s\n', tm.sim_correlation, tm.sim_match_time);
        fprintf('\n');
    end
    
    % Overall summary
    fprintf('ANALYSIS SUMMARY:\n');
    fprintf('  Sampling rate: %.0f Hz\n', results.metadata.sampling_rate);
    fprintf('  IR duration: %.2f seconds (%d samples)\n', results.metadata.ir_length_seconds, results.metadata.ir_length_samples);
    fprintf('  Analysis timestamp: %s\n', char(results.metadata.analysis_timestamp));
    fprintf('\n');
    
    % Statistical analysis
    fprintf('STATISTICAL ANALYSIS:\n');
    fprintf('  Correlation coefficient: %.3f\n', results.statistical_analysis.correlation_coefficient);
    fprintf('  NRMSE: %.4f\n', results.statistical_analysis.nrmse);
    fprintf('  SNR: %.1f dB\n', results.statistical_analysis.snr_db);
    fprintf('  RMSE: %.4f\n', results.statistical_analysis.rmse);
    fprintf('  MAE: %.4f\n', results.statistical_analysis.mae);
    fprintf('\n');
    
    % Quality assessment
    corr = results.statistical_analysis.correlation_coefficient;
    nrmse = results.statistical_analysis.nrmse;
    
    fprintf('QUALITY ASSESSMENT:\n');
    if corr > 0.9 && nrmse < 0.1
        quality = 'EXCELLENT';
        quality_symbol = '✓';
    elseif corr > 0.8 && nrmse < 0.2
        quality = 'GOOD';
        quality_symbol = '✓';
    elseif corr > 0.6 && nrmse < 0.4
        quality = 'FAIR';
        quality_symbol = '○';
    else
        quality = 'POOR';
        quality_symbol = '✗';
    end
    fprintf('  Overall simulation quality: %s %s\n', quality_symbol, quality);
    fprintf('\n');
    
    % Time domain analysis
    fprintf('TIME DOMAIN ANALYSIS:\n');
    fprintf('  Max cross-correlation: %.3f\n', results.time_domain.cross_correlation.max_correlation);
    fprintf('  Optimal time delay: %.3f ms\n', results.time_domain.cross_correlation.time_delay * 1000);
    fprintf('  RMS ratio (sim/meas): %.3f (%.1f dB)\n', results.time_domain.rms.ratio, results.time_domain.rms.difference_db);
    fprintf('  Peak ratio (sim/meas): %.3f\n', results.time_domain.peak.ratio);
    fprintf('\n');
    
    % Frequency domain analysis
    fprintf('FREQUENCY DOMAIN ANALYSIS:\n');
    fprintf('  Low band (20-200 Hz): %.1f dB difference\n', results.frequency_domain.frequency_bands.low.difference_db);
    fprintf('  Mid band (200-2000 Hz): %.1f dB difference\n', results.frequency_domain.frequency_bands.mid.difference_db);
    fprintf('  High band (2000-8000 Hz): %.1f dB difference\n', results.frequency_domain.frequency_bands.high.difference_db);
    fprintf('\n');
    
    % RT60 analysis
    fprintf('RT60 ANALYSIS:\n');
    freq_bands = results.metadata.frequency_bands;
    for i = 1:length(freq_bands)
        band_name = sprintf('f_%dHz', freq_bands(i));
        if isfield(results.acoustic_metrics.rt60, band_name)
            rt60_data = results.acoustic_metrics.rt60.(band_name);
            fprintf('  %4d Hz: Measured=%.2fs, Simulated=%.2fs, Error=%+.1f%%\n', ...
                freq_bands(i), rt60_data.measured, rt60_data.simulated, rt60_data.relative_error);
        end
    end
    fprintf('\n');
    
    % Clarity metrics
    if isfield(results.acoustic_metrics, 'clarity')
        fprintf('CLARITY METRICS:\n');
        fprintf('  C50: Measured=%.1f dB, Simulated=%.1f dB, Difference=%.1f dB\n', ...
            results.acoustic_metrics.clarity.c50.measured, ...
            results.acoustic_metrics.clarity.c50.simulated, ...
            results.acoustic_metrics.clarity.c50.difference);
        fprintf('  C80: Measured=%.1f dB, Simulated=%.1f dB, Difference=%.1f dB\n', ...
            results.acoustic_metrics.clarity.c80.measured, ...
            results.acoustic_metrics.clarity.c80.simulated, ...
            results.acoustic_metrics.clarity.c80.difference);
        fprintf('\n');
    end
    
    % Early-to-late ratio
    if isfield(results.acoustic_metrics, 'early_late_ratio')
        fprintf('EARLY-TO-LATE ENERGY RATIO:\n');
        fprintf('  Measured: %.3f\n', results.acoustic_metrics.early_late_ratio.measured);
        fprintf('  Simulated: %.3f\n', results.acoustic_metrics.early_late_ratio.simulated);
        fprintf('  Difference: %+.3f\n', results.acoustic_metrics.early_late_ratio.difference);
        fprintf('\n');
    end
    
    % Recommendations
    fprintf('RECOMMENDATIONS:\n');
    if corr < 0.7
        fprintf('  • Low correlation suggests significant waveform differences\n');
    end
    if nrmse > 0.3
        fprintf('  • High NRMSE indicates substantial prediction errors\n');
    end
    if abs(results.time_domain.rms.difference_db) > 3
        fprintf('  • RMS level difference >3dB may indicate gain mismatch\n');
    end
    
    % Check RT60 errors
    rt60_errors = [];
    for i = 1:length(freq_bands)
        band_name = sprintf('f_%dHz', freq_bands(i));
        if isfield(results.acoustic_metrics.rt60, band_name)
            rt60_errors(end+1) = abs(results.acoustic_metrics.rt60.(band_name).relative_error);
        end
    end
    if ~isempty(rt60_errors) && mean(rt60_errors) > 20
        fprintf('  • Average RT60 error >20%% suggests reverberation modeling issues\n');
    end
    
    if corr > 0.8 && nrmse < 0.2
        fprintf('  • Good overall agreement between measured and simulated data\n');
    end
    
    fprintf('\n');
    fprintf('SUCCESS: Combined template-aligned analysis completed!\n');
    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('Combined template-aligned IR comparison completed successfully\n');
    fprintf('========================================================\n');
end

