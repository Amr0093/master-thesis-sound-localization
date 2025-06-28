% SC-DAMAS Implementation for Acoustic Imaging in Reverberant Field
% This script implements the SC-DAMAS algorithm described in the paper
% "An improved deconvolution beamforming algorithm for acoustic imaging
% of low signal-to-noise ratio sound sources in reverberant field"

%% Check for external call and prepare
if exist('EXTERNAL_CALL_FLAG', 'var') && EXTERNAL_CALL_FLAG
    fprintf('Beamforming: Received external call, cleaning workspace...\n');
    % clear;
    close all;
    % clc;
    % Send "ready" signal by throwing a specific error
    error('Beamforming_check_completed'); % This is our "handshake" signal
end

%% Normal operation - check if external variables exist
if exist('final_signals', 'var') || exist('batch_mode', 'var')
    % External call (second time) - external variables are set, don't clear
    fprintf('Beamforming is Ready for normal external usage\n');
else
    % Independent run - clear workspace
    clear;
    close all;
    clc;
    fprintf('Beamforming is Ready for internal usage\n');
end
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

%% Get simulated microphone data from simulator workspace
% Variables should already exist from FastHybridDopplerReverbSimulator_enhanced_2
if ~exist('final_signals', 'var')
    error('Beamforming requires simulator data. Run FastHybridDopplerReverbSimulator_enhanced_2 first.');
end

rng(42);

% Map simulator variables to beamforming expected names
noisy_signals = final_signals;  % Main signals from simulator
numMics = size(final_signals, 1);
t = (0:size(final_signals,2)-1)/fs;  % Time vector

% Convert impulse responses to frequency domain for H
if exist('simulated_impulse_responses', 'var')
    H = fft(simulated_impulse_responses, [], 2);
    freq_vector = linspace(0, fs/2, size(H, 2)/2+1);
    H = H(:, 1:length(freq_vector)); % Take positive frequencies only
else
    % Fallback: create basic frequency response
    freq_vector = linspace(0, fs/2, floor(size(final_signals,2)/2)+1);
    H = ones(numMics, length(freq_vector));
    warning('No impulse responses found, using simplified frequency response');
end

% sourcePos and micPositions should already exist from simulator
% roomDim should already exist from simulator
% SNR_dB should already exist from simulator

%% Parameters
c = 343;                      % Speed of sound (m/s)
% freq_to_analyze = [1000, 2000, 4000];  % Frequencies to analyze (Hz)
freq_to_analyze = 1000;
grid_size = 15;               % Size of scanning grid (odd number preferred)
scan_area_size = 1;           % Size of scanning area in meters
sparsity_step = 3;            % Step size for adaptive sparsity (parameter s in SAMP-CoSaMP)
delta_criteria = 1e-3;        % Relaxed convergence criteria (was 1e-6)
max_iterations = 20;          % Reduced max iterations (was 50)
plot_results = true;          % Set to true to plot results
save_results = true;          % Set to true to save results

%% Create scanning grid
% Grid centered at actual source position for testing/validation
grid_center = [sourcePos(1), sourcePos(3)];  % x, z coordinates

% Define scanning grid (xz-plane at source's y position)
% Define scanning grid (xz-plane at source's y position)
[x_grid, z_grid] = meshgrid(linspace(grid_center(1)-scan_area_size/2, grid_center(1)+scan_area_size/2, grid_size), ...
                            linspace(grid_center(2)-scan_area_size/2, grid_center(2)+scan_area_size/2, grid_size));
grid_points = [x_grid(:), repmat(sourcePos(2), grid_size^2, 1), z_grid(:)];

% Safety check: ensure grid points aren't too close to microphones
min_distance = 0.05; % 5cm minimum
for i = 1:size(grid_points,1)
    for j = 1:numMics
        if norm(grid_points(i,:) - micPositions(j,:)) < min_distance
            warning('Grid point %d too close to microphone %d, results may be unstable', i, j);
        end
    end
end

disp(['Created scanning grid with ', num2str(grid_size), 'x', num2str(grid_size), ' = ', num2str(grid_size^2), ' points']);
fprintf('Processing %d grid points...\n', grid_size^2);

%% Initialize results containers
cbf_results = zeros(length(freq_to_analyze), grid_size, grid_size);
omp_damas_results = zeros(length(freq_to_analyze), grid_size, grid_size);
sc_damas_results = zeros(length(freq_to_analyze), grid_size, grid_size);

% Loop through frequencies to analyze
for freq_idx = 1:length(freq_to_analyze)
    freq = freq_to_analyze(freq_idx);
    disp(['Processing frequency: ', num2str(freq), ' Hz']);
    
    % Find closest frequency bin index
    [~, freq_bin_idx] = min(abs(freq_vector - freq));
    
    %% Extract signals at target frequency (using FFT)
    % Get FFT of microphone signals
    signal_length = size(noisy_signals, 2);
    mic_fft = fft(noisy_signals, signal_length, 2);
    
    % Extract the complex pressure at the target frequency
    p = mic_fft(:, freq_bin_idx);
    
    %% Calculate Cross-Spectral Matrix (CSM)
    CSM = p * p';
    
    % Remove diagonal for noise reduction (optional)
    CSM = CSM - diag(diag(CSM));
    
    %% Conventional Beamforming (CBF) - for comparison
    disp('Calculating Conventional Beamforming (CBF)...');
    b_cbf = zeros(grid_size^2, 1);
    
    % For each scanning point
    for i = 1:grid_size^2
        grid_point = grid_points(i, :);
        
        % Calculate steering vector using room impulse response (not free-field Green's function)
        % This is the key difference for reverberant environments
        a = zeros(numMics, 1);
        
        for mic = 1:numMics
            % Extract frequency response H for this microphone at the target frequency
            h_mic_freq = H(mic, freq_bin_idx);
            
            % Calculate distance from grid point to microphone
            r = norm(grid_point - micPositions(mic, :));
            
            % Phase adjustment for steering
            phase_factor = exp(-1j * 2 * pi * freq * r / c);
            
            % Combine with room impulse response frequency component
            a(mic) = h_mic_freq * phase_factor;
        end
        
        % Normalize steering vector
        a = a / norm(a);
        
        % Calculate beamformer output
        b_cbf(i) = abs(a' * CSM * a);
    end
    
    % Reshape and store CBF results
    cbf_results(freq_idx, :, :) = reshape(b_cbf, grid_size, grid_size);
    
    %% OMP-DAMAS (for comparison)
    disp('Calculating OMP-DAMAS...');
    
    % Initialize PSF (Point Spread Function) matrix
    PSF = zeros(grid_size^2, grid_size^2);
    
    % Calculate PSF
    for i = 1:grid_size^2
        % Steering vector for point i
        a_i = zeros(numMics, 1);
        grid_point_i = grid_points(i, :);
        
        for mic = 1:numMics
            h_mic_freq = H(mic, freq_bin_idx);
            r = norm(grid_point_i - micPositions(mic, :));
            phase_factor = exp(-1j * 2 * pi * freq * r / c);
            a_i(mic) = h_mic_freq * phase_factor;
        end
        a_i = a_i / norm(a_i);
        
        for j = 1:grid_size^2
            % Steering vector for point j
            a_j = zeros(numMics, 1);
            grid_point_j = grid_points(j, :);
            
            for mic = 1:numMics
                h_mic_freq = H(mic, freq_bin_idx);
                r = norm(grid_point_j - micPositions(mic, :));
                phase_factor = exp(-1j * 2 * pi * freq * r / c);
                a_j(mic) = h_mic_freq * phase_factor;
            end
            a_j = a_j / norm(a_j);
            
            % PSF(i,j) is the spatial response of focusing at point i when source is at point j
            PSF(i, j) = abs(a_i' * a_j)^2;
        end
    end
    
    % OMP-DAMAS algorithm
    x_omp = zeros(grid_size^2, 1);
    r = b_cbf;  % Initial residual is the CBF result
    max_iter = 20;  % Maximum iterations for OMP
    
    for iter = 1:max_iter
        % Find index with maximum correlation
        [~, idx] = max(abs(PSF' * r));
        
        % Update solution
        x_omp(idx) = x_omp(idx) + r(idx);
        
        % Update residual
        r = b_cbf - PSF * x_omp;
        
        % Check convergence
        if norm(r) < delta_criteria
            break;
        end
    end
    
    % Reshape and store OMP-DAMAS results
    omp_damas_results(freq_idx, :, :) = reshape(x_omp, grid_size, grid_size);
    
    %% SC-DAMAS algorithm (SAMP-CoSaMP based DAMAS)
    disp('Calculating SC-DAMAS...');
    
    % Initialize SC-DAMAS variables
    %x = zeros(grid_size^2, 1);  % Solution vector
    r = b_cbf;                  % Initial residual
    stage = 1;                  % Initial stage for sparsity adaptation
    L = sparsity_step;          % Initial sparsity estimate (step size)
    F = [];                     % Support set
    k = 1;                      % Iteration counter
    
    % Step 1: Initial sparsity estimation via SAMP
    % Calculate initial correlation
    u = abs(PSF' * r);
    [~, indices] = sort(u, 'descend');
    K = 1;  % Start with minimal sparsity
    
    % Find initial support set
    while true
        F = indices(1:K);
        
        % Calculate projection
        PSF_F = PSF(:, F);
        x_F = zeros(grid_size^2, 1);
        % Simple singularity protection
       % Improved singularity protection
        if cond(PSF_F) > 1e8 || size(PSF_F,1) < size(PSF_F,2)
            x_F(F) = (PSF_F' * PSF_F + 1e-4 * eye(size(PSF_F,2))) \ (PSF_F' * b_cbf);
        else
            x_F(F) = PSF_F \ b_cbf;  % Use backslash instead of pinv
        end
        
        % Check if we need to increase sparsity
        residual_norm = norm(b_cbf - PSF * x_F);
        
        % Compute delta_K (we'll use a simplified version)
        delta_K = 0.1;
        
        if residual_norm <= ((1-delta_K)/sqrt(1+delta_K)) * norm(b_cbf)
            break;
        else
            K = K + L;
            if K > grid_size^2
                K = grid_size^2;
                break;
            end
        end
    end
    
    % Initial residual
    r = b_cbf - PSF * x_F;
    x = x_F;
    
    % Step 2: Main SC-DAMAS iterations (combines SAMP and CoSaMP)
    % Step 2: Main SC-DAMAS iterations (combines SAMP and CoSaMP)
    tic; % Start timer
    max_time = 60; % 60 second timeout
    
    while k <= max_iterations
        if toc > max_time
            warning('SC-DAMAS timeout reached, stopping iterations');
            break;
        end
        
        current_residual = norm(r);
        if mod(k, 5) == 0
            fprintf('SC-DAMAS iteration %d/%d, residual norm: %.6f\n', k, max_iterations, current_residual);
        end
        
        % Check for stagnation
        if k > 10 && exist('prev_residual', 'var')
            improvement_rate = (prev_residual - current_residual) / prev_residual;
            if improvement_rate < 0.01  % Less than 1% improvement
                fprintf('SC-DAMAS stagnating (%.1f%% improvement), stopping early\n', improvement_rate*100);
                break;
            end
        end
        prev_residual = current_residual;
        
        % Calculate correlation
        u = abs(PSF' * r);
        
        % Identify significant components (using a threshold)
        threshold = max(0.1 * max(u), 0.5 * max(u) * exp(-k/10));
        candidate_set = find(u >= threshold);
        
        % Choose the top K components
        [~, indices] = sort(u(candidate_set), 'descend');
        top_indices = candidate_set(indices(1:min(K, length(indices))));
        
        % Update support set
        F_new = union(F, top_indices');
        
        % Calculate new estimate
        PSF_F_new = PSF(:, F_new);
        x_new = zeros(grid_size^2, 1);
        % Simple singularity protection
        % Improved singularity protection
        if cond(PSF_F_new) > 1e8 || size(PSF_F_new,1) < size(PSF_F_new,2)
            x_new(F_new) = (PSF_F_new' * PSF_F_new + 1e-4 * eye(size(PSF_F_new,2))) \ (PSF_F_new' * b_cbf);
        else
            x_new(F_new) = PSF_F_new \ b_cbf;  % Use backslash instead of pinv
        end
        
        % Calculate new residual
        r_new = b_cbf - PSF * x_new;
        
        % Check convergence
        convergence_norm = norm(x_new - x);
        if convergence_norm < delta_criteria
            fprintf('SC-DAMAS converged at iteration %d (norm: %.2e)\n', k, convergence_norm);
            x = x_new;
            break;
        end
        
        % Additional safety: check for NaN/Inf
        if any(~isfinite(x_new))
            warning('SC-DAMAS produced non-finite values, stopping iterations');
            break;
        end
        
        % Check if we need to increase sparsity
        
        if norm(r_new) >= norm(r)
            stage = stage + 1;
            L = stage * sparsity_step;
            K = min(K + L, grid_size^2);
            
            % If sparsity gets too high, just exit
            if K > grid_size^2 * 0.8
                fprintf('SC-DAMAS: High sparsity reached, stopping\n');
                break;
            end
        else
            F = F_new;
            r = r_new;
            x = x_new;
            k = k + 1;
        end
    end
    
    % Ensure non-negativity (sound sources should have positive power)
    x(x < 0) = 0;
    
    % Reshape and store SC-DAMAS results
    sc_damas_results(freq_idx, :, :) = reshape(x, grid_size, grid_size);
end

%% Plot results for comparison with paper figures
if plot_results
    % Create figure similar to Figs. 4-6 in the paper
    for freq_idx = 1:length(freq_to_analyze)
        freq = freq_to_analyze(freq_idx);
        
        figure('Position', [100, 100, 1200, 400]);
        
        % Plot CBF results
        subplot(1, 3, 1);
        imagesc(squeeze(cbf_results(freq_idx, :, :)));
        title(['CBF at ', num2str(freq), ' Hz']);
        colorbar;
        axis equal tight;
        xlabel('X grid index');
        ylabel('Z grid index');
        
        % Plot OMP-DAMAS results
        subplot(1, 3, 2);
        imagesc(squeeze(omp_damas_results(freq_idx, :, :)));
        title(['OMP-DAMAS at ', num2str(freq), ' Hz']);
        colorbar;
        axis equal tight;
        xlabel('X grid index'); 
        ylabel('Z grid index');
        
        % Plot SC-DAMAS results
        subplot(1, 3, 3);
        imagesc(squeeze(sc_damas_results(freq_idx, :, :)));
        title(['SC-DAMAS at ', num2str(freq), ' Hz']);
        colorbar;
        axis equal tight;
        xlabel('X grid index');
        ylabel('Z grid index');
        
        sgtitle(['Acoustic Imaging Results at ', num2str(freq), ' Hz, SNR = ', num2str(SNR_dB), ' dB']);
        % Save 2D figure immediately after creation
        if save_results
            if ~exist('figures_dir', 'var')
                figures_dir = fullfile('beamforming_results', 'figures');
                if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
            end
            figure_name = sprintf('beamforming_2D_%dHz_%s', freq, timestamp);
            saveas(gcf, fullfile(figures_dir, [figure_name '.png']), 'png');
        end
    end
    
    % Plot 3D surface plots for better visualization
    for freq_idx = 1:length(freq_to_analyze)
        freq = freq_to_analyze(freq_idx);
        
        figure('Position', [100, 100, 1200, 400]);
        
        % Extract results for this frequency
        cbf_result = squeeze(cbf_results(freq_idx, :, :));
        omp_result = squeeze(omp_damas_results(freq_idx, :, :));
        sc_result = squeeze(sc_damas_results(freq_idx, :, :));
        
        % Normalize for better visualization
        cbf_result = cbf_result / max(cbf_result(:));
        omp_result = omp_result / max(omp_result(:));
        sc_result = sc_result / max(sc_result(:));
        
        % Plot CBF results
        subplot(1, 3, 1);
        surf(cbf_result);
        title(['CBF at ', num2str(freq), ' Hz']);
        zlim([0 1]);
        xlabel('X grid index');
        ylabel('Z grid index');
        
        % Plot OMP-DAMAS results
        subplot(1, 3, 2);
        surf(omp_result);
        title(['OMP-DAMAS at ', num2str(freq), ' Hz']);
        zlim([0 1]);
        xlabel('X grid index');
        ylabel('Z grid index');
        
        % Plot SC-DAMAS results
        subplot(1, 3, 3);
        surf(sc_result);
        title(['SC-DAMAS at ', num2str(freq), ' Hz']);
        zlim([0 1]);
        xlabel('X grid index');
        ylabel('Z grid index');
        
        sgtitle(['3D Visualization at ', num2str(freq), ' Hz, SNR = ', num2str(SNR_dB), ' dB']);
        % Save 3D figure immediately after creation
        if save_results
            figure_name = sprintf('beamforming_3D_%dHz_%s', freq, timestamp);
            saveas(gcf, fullfile(figures_dir, [figure_name '.png']), 'png');
            fprintf('Figures saved for %d Hz with timestamp %s\n', freq, timestamp);
        end
    end
end

%% Find peak positions for comparison with Table 3
% Initialize arrays to store estimated source positions
cbf_peak_pos = zeros(length(freq_to_analyze), 2);
omp_peak_pos = zeros(length(freq_to_analyze), 2);
sc_peak_pos = zeros(length(freq_to_analyze), 2);

% Find peak positions for each method at each frequency
for freq_idx = 1:length(freq_to_analyze)
    % CBF
    cbf_result = squeeze(cbf_results(freq_idx, :, :));
    [~, max_idx] = max(cbf_result(:));
    [row, col] = ind2sub(size(cbf_result), max_idx);
    cbf_peak_pos(freq_idx, :) = [x_grid(row, col), z_grid(row, col)];
    
    % OMP-DAMAS
    omp_result = squeeze(omp_damas_results(freq_idx, :, :));
    [~, max_idx] = max(omp_result(:));
    [row, col] = ind2sub(size(omp_result), max_idx);
    omp_peak_pos(freq_idx, :) = [x_grid(row, col), z_grid(row, col)];
    
    % SC-DAMAS
    sc_result = squeeze(sc_damas_results(freq_idx, :, :));
    [~, max_idx] = max(sc_result(:));
    [row, col] = ind2sub(size(sc_result), max_idx);
    sc_peak_pos(freq_idx, :) = [x_grid(row, col), z_grid(row, col)];
end

% Display results in a table format (similar to Table 3 in the paper)
disp('Sound Source Localization Results:');
disp('------------------------------------------------------');
disp('Frequency (Hz) | CBF (x,z) | OMP-DAMAS (x,z) | SC-DAMAS (x,z)');
disp('------------------------------------------------------');
for freq_idx = 1:length(freq_to_analyze)
    fprintf('%13d | (%.2f, %.2f) | (%.2f, %.2f) | (%.2f, %.2f)\n', ...
        freq_to_analyze(freq_idx), ...
        cbf_peak_pos(freq_idx, 1), cbf_peak_pos(freq_idx, 2), ...
        omp_peak_pos(freq_idx, 1), omp_peak_pos(freq_idx, 2), ...
        sc_peak_pos(freq_idx, 1), sc_peak_pos(freq_idx, 2));
end
disp('------------------------------------------------------');
disp(['Actual source position: (', num2str(sourcePos(1)), ', ', num2str(sourcePos(3)), ')']);

% Calculate localization errors
cbf_errors = sqrt(sum((cbf_peak_pos - [sourcePos(1), sourcePos(3)]).^2, 2));
omp_errors = sqrt(sum((omp_peak_pos - [sourcePos(1), sourcePos(3)]).^2, 2));
sc_errors = sqrt(sum((sc_peak_pos - [sourcePos(1), sourcePos(3)]).^2, 2));

disp('Localization Errors (Euclidean distance):');
disp('------------------------------------------------------');
disp('Frequency (Hz) | CBF | OMP-DAMAS | SC-DAMAS');
disp('------------------------------------------------------');
for freq_idx = 1:length(freq_to_analyze)
    fprintf('%13d | %.3f | %.3f | %.3f\n', ...
        freq_to_analyze(freq_idx), ...
        cbf_errors(freq_idx), ...
        omp_errors(freq_idx), ...
        sc_errors(freq_idx));
end
disp('------------------------------------------------------');

%% Save results with timestamp to prevent overwriting
if save_results
    % Create timestamp for unique filenames
    
    
    % Create beamforming results directory if it doesn't exist
    results_dir = 'beamforming_results';
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    
    % Save .mat file with timestamp
    results_filename = fullfile(results_dir, sprintf('sc_damas_results_%s.mat', timestamp));
    
    % Add metadata to saved results
    beamforming_metadata = struct();
    beamforming_metadata.sourcePos = sourcePos;
    beamforming_metadata.micPositions = micPositions;
    beamforming_metadata.timestamp = timestamp;
    beamforming_metadata.SNR_dB = SNR_dB;
    if exist('pos_idx', 'var')
        beamforming_metadata.position_index = pos_idx;
    end
    
    save(results_filename, 'cbf_results', 'omp_damas_results', 'sc_damas_results', ...
        'cbf_peak_pos', 'omp_peak_pos', 'sc_peak_pos', ...
        'cbf_errors', 'omp_errors', 'sc_errors', ...
        'freq_to_analyze', 'grid_points', 'x_grid', 'z_grid', ...
        'beamforming_metadata');
    
    fprintf('Results saved to %s\n', results_filename);
end

disp('SC-DAMAS implementation complete.');