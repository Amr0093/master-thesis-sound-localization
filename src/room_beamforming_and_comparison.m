% SC-DAMAS Implementation for Acoustic Imaging in Reverberant Field
% This script implements the SC-DAMAS algorithm described in the paper
% "An improved deconvolution beamforming algorithm for acoustic imaging
% of low signal-to-noise ratio sound sources in reverberant field"
% VERSION 2 - WITH ROBUST NUMERICAL STABILITY

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

%% Get simulated microphone data from simulator workspace
% Variables should already exist from FastHybridDopplerReverbSimulator_enhanced_2
if ~exist('final_signals', 'var')
    error('Beamforming requires simulator data. Run FastHybridDopplerReverbSimulator_enhanced_2 first.');
end

%% DEBUG - Check array dimensions
fprintf('=== BEAMFORMING DEBUG INFO ===\n');
fprintf('final_signals size: [%d x %d]\n', size(final_signals, 1), size(final_signals, 2));
fprintf('micPositions size: [%d x %d]\n', size(micPositions, 1), size(micPositions, 2));
fprintf('numMics: %d\n', numMics);
fprintf('sourcePos: [%.2f, %.2f, %.2f]\n', sourcePos(1), sourcePos(2), sourcePos(3));
fprintf('grid_size: %d\n', grid_size);
fprintf('==============================\n');

% Safety checks
if numMics < 2
    error('Beamforming requires at least 2 microphones, found %d', numMics);
end

if size(micPositions, 2) ~= 3
    error('micPositions must be [N x 3] array, found [%d x %d]', size(micPositions, 1), size(micPositions, 2));
end

if length(sourcePos) ~= 3
    error('sourcePos must be [1 x 3] array, found [1 x %d]', length(sourcePos));
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
delta_criteria = 1e-6;        % Convergence criteria for iterations
max_iterations = 50;          % Maximum iterations for DAMAS
plot_results = true;          % Set to true to plot results
save_results = true;          % Set to true to save results

%% Robust Parameters & Checks - NEW IN V2
% Add these parameters for numerical stability
regularization_param = 1e-8;     % For matrix conditioning
max_condition_number = 1e12;     % Maximum allowed condition number
min_eigenvalue_ratio = 1e-10;    % Minimum eigenvalue threshold
timeout_seconds = 30;            % Maximum processing time per frequency
enable_progress_display = true;   % Show progress updates

% Create timestamp for unique filenames - MOVED UP IN V2
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% Start timer for overall processing - NEW IN V2
total_start_time = tic;

%% Create scanning grid
% Grid centered at actual source position for testing/validation
grid_center = [sourcePos(1), sourcePos(3)];  % x, z coordinates

% Define scanning grid (xz-plane at source's y position)
[x_grid, z_grid] = meshgrid(linspace(grid_center(1)-scan_area_size/2, grid_center(1)+scan_area_size/2, grid_size), ...
                            linspace(grid_center(2)-scan_area_size/2, grid_center(2)+scan_area_size/2, grid_size));
grid_points = [x_grid(:), repmat(sourcePos(2), grid_size^2, 1), z_grid(:)];

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
    
    %% SC-DAMAS algorithm (SAMP-CoSaMP based DAMAS) - ROBUST VERSION V2
    disp('Calculating SC-DAMAS...');
    freq_start_time = tic;

    % Initialize SC-DAMAS variables
    r = b_cbf;                  % Initial residual
    stage = 1;                  % Initial stage for sparsity adaptation
    L = sparsity_step;          % Initial sparsity estimate (step size)
    F = zeros(0,1);            % Support set (empty column vector)
    k = 1;                      % Iteration counter

    % Check for timeout - NEW IN V2
    if toc(total_start_time) > timeout_seconds * freq_idx
        warning('Timeout approaching, using simplified solution for frequency %d Hz', freq);
        x = b_cbf / max(b_cbf(:)) * 0.1;  % Simplified fallback
        sc_damas_results(freq_idx, :, :) = reshape(x, grid_size, grid_size);
        continue;  % Skip to next frequency
    end

    % Step 1: Initial sparsity estimation via SAMP (with checks) - ENHANCED IN V2
    try
        fprintf('DEBUG: Starting SC-DAMAS, F initial size: %s\n', mat2str(size(F)));
        % Calculate initial correlation with bounds checking - NEW IN V2
        u = abs(PSF' * r);
        if any(~isfinite(u))
            warning('Non-finite correlation values detected, using fallback');
            u = abs(randn(size(u)) * 0.01);  % Small random fallback
        end
        
        [~, indices] = sort(u, 'descend');
        K = max(1, min(10, floor(grid_size^2/10)));  % Safer initial sparsity - CHANGED IN V2
        
        % Find initial support set with safety checks - ENHANCED IN V2
        iteration_count = 0;
        while iteration_count < 20  % Prevent infinite loops - NEW IN V2
            iteration_count = iteration_count + 1;
            
            fprintf('DEBUG: Before F assignment - indices size: %s, K=%d\n', mat2str(size(indices)), K);
            F = indices(1:min(K, length(indices)));
            F = F(:);  % Ensure column vector
            fprintf('DEBUG: After F assignment - F size: %s\n', mat2str(size(F)));
            F = F(:);  % Ensure F is always a column vector
            
            % Use robust matrix solving - NEW IN V2
            PSF_F = PSF(:, F);
            x_F = zeros(grid_size^2, 1);
            
            % Robust solve instead of direct pinv - NEW IN V2
            x_F_solution = robust_matrix_solve(PSF_F, b_cbf, regularization_param, max_condition_number, 'SC-DAMAS Init');
            x_F(F) = x_F_solution;
            
            % Check if we need to increase sparsity
            residual_norm = norm(b_cbf - PSF * x_F);
            
            % More conservative delta_K - CHANGED IN V2
            delta_K = 0.3;  % Increased from 0.1 for stability
            
            if residual_norm <= ((1-delta_K)/sqrt(1+delta_K)) * norm(b_cbf) || K >= grid_size^2/2
                break;
            else
                K = min(K + L, grid_size^2);
            end
            
            % Show progress if enabled - NEW IN V2
            if enable_progress_display && mod(iteration_count, 5) == 0
                fprintf('  Init iteration %d, K=%d, residual=%.2e\n', iteration_count, K, residual_norm);
            end
        end
        
        % Initial residual
        r = b_cbf - PSF * x_F;
        x = x_F;
        
        % Step 2: Main SC-DAMAS iterations (with timeout and stability checks) - ENHANCED IN V2
        while k <= max_iterations
            % Check for timeout - NEW IN V2
            if toc(freq_start_time) > timeout_seconds
                warning('SC-DAMAS timeout for frequency %d Hz, using current solution', freq);
                break;
            end
            
            % Calculate correlation with bounds checking - NEW IN V2
            u = abs(PSF' * r);
            if any(~isfinite(u))
                warning('Non-finite correlation in iteration %d, breaking', k);
                break;
            end
            
            % More conservative threshold - CHANGED IN V2
            threshold = max(0.1 * max(u), 1e-10);  % Prevent threshold from being too small
            candidate_set = find(u >= threshold);
            
            if isempty(candidate_set)
                if enable_progress_display
                    fprintf('  No candidates found, breaking at iteration %d\n', k);
                end
                break;
            end
            
            % Choose the top K components
            [~, indices] = sort(u(candidate_set), 'descend');
            num_to_select = min(K, length(indices), length(candidate_set));
            top_indices = candidate_set(indices(1:num_to_select));
            
            % Debug the dimensions before union
            fprintf('DEBUG: F size: %s, top_indices size: %s\n', mat2str(size(F)), mat2str(size(top_indices)));
            fprintf('DEBUG: candidate_set size: %s, indices size: %s\n', mat2str(size(candidate_set)), mat2str(size(indices)));
            
            % Update support set
            % Update support set - with extra safety
            F = F(:);  % Force F to column vector  
            top_indices = top_indices(:);  % Force top_indices to column vector
            % Ultra-safe union operation
            if isempty(F)
                F_new = top_indices(:);  % Force column
            elseif isempty(top_indices)
                F_new = F(:);           % Force column  
            else
                % Debug before union call
                fprintf('DEBUG: About to call union - F size: %s, top_indices size: %s\n', mat2str(size(F)), mat2str(size(top_indices)));
                F_new = union(F(:), top_indices(:)); % Force both to column
                fprintf('DEBUG: Union completed successfully\n');
            end
            fprintf('DEBUG: F_new size: %s\n', mat2str(size(F_new)));
            
            % Limit support set size to prevent overfitting - NEW IN V2
            if length(F_new) > grid_size^2/3
                [~, best_indices] = sort(u(F_new), 'descend');
                F_new = F_new(best_indices(1:floor(grid_size^2/3)));
            end
            
            % Calculate new estimate using robust solving - NEW IN V2
            PSF_F_new = PSF(:, F_new);
            x_new = zeros(grid_size^2, 1);
            
            % Robust solve instead of direct pinv - NEW IN V2
            x_new_solution = robust_matrix_solve(PSF_F_new, b_cbf, regularization_param, max_condition_number, 'SC-DAMAS Main');
            x_new(F_new) = x_new_solution;
            
            % Calculate new residual
            r_new = b_cbf - PSF * x_new;
            
            % Check convergence with relative criteria - ENHANCED IN V2
            relative_change = norm(x_new - x) / (norm(x) + 1e-10);
            if relative_change < delta_criteria
                x = x_new;
                if enable_progress_display
                    fprintf('  Converged at iteration %d (relative change: %.2e)\n', k, relative_change);
                end
                break;
            end
            
            % Adaptive sparsity adjustment - ENHANCED IN V2
            if norm(r_new) >= 0.95 * norm(r)  % Less aggressive than >= norm(r)
                stage = stage + 1;
                L = min(stage * sparsity_step, 5);  % Limit step growth
                K = min(K + L, grid_size^2/2);     % Don't exceed half the grid
            else
                F = F_new;
                r = r_new;
                x = x_new;
                k = k + 1;
            end
            
            % Show progress - NEW IN V2
            if enable_progress_display && mod(k, 10) == 0
                fprintf('  Iteration %d: residual=%.2e, support_size=%d\n', k, norm(r), length(F));
            end
        end
        
    catch ME
        warning('SC-DAMAS failed for frequency %d Hz: %s. Using CBF result.', freq, ME.message);
        x = b_cbf / max(abs(b_cbf(:))) * 0.1;  % Normalized CBF fallback
    end

    % Ensure non-negativity and finite values - ENHANCED IN V2
    x(x < 0) = 0;
    x(~isfinite(x)) = 0;

    % Reshape and store SC-DAMAS results
    sc_damas_results(freq_idx, :, :) = reshape(x, grid_size, grid_size);

    fprintf('SC-DAMAS completed for %d Hz (%.1f seconds)\n', freq, toc(freq_start_time));
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
        % Create comprehensive figure directory structure
        if save_results
            % Create beamforming results directory if it doesn't exist
            results_dir = 'beamforming_results';  % ← THIS LINE MUST COME FIRST
            if ~exist(results_dir, 'dir')
                mkdir(results_dir);
                
            end
            
            % Now you can use results_dir safely
            figures_dir = fullfile(results_dir, 'figures');
            if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
            
            % Save figure list for GUI access
            figure_list = {};
            for freq_idx = 1:length(freq_to_analyze)
                freq = freq_to_analyze(freq_idx);
                figure_list{end+1} = sprintf('beamforming_2D_%dHz_%s.png', freq, timestamp);
                figure_list{end+1} = sprintf('beamforming_3D_%dHz_%s.png', freq, timestamp);
            end
            
            % Save figure list in the results directory
            save(fullfile(results_dir, 'figure_list.mat'), 'figure_list', 'figures_dir');
            fprintf('✅ Figure directory created: %s\n', figures_dir);
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
    beamforming_metadata.version = 'v2_robust';  % NEW IN V2
    beamforming_metadata.robust_params = struct('regularization_param', regularization_param, ...
                                               'max_condition_number', max_condition_number, ...
                                               'timeout_seconds', timeout_seconds);  % NEW IN V2

    
    % Create structured results like IR analysis does
    beamforming_results = struct();
    beamforming_results.metadata = beamforming_metadata;
    beamforming_results.algorithms = struct();
    beamforming_results.algorithms.cbf = struct('results', cbf_results, 'peak_positions', cbf_peak_pos, 'errors', cbf_errors);
    beamforming_results.algorithms.omp_damas = struct('results', omp_damas_results, 'peak_positions', omp_peak_pos, 'errors', omp_errors);
    beamforming_results.algorithms.sc_damas = struct('results', sc_damas_results, 'peak_positions', sc_peak_pos, 'errors', sc_errors);
    beamforming_results.grid = struct('x_grid', x_grid, 'z_grid', z_grid, 'points', grid_points);
    beamforming_results.frequencies = freq_to_analyze;
    
    % Create unified structured results (like IR analysis)
    beamforming_results = struct();
    beamforming_results.metadata = beamforming_metadata;
    beamforming_results.metadata.analysis_timestamp = datetime('now');
    beamforming_results.algorithms = struct();
    beamforming_results.algorithms.cbf = struct('results', cbf_results, 'peak_positions', cbf_peak_pos, 'errors', cbf_errors);
    beamforming_results.algorithms.omp_damas = struct('results', omp_damas_results, 'peak_positions', omp_peak_pos, 'errors', omp_errors);
    beamforming_results.algorithms.sc_damas = struct('results', sc_damas_results, 'peak_positions', sc_peak_pos, 'errors', sc_errors);
    beamforming_results.grid = struct('x_grid', x_grid, 'z_grid', z_grid, 'points', grid_points);
    beamforming_results.frequencies = freq_to_analyze;
    
    % Save both formats for compatibility
    save(results_filename, 'cbf_results', 'omp_damas_results', 'sc_damas_results', ...
        'cbf_peak_pos', 'omp_peak_pos', 'sc_peak_pos', ...
        'cbf_errors', 'omp_errors', 'sc_errors', ...
        'freq_to_analyze', 'grid_points', 'x_grid', 'z_grid', ...
        'beamforming_metadata', 'beamforming_results');
    
    % Store in base workspace for immediate GUI access
    assignin('base', 'beamforming_results', beamforming_results);
    assignin('base', 'latest_beamforming_folder', results_dir);
    
    fprintf('✅ Unified results structure created and saved\n');
    
   
    
    fprintf('Results saved to %s\n', results_filename);
end

disp('SC-DAMAS implementation complete.');

%% Robust Matrix Solving Function - NEW IN V2
function x_result = robust_matrix_solve(A, b, reg_param, max_cond, method_name)
    % Robust matrix solving with multiple fallback strategies
    
    try
        % Check matrix dimensions
        [m, n] = size(A);
        if m == 0 || n == 0
            warning('%s: Empty matrix encountered', method_name);
            x_result = zeros(n, 1);
            return;
        end
        
        % Strategy 1: Check condition number
        cond_num = cond(A);
        if cond_num > max_cond || ~isfinite(cond_num)
            fprintf('%s: Poor conditioning (cond=%.2e), using regularization\n', method_name, cond_num);
            % Add regularization
            A_reg = A' * A + reg_param * eye(size(A, 2));
            x_result = A_reg \ (A' * b);
        else
            % Strategy 2: Standard pseudo-inverse
            x_result = pinv(A) * b;
        end
        
        % Check for NaN or Inf results
        if any(~isfinite(x_result))
            warning('%s: Non-finite result, using fallback', method_name);
            % Fallback: Simple least squares with heavy regularization
            A_safe = A' * A + (reg_param * 100) * eye(size(A, 2));
            x_result = A_safe \ (A' * b);
        end
        
    catch ME
        warning('%s: Matrix solve failed (%s), using zero solution', method_name, ME.message);
        x_result = zeros(size(A, 2), 1);
    end
    % Set the results folder for GUI access
    assignin('base', 'latest_beamforming_folder', results_dir);
    fprintf('✅ Results folder set for GUI: %s\n', results_dir);
end