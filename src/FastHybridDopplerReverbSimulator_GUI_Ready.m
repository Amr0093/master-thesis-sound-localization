function results = FastHybridDopplerReverbSimulator_GUI_Ready(params)
%% Enhanced Fast Hybrid Doppler-Reverb Simulator - GUI Ready Version
% This function implements a hybrid approach that combines:
% 1. NOTAP-inspired source-time dominant approach for Doppler effects
% 2. Time/Frequency processing for efficient room reflections with accurate Doppler 
% 3. Complete source-receiver mapping data structures for localization
% 
% Enhanced version with:
% - Multiple room shape options (rectangular room or cylindrical tunnel)
% - Multiple microphone array geometries (rectangular grid, linear, circular, spiral, spherical)
% - Advanced noise modeling options
% - F-parameter calculation for array quality assessment
% - Support for even/odd microphone counts
% - Optional simplified vehicle body
% - GUI-ready parameter structure and organized outputs
%
% INPUT:
%   params - Structure containing all simulation parameters from GUI
%
% OUTPUT:
%   results - Structure containing all simulation results and file paths
fprintf('Debug: Simulator received these fields:\n');
disp(fieldnames(params));

if isfield(params, 'array_center_x')
    fprintf('Debug: params.array_center_x = %g\n', params.array_center_x);
else
    fprintf('Debug: params.array_center_x NOT FOUND!\n');
end


% Check start_pos
if isfield(params, 'start_pos')
    fprintf('Debug: start_pos = [%g, %g, %g]\n', params.start_pos(1), params.start_pos(2), params.start_pos(3));
else
    fprintf('Debug: start_pos NOT FOUND!\n');
end

% Check end_pos  
if isfield(params, 'end_pos')
    fprintf('Debug: end_pos = [%g, %g, %g]\n', params.end_pos(1), params.end_pos(2), params.end_pos(3));
else
    fprintf('Debug: end_pos NOT FOUND!\n');
end

% Check alternative field names (in case they're different)
if isfield(params, 'start_x')
    fprintf('Debug: Found start_x instead = %g\n', params.start_x);
end

if isfield(params, 'startPosition')
    fprintf('Debug: Found startPosition instead = %s\n', params.startPosition);
end

%% Initialize Results Structure
results = struct();
results.success = false;
results.error_message = '';
results.execution_time = 0;
results.files_created = {};
results.data = struct();

try
    start_time = tic;
    
    %% Parameter Processing and Validation
    fprintf('\n=== GUI-Ready Fast Hybrid Doppler-Reverb Simulator ===\n');
    fprintf('Processing input parameters...\n');
    
    % Process and validate input parameters
    [processed_params, validation_result] = processInputParameters(params);
    
    if ~validation_result.valid
        results.error_message = sprintf('Parameter validation failed: %s', validation_result.message);
        return;
    end
    
    % Store processed parameters in results
    results.parameters = processed_params;
    
    %% Check for external call and prepare
    if exist('EXTERNAL_CALL_FLAG', 'var') && EXTERNAL_CALL_FLAG
        fprintf('Simulator: Received external call, cleaning workspace...\n');
        clear;
        close all;
        clc;
        % Send "ready" signal by throwing a specific error
        error('Simulator_check_completed'); % This is our "handshake" signal
    end

    %% Extract processed parameters for simulation
    % Room Configuration
    room_type = processed_params.room_type;
    roomDim = processed_params.roomDim;
    c = processed_params.c;
    fs = processed_params.fs;
    duration = processed_params.duration;
    
    % Source Configuration
    start_pos = processed_params.start_pos;
    end_pos = processed_params.end_pos;
    

    signal_choice = processed_params.signal_choice;
    source_power = processed_params.source_power;
    
    % Array Configuration
    array_type = processed_params.array_type;
    array_orientation = processed_params.array_orientation;
    numMics = processed_params.numMics;
    arrayCenter = processed_params.arrayCenter;
    
    % Simulation Options
    noise_type = processed_params.noise_type;
    enable_noise = processed_params.enable_noise;
    use_vehicle_body = processed_params.use_vehicle_body;
    enable_visualization = processed_params.enable_visualization;
    use_freq_domain_reflections = processed_params.use_freq_domain_reflections;
    use_batch_processing = processed_params.use_batch_processing;
    max_reflection_order = processed_params.max_reflection_order;
    % batch_size = processed_params.batch_size;
    oversample_ratio = processed_params.oversample_ratio;
    save_mappings_for_localization = processed_params.save_mappings_for_localization;
    reg = processed_params.reg;
    
    % Additional room parameters for cylindrical
    if strcmp(room_type, 'cylindrical')
        tunnel_length = processed_params.tunnel_length;
        tunnel_radius = processed_params.tunnel_radius;
    end

    %% Display Configuration Summary
    fprintf('\nSimulation Configuration:\n');
    fprintf('  Room type: %s\n', room_type);
    fprintf('  Array type: %s (%d microphones)\n', array_type, numMics);
    fprintf('  Array orientation: %s\n', array_orientation);
    fprintf('  Enable noise: %s', string(enable_noise));
    if enable_noise
        fprintf(' (%s)\n', noise_type);
    else
        fprintf('\n');
    end
    fprintf('  Vehicle body: %s\n', string(use_vehicle_body));
    fprintf('  Visualization: %s\n', string(enable_visualization));
    fprintf('  Freq domain reflections: %s\n', string(use_freq_domain_reflections));
    fprintf('  Batch processing: %s\n', string(use_batch_processing));
    fprintf('  Max reflection order: %d\n', max_reflection_order);
    fprintf('  Save localization mappings: %s\n', string(save_mappings_for_localization));

    %% Motion Analysis
    is_moving = ~isequal(start_pos, end_pos);
    
    if is_moving
        v_vec = (end_pos - start_pos) / duration;
        source_trajectory = @(t) start_pos + v_vec * t;
        fprintf('\nSource Motion:\n');
        fprintf('  From: [%.2f, %.2f, %.2f] to [%.2f, %.2f, %.2f]\n', ...
            start_pos(1), start_pos(2), start_pos(3), end_pos(1), end_pos(2), end_pos(3));
        fprintf('  Velocity: [%.2f, %.2f, %.2f] m/s (%.2f m/s magnitude)\n', ...
            v_vec(1), v_vec(2), v_vec(3), norm(v_vec));
    else
        source_trajectory = @(t) start_pos;
        fprintf('\nSource Motion: Stationary at [%.2f, %.2f, %.2f]\n', ...
            start_pos(1), start_pos(2), start_pos(3));
    end

    % Vehicle body parameters (if enabled)
    if use_vehicle_body
        vehicle_length = 4.5;
        vehicle_width = 1.8;
        vehicle_height = 1.5;
        vehicle_trajectory = @(t) source_trajectory(t);
        fprintf('  Vehicle body enabled: %.1f×%.1f×%.1f m\n', vehicle_length, vehicle_width, vehicle_height);
    end

    %% Define Sensor Array
    fprintf('\nConfiguring %d-microphone %s array...\n', numMics, array_type);
    
    array_size = 0.1;
    array_elements = 2;
    element_spacing = array_size/(array_elements-1);
    
    % Initialize microphone positions
    micPositions = zeros(numMics, 3);
    
    % Set microphone positions based on array type
    switch array_type
        case 'rectangular_grid'
            [x_grid, z_grid] = meshgrid(-(array_size/2):element_spacing:(array_size/2), ...
                                       -(array_size/2):element_spacing:(array_size/2));
            
            switch array_orientation
                case 'x-axis'
                    array_x = repmat(arrayCenter(1), size(x_grid));
                    array_y = arrayCenter(2) + x_grid;
                    array_z = arrayCenter(3) + z_grid;
                case 'y-axis'
                    array_x = arrayCenter(1) + x_grid;
                    array_y = repmat(arrayCenter(2), size(x_grid));
                    array_z = arrayCenter(3) + z_grid;
                case 'z-axis'
                    array_x = arrayCenter(1) + x_grid;
                    array_y = arrayCenter(2) + z_grid;
                    array_z = repmat(arrayCenter(3), size(x_grid));
            end
            
            counter = 1;
            for i = 1:size(array_x, 1)
                for j = 1:size(array_x, 2)
                    micPositions(counter, :) = [array_x(i,j), array_y(i,j), array_z(i,j)];
                    counter = counter + 1;
                end
            end
        case 'single'  % NEW: Explicit single microphone array type
            micPositions(1, :) = arrayCenter;
            fprintf('Debug: Explicit single microphone array at [%.2f, %.2f, %.2f]\n', ...
                arrayCenter(1), arrayCenter(2), arrayCenter(3));
            
        case 'linear'
            if numMics == 1
                % Single microphone at array center
                micPositions(1, :) = arrayCenter;
                fprintf('Debug: Single microphone positioned at [%.2f, %.2f, %.2f]\n', ...
                    arrayCenter(1), arrayCenter(2), arrayCenter(3));
            else
                % Multiple microphones in linear array
                positions = linspace(-array_size/2, array_size/2, numMics);
                for i = 1:numMics
                    switch array_orientation
                        case 'x-axis'
                            micPositions(i, :) = [arrayCenter(1), arrayCenter(2) + positions(i), arrayCenter(3)];
                        case 'y-axis'
                            micPositions(i, :) = [arrayCenter(1) + positions(i), arrayCenter(2), arrayCenter(3)];
                        case 'z-axis'
                            micPositions(i, :) = [arrayCenter(1), arrayCenter(2), arrayCenter(3) + positions(i)];
                    end
                end
                
            end
            
        case 'circular'
            if numMics == 1
                % Single microphone at array center
                micPositions(1, :) = arrayCenter;
            else
                % Multiple microphones in circular array
                radius = array_size/2;
                angles = linspace(0, 2*pi, numMics+1);
                angles = angles(1:end-1);
                
                for i = 1:numMics
                    switch array_orientation
                        case 'x-axis'
                            micPositions(i, :) = [arrayCenter(1), ...
                                                 arrayCenter(2) + radius*cos(angles(i)), ...
                                                 arrayCenter(3) + radius*sin(angles(i))];
                        case 'y-axis'
                            micPositions(i, :) = [arrayCenter(1) + radius*cos(angles(i)), ...
                                                 arrayCenter(2), ...
                                                 arrayCenter(3) + radius*sin(angles(i))];
                        case 'z-axis'
                            micPositions(i, :) = [arrayCenter(1) + radius*cos(angles(i)), ...
                                                 arrayCenter(2) + radius*sin(angles(i)), ...
                                                 arrayCenter(3)];
                    end
                end
            end
            
        case 'spiral'
            if numMics == 1
                % Single microphone at array center
                micPositions(1, :) = arrayCenter;
            else
                % Multiple microphones in spiral array
                radius = array_size/2;
                num_arms = 4;
                points_per_arm = ceil(numMics / num_arms);
                
                mic_idx = 1;
                for arm = 1:num_arms
                    phi0 = 2*pi*(arm-1)/num_arms;
                    for pt = 1:points_per_arm
                        if mic_idx <= numMics
                            r = radius * sqrt(pt/points_per_arm);
                            phi = phi0 + 0.1*r;
                            
                            switch array_orientation
                                case 'x-axis'
                                    micPositions(mic_idx, :) = [arrayCenter(1), ...
                                                              arrayCenter(2) + r*cos(phi), ...
                                                              arrayCenter(3) + r*sin(phi)];
                                case 'y-axis'
                                    micPositions(mic_idx, :) = [arrayCenter(1) + r*cos(phi), ...
                                                              arrayCenter(2), ...
                                                              arrayCenter(3) + r*sin(phi)];
                                case 'z-axis'
                                    micPositions(mic_idx, :) = [arrayCenter(1) + r*cos(phi), ...
                                                              arrayCenter(2) + r*sin(phi), ...
                                                              arrayCenter(3)];
                            end
                            
                            mic_idx = mic_idx + 1;
                        end
                    end
                end
            end
            
        case 'spherical'
            if numMics == 1
                % Single microphone at array center
                micPositions(1, :) = arrayCenter;
            else
                % Multiple microphones in spherical array
                radius = array_size/2;
                golden_ratio = (1 + sqrt(5)) / 2;
                indices = 0:numMics-1;
                
                phi = 2*pi * indices / golden_ratio;
                costheta = 1 - (2*indices + 1) / numMics;
                sintheta = sqrt(1 - costheta.^2);
                
                for i = 1:numMics
                    micPositions(i, :) = [arrayCenter(1) + radius * sintheta(i) * cos(phi(i)), ...
                                         arrayCenter(2) + radius * sintheta(i) * sin(phi(i)), ...
                                         arrayCenter(3) + radius * costheta(i)];
                end
            end
    end

    %% Calculate F-parameter (array quality metric)
    fprintf('Calculating array F-parameter...\n');
    
    co_array = zeros(numMics^2, 3);
    counter = 1;
    for i = 1:numMics
        for j = 1:numMics
            co_array(counter, :) = micPositions(i, :) - micPositions(j, :);
            counter = counter + 1;
        end
    end

    unique_vectors = 0;
    tolerance = 1e-10;
    is_unique = false(size(co_array, 1), 1);

    for i = 1:size(co_array, 1)
        if ~any(is_unique)
            is_unique(i) = true;
            unique_vectors = unique_vectors + 1;
        else
            unique_indices = find(is_unique);
            is_new = true;
            
            for j = 1:length(unique_indices)
                if norm(co_array(i, :) - co_array(unique_indices(j), :)) < tolerance
                    is_new = false;
                    break;
                end
            end
            
            if is_new
                is_unique(i) = true;
                unique_vectors = unique_vectors + 1;
            end
        end
    end

    max_unique_vectors = numMics^2 - (numMics - 1);
    F_parameter = unique_vectors / max_unique_vectors;
    fprintf('Array F-parameter (unicity): %.4f\n', F_parameter);

    %% Define Acoustic Properties
    freq = [125, 250, 500, 1000, 2000, 4000];
    
    if strcmp(room_type, 'rectangular')
        absorption_coef = zeros(6, length(freq));
        absorption_coef(1:2, :) = repmat([0.18, 0.06, 0.04, 0.03, 0.02, 0.02], 2, 1);
        absorption_coef(3:4, :) = repmat([0.1, 0.05, 0.06, 0.07, 0.09, 0.08], 2, 1);
        absorption_coef(5, :) = [0.01, 0.01, 0.01, 0.01, 0.02, 0.02];
        absorption_coef(6, :) = [0.35, 0.25, 0.18, 0.12, 0.07, 0.04];
    else
        absorption_coef = zeros(6, length(freq));
        absorption_coef(1, :) = [0.03, 0.03, 0.03, 0.04, 0.05, 0.07];
        absorption_coef(2:3, :) = repmat([0.01, 0.01, 0.01, 0.01, 0.02, 0.02], 2, 1);
        absorption_coef(5, :) = [0.02, 0.02, 0.03, 0.03, 0.05, 0.05];
    end
    
    % Calculate reflection coefficients
    reflection_coef = sqrt(1 - absorption_coef);
    freq_analysis = 1000;
    analysis_frequencies = [125, 250, 500, 1000, 2000, 4000];
    num_freq_bands = length(analysis_frequencies);
    [~, freq_idx] = min(abs(freq - freq_analysis));
    
    % ENHANCED: Create frequency-dependent beta coefficients structure
    if strcmp(room_type, 'rectangular')
        beta_coefficients = struct();
        beta_coefficients.x1 = reflection_coef(1, :); % Front wall - ALL frequencies
        beta_coefficients.x2 = reflection_coef(3, :); % Right wall - ALL frequencies  
        beta_coefficients.y1 = reflection_coef(6, :); % Back wall - ALL frequencies
        beta_coefficients.y2 = reflection_coef(2, :); % Left wall - ALL frequencies
        beta_coefficients.z1 = reflection_coef(5, :); % Ground - ALL frequencies
        beta_coefficients.z2 = reflection_coef(4, :); % Ceiling - ALL frequencies
        

        % Also keep single coefficients for backwards compatibility
        beta_x1 = reflection_coef(1, freq_idx);
        beta_x2 = reflection_coef(3, freq_idx);
        beta_y1 = reflection_coef(6, freq_idx);
        beta_y2 = reflection_coef(2, freq_idx);
        beta_z1 = reflection_coef(5, freq_idx);
        beta_z2 = reflection_coef(4, freq_idx);
    else
        % ENHANCED: Create frequency-dependent coefficients for cylindrical room
        beta_tunnel_wall = reflection_coef(1, :); % ALL frequencies
        beta_entrance = reflection_coef(2, :);    % ALL frequencies
        beta_exit = reflection_coef(3, :);        % ALL frequencies
        beta_floor = reflection_coef(5, :);       % ALL frequencies
        
        fprintf('Debug: Cylindrical tunnel reflection coefficients loaded for %d frequency bands\n', num_freq_bands);
    end
    %% Load or Generate Source Signal
    fprintf('Generating source signal (%s)...\n', signal_choice);
    
    switch signal_choice
        case 'siren'
            try
                [source_signal, src_fs] = audioread('siren_sound/original.wav');
                if size(source_signal, 2) > 1
                    source_signal = mean(source_signal, 2);
                end
                if src_fs ~= fs
                    source_signal = resample(source_signal, fs, src_fs);
                end
                fprintf('Loaded siren sound, length: %.2f seconds\n', length(source_signal)/fs);
            catch
                error('Siren sound file not found!');
            end
            
        case 'chirp'
            Tsig = 0.5;
            t_sig = 0:1/fs:Tsig-1/fs;
            source_signal = chirp(t_sig, 500, Tsig, 2000)';
            
        case 'measurement_sweep'
            signalDuration = duration;
            freqStart = 20;
            freqEnd = 20000;
            source_signal = generateExpSweep(fs, signalDuration, freqStart, freqEnd);
            
        otherwise
            error('Unknown signal choice: %s', signal_choice);
    end

    reference_distance = 1.0;
    scaling_factor = sqrt(source_power/(4*pi*reference_distance^2));
    source_signal = source_signal * scaling_factor;

    %% Timing Setup
    source_duration = length(source_signal)/fs;
    t_s = (0:length(source_signal)-1)'/fs;
    
    output_duration = max(duration, source_duration);
    t_r = (0:1/fs:output_duration-1/fs)';
    num_samples = length(t_r);


    
    fs_over = fs * oversample_ratio;
    
    % attenuation_beta = 0.1;
    % attenuation_gamma = 0.001;

    %% Pre-compute Reflection Paths
    fprintf('Pre-computing reflection paths...\n');
    reflection_paths_start = tic;
    
    reflection_paths = struct();
    path_counter = 1;

    if strcmp(room_type, 'rectangular')
        for nx = -max_reflection_order:max_reflection_order
            for ny = -max_reflection_order:max_reflection_order
                for nz = -max_reflection_order:max_reflection_order
                    if nx == 0 && ny == 0 && nz == 0
                        continue;
                    end
                    
                    if sum(abs([nx ny nz])) > max_reflection_order
                        continue;
                    end
                    
                    % ENHANCED: Calculate reflection coefficients for ALL frequencies
                    total_reflection_all_freq = zeros(1, num_freq_bands);
                    
                    for freq_idx = 1:num_freq_bands
                        % Calculate for each frequency band using proper wall assignments
                        if nx < 0
                            refl_x = beta_coefficients.x1(freq_idx)^abs(nx);
                        elseif nx > 0
                            refl_x = beta_coefficients.x2(freq_idx)^abs(nx);
                        else
                            refl_x = 1.0; % No reflection in this direction
                        end
                        
                        if ny < 0
                            refl_y = beta_coefficients.y1(freq_idx)^abs(ny);
                        elseif ny > 0
                            refl_y = beta_coefficients.y2(freq_idx)^abs(ny);
                        else
                            refl_y = 1.0; % No reflection in this direction
                        end
                        
                        if nz < 0
                            refl_z = beta_coefficients.z1(freq_idx)^abs(nz);
                        elseif nz > 0
                            refl_z = beta_coefficients.z2(freq_idx)^abs(nz);
                        else
                            refl_z = 1.0; % No reflection in this direction
                        end
                        
                        total_reflection_all_freq(freq_idx) = refl_x * refl_y * refl_z;
                    end
                    
                    % ENHANCED: Store comprehensive frequency-dependent information
                    reflection_paths(path_counter).indices = [nx, ny, nz];
                    reflection_paths(path_counter).coefficients_freq = total_reflection_all_freq; % Array of coefficients for all frequencies
                    reflection_paths(path_counter).frequencies = analysis_frequencies; % Frequency vector
                    reflection_paths(path_counter).offset = [2*nx*roomDim(1), 2*ny*roomDim(2), 2*nz*roomDim(3)];
                    reflection_paths(path_counter).type = 'rectangular';
                    
                    % ENHANCED: Also store single coefficient for backwards compatibility (use 1kHz)
                    [~, ref_freq_idx] = min(abs(analysis_frequencies - 1000));
                    reflection_paths(path_counter).coefficient = total_reflection_all_freq(ref_freq_idx);
                    
                    path_counter = path_counter + 1;
                end
            end
        end
    else
        num_segments = 16;
        [~, ref_freq_idx] = min(abs(analysis_frequencies - 1000));
        % Generate primary reflections from tunnel wall
        for seg = 1:num_segments
            angle = 2*pi * (seg-1) / num_segments;
            normal = [0, cos(angle), sin(angle)];
            
            reflection_paths(path_counter).indices = [0, seg, 0];
            reflection_paths(path_counter).coefficients_freq = beta_tunnel_wall; % ALL frequencies
            reflection_paths(path_counter).frequencies = analysis_frequencies;
            reflection_paths(path_counter).coefficient = beta_tunnel_wall(ref_freq_idx); % 1kHz reference
            reflection_paths(path_counter).normal = normal;
            reflection_paths(path_counter).type = 'cylinder_wall';
            path_counter = path_counter + 1;
        end
        
        % Add floor reflection
        reflection_paths(path_counter).indices = [0, 0, -1];
        reflection_paths(path_counter).coefficients_freq = beta_floor; % ALL frequencies
        reflection_paths(path_counter).frequencies = analysis_frequencies;
        reflection_paths(path_counter).coefficient = beta_floor(ref_freq_idx); % 1kHz reference
        reflection_paths(path_counter).normal = [0, 0, 1];
        reflection_paths(path_counter).type = 'cylinder_floor';
        path_counter = path_counter + 1;
        
        % Add entrance/exit reflections
        reflection_paths(path_counter).indices = [-1, 0, 0];
        reflection_paths(path_counter).coefficients_freq = beta_entrance; % ALL frequencies
        reflection_paths(path_counter).frequencies = analysis_frequencies;
        reflection_paths(path_counter).coefficient = beta_entrance(ref_freq_idx); % 1kHz reference
        reflection_paths(path_counter).normal = [1, 0, 0];
        reflection_paths(path_counter).type = 'cylinder_entrance';
        path_counter = path_counter + 1;
        
        reflection_paths(path_counter).indices = [1, 0, 0];
        reflection_paths(path_counter).coefficients_freq = beta_exit; % ALL frequencies
        reflection_paths(path_counter).frequencies = analysis_frequencies;
        reflection_paths(path_counter).coefficient = beta_exit(ref_freq_idx); % 1kHz reference
        reflection_paths(path_counter).normal = [-1, 0, 0];
        reflection_paths(path_counter).type = 'cylinder_exit';
        path_counter = path_counter + 1; 
    end

    reflection_paths_time = toc(reflection_paths_start);
    fprintf('Pre-computed %d reflection paths in %.2f seconds\n', length(reflection_paths), reflection_paths_time);

    %% Vehicle Body Setup
    if use_vehicle_body
        fprintf('Setting up vehicle body reflections...\n');
        
        vehicle_surfaces = struct();
        
        vehicle_surfaces(1).normal = [0, 0, 1];
        vehicle_surfaces(1).offset = vehicle_height/2;
        vehicle_surfaces(1).reflection = 0.8;
        
        vehicle_surfaces(2).normal = [0, 0, -1];
        vehicle_surfaces(2).offset = vehicle_height/2;
        vehicle_surfaces(2).reflection = 0.7;
        
        vehicle_surfaces(3).normal = [1, 0, 0];
        vehicle_surfaces(3).offset = vehicle_length/2;
        vehicle_surfaces(3).reflection = 0.8;
        
        vehicle_surfaces(4).normal = [-1, 0, 0];
        vehicle_surfaces(4).offset = vehicle_length/2;
        vehicle_surfaces(4).reflection = 0.8;
        
        vehicle_surfaces(5).normal = [0, 1, 0];
        vehicle_surfaces(5).offset = vehicle_width/2;
        vehicle_surfaces(5).reflection = 0.8;
        
        vehicle_surfaces(6).normal = [0, -1, 0];
        vehicle_surfaces(6).offset = vehicle_width/2;
        vehicle_surfaces(6).reflection = 0.8;
    end

    %% NOTAP Source-Time Dominant Processing
    fprintf('Starting Source-Time Dominant processing...\n');
    
    output_signals = zeros(numMics, num_samples);
    direct_signals = zeros(numMics, num_samples);
    reflection_signals = zeros(numMics, num_samples);
    vehicle_reflection_signals = zeros(numMics, num_samples);
    doppler_factors = zeros(numMics, num_samples);

    if save_mappings_for_localization
        localization_data = struct();
        localization_data.direct_path = struct('mic_index', {}, 'source_times', {}, ...
            'receiver_times', {}, 'source_positions', {}, 'distances', {}, 'doppler_factors', {});
        localization_data.reflections = struct('mic_index', {}, 'path_index', {}, 'indices', {}, ...
            'source_times', {}, 'receiver_times', {}, 'mirror_positions', {}, ...
            'distances', {}, 'reflection_coef', {}, 'doppler_factors', {});
    end

    %% Visualization Setup
    if enable_visualization
        fprintf('Setting up visualizations...\n');
        figure(1);
        hold on;
        
        if strcmp(room_type, 'rectangular')
            x = [0 0 0 0 roomDim(1) roomDim(1) roomDim(1) roomDim(1)];
            y = [0 0 roomDim(2) roomDim(2) 0 0 roomDim(2) roomDim(2)];
            z = [0 roomDim(3) 0 roomDim(3) 0 roomDim(3) 0 roomDim(3)];
            faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
            patch('Vertices', [x(:) y(:) z(:)], 'Faces', faces, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
        else
            tunnel_center = [tunnel_length/2, tunnel_radius, tunnel_radius];
            
            [X, Y, Z] = cylinder(tunnel_radius, 36);
            X = X * tunnel_length;
            Y = Y + tunnel_radius;
            Z = Z + tunnel_radius;
            surf(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', ':');
            
            theta = linspace(0, 2*pi, 37);
            x_circle = tunnel_radius * cos(theta) + tunnel_radius;
            y_circle = tunnel_radius * sin(theta) + tunnel_radius;
            
            plot3(zeros(size(theta)), x_circle, y_circle, 'k-', 'LineWidth', 1.5);
            plot3(repmat(tunnel_length, size(theta)), x_circle, y_circle, 'k-', 'LineWidth', 1.5);
            
            floor_x = [0, tunnel_length, tunnel_length, 0];
            floor_y = [0, 0, 2*tunnel_radius, 2*tunnel_radius];
            floor_z = zeros(1, 4);
            patch('Vertices', [floor_x(:) floor_y(:) floor_z(:)], 'Faces', [1 2 3 4], ...
                'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'k');
        end
        
        if is_moving
            t_plot = linspace(0, duration, 100);
            source_path = zeros(length(t_plot), 3);
            for i = 1:length(t_plot)
                source_path(i,:) = source_trajectory(t_plot(i));
            end
            plot3(source_path(:,1), source_path(:,2), source_path(:,3), 'r-', 'LineWidth', 2);
            scatter3(start_pos(1), start_pos(2), start_pos(3), 100, 'r', 'filled');
            text(start_pos(1)+0.1, start_pos(2), start_pos(3), 'Start', 'Color', 'r');
            scatter3(end_pos(1), end_pos(2), end_pos(3), 100, 'r', 'filled');
            text(end_pos(1)+0.1, end_pos(2), end_pos(3), 'End', 'Color', 'r');
        else
            scatter3(start_pos(1), start_pos(2), start_pos(3), 100, 'r', 'filled');
            text(start_pos(1)+0.1, start_pos(2), start_pos(3), 'Sound Source', 'Color', 'r');
        end
        
        if use_vehicle_body
            vehicle_pos = source_trajectory(0);
            vehicle_corners = [
                vehicle_pos(1) - vehicle_length/2, vehicle_pos(2) - vehicle_width/2, vehicle_pos(3) - vehicle_height/2;
                vehicle_pos(1) - vehicle_length/2, vehicle_pos(2) + vehicle_width/2, vehicle_pos(3) - vehicle_height/2;
                vehicle_pos(1) + vehicle_length/2, vehicle_pos(2) + vehicle_width/2, vehicle_pos(3) - vehicle_height/2;
                vehicle_pos(1) + vehicle_length/2, vehicle_pos(2) - vehicle_width/2, vehicle_pos(3) - vehicle_height/2;
                vehicle_pos(1) - vehicle_length/2, vehicle_pos(2) - vehicle_width/2, vehicle_pos(3) + vehicle_height/2;
                vehicle_pos(1) - vehicle_length/2, vehicle_pos(2) + vehicle_width/2, vehicle_pos(3) + vehicle_height/2;
                vehicle_pos(1) + vehicle_length/2, vehicle_pos(2) + vehicle_width/2, vehicle_pos(3) + vehicle_height/2;
                vehicle_pos(1) + vehicle_length/2, vehicle_pos(2) - vehicle_width/2, vehicle_pos(3) + vehicle_height/2
            ];
            
            vehicle_faces = [
                1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6; 5 6 7 8; 1 2 3 4
            ];
            
            patch('Vertices', vehicle_corners, 'Faces', vehicle_faces, ...
                'FaceColor', [0.7 0.7 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
        end
        
        scatter3(micPositions(:,1), micPositions(:,2), micPositions(:,3), 25, 'b', 'filled');
        text(arrayCenter(1), arrayCenter(2)-0.2, arrayCenter(3), 'Microphone Array', 'Color', 'b');
        
        xlabel('X (m)');
        ylabel('Y (m)');
        zlabel('Z (m)');
        title(['Room Geometry: ' strrep(room_type, '_', ' ') ' - Array Type: ' strrep(array_type, '_', ' ')]);
        grid on;
        axis equal;
        
        if strcmp(room_type, 'rectangular')
            axis([-0.5 roomDim(1)+0.5 -0.5 roomDim(2)+0.5 -0.5 roomDim(3)+0.5]);
        else
            axis([-0.5 tunnel_length+0.5 -0.5 2*tunnel_radius+0.5 -0.5 2*tunnel_radius+0.5]);
        end
        
        view(30, 20);
        drawnow;
    end

    %% Main Processing Loop
    fprintf('Starting Fast Hybrid Doppler-Reverb simulation...\n');
    total_start_time = tic;

    for mic = 1:numMics
        mic_pos = micPositions(mic, :);
        fprintf('Processing microphone %d/%d...\n', mic, numMics);
        
        %% SOURCE-TIME DOMINANT APPROACH
        fprintf('  Calculating transmission times using source-time dominant approach...\n');
        start_time = tic;
        
        source_positions = zeros(length(t_s), 3);
        for i = 1:length(t_s)
            source_positions(i,:) = source_trajectory(t_s(i));
        end
        
        distances = sqrt(sum((source_positions - repmat(mic_pos, length(t_s), 1)).^2, 2));
        t_x_direct = t_s + distances / c;
        
        if is_moving
            % dr = diff(distances);
            dt_x = diff(t_x_direct);
            
            f_s_x = 1 ./ dt_x;
            doppler_raw = fs ./ f_s_x;
            doppler_raw = [doppler_raw(1); doppler_raw];
            
            window_size = 11;
            doppler_smooth = smoothdata(doppler_raw, 'movmean', window_size);
        else
            doppler_smooth = ones(length(t_s), 1);
        end
        
        doppler_factors(mic, 1:length(doppler_smooth)) = doppler_smooth(1:min(length(doppler_smooth), size(doppler_factors, 2)));
        
        if save_mappings_for_localization
            direct_path_data = struct();
            direct_path_data.mic_index = mic;
            direct_path_data.source_times = t_s;
            direct_path_data.receiver_times = t_x_direct;
            direct_path_data.source_positions = source_positions;
            direct_path_data.distances = distances;
            direct_path_data.doppler_factors = doppler_smooth;
            
            localization_data.direct_path(mic) = direct_path_data;
        end
        
        %% Process Direct Path
        fprintf('  Processing direct path using NOTAP approach...\n');
        
        cutoff_freq = 0.8 * fs/2;
        [b_filter, a_filter] = butter(6, cutoff_freq/(fs_over/2));
        
        if is_moving
            t_over = min(t_x_direct):1/fs_over:max(t_x_direct);
            p_over = zeros(size(t_over));
            
            for i = 1:length(t_over)
                idx = find(t_x_direct <= t_over(i), 1, 'last');
                
                if ~isempty(idx) && idx <= length(source_signal)
                    distance = distances(idx);
                    attenuation = 1 / (4 * pi * distance^2) * exp(-0.01 * distance);
                    p_over(i) = source_signal(idx) * attenuation;
                end
            end
            
            p_filtered = filtfilt(b_filter, a_filter, p_over);
            p_direct = p_filtered(1:oversample_ratio:end);
            
        else
            distance = distances(1);
            delay_samples = round(distance/c * fs);
            
            p_direct = zeros(num_samples, 1);
            
            if delay_samples < num_samples
                attenuation = 1 / (4 * pi * distance^2) * exp(-0.01 * distance);
                
                signal_end = min(delay_samples + length(source_signal), num_samples);
                signal_length = signal_end - delay_samples;
                
                if signal_length > 0
                    p_direct(delay_samples+1:signal_end) = source_signal(1:signal_length) * attenuation;
                end
            end
        end
        
        if length(p_direct) > num_samples
            p_direct = p_direct(1:num_samples);
        else
            p_direct(end+1:num_samples) = 0;
        end
        
        direct_signals(mic, :) = p_direct;
        
        fprintf('  Direct path processing completed in %.2f seconds\n', toc(start_time));
        
        %% Process Reflections
        fprintf('  Processing reflections...\n');
        start_time = tic;
        
        if use_freq_domain_reflections
            fprintf('  Using efficient frequency-domain approach for reflections...\n');
            
            NFFT = 2^nextpow2(max(length(source_signal), num_samples));
            f = fs/2 * linspace(0, 1, NFFT/2+1);
            
            reflection_signal_freq = zeros(NFFT, 1);
            
            for path_idx = 1:length(reflection_paths)
                current_path = reflection_paths(path_idx);
                
                if mod(path_idx, 5) == 0 || path_idx == length(reflection_paths)
                    fprintf('    Processing reflection path %d/%d\n', path_idx, length(reflection_paths));
                end
                
                if strcmp(room_type, 'rectangular') || strcmp(current_path.type, 'rectangular')
                    mirror_offset = current_path.offset;
                    mirror_trajectory = @(t) source_trajectory(t) + mirror_offset;
                else
                    switch current_path.type
                        case 'cylinder_wall'
                            normal = current_path.normal;
                            mirror_trajectory = @(t) computeCylindricalMirrorPos(source_trajectory(t), normal, tunnel_radius, tunnel_center);
                        case 'cylinder_floor'
                            mirror_trajectory = @(t) reflect_across_floor(source_trajectory(t));
                        case 'cylinder_entrance'
                            mirror_trajectory = @(t) reflect_across_entrance(source_trajectory(t));
                        case 'cylinder_exit'
                            mirror_trajectory = @(t) reflect_across_exit(source_trajectory(t), tunnel_length);
                        otherwise
                            continue;
                    end
                end
                
                num_points = 20;
                sample_times = linspace(0, min(duration, source_duration), num_points);
                
                mirror_positions = zeros(num_points, 3);
                for i = 1:num_points
                    mirror_positions(i,:) = mirror_trajectory(sample_times(i));
                end
                
                mirror_distances = sqrt(sum((mirror_positions - repmat(mic_pos, num_points, 1)).^2, 2));
                delays = mirror_distances / c;
                % attenuations = current_path.coefficient ./ ...
                %     (1 + attenuation_beta * mirror_distances + attenuation_gamma * mirror_distances.^2);
                attenuations = current_path.coefficient ./ ...
                    (4 * pi * distance^2) * exp(-0.01 * distance);
                
                if is_moving
                    velocity_term = diff(mirror_distances) ./ diff(sample_times');
                    doppler_factors_mirror = c ./ (c - velocity_term);
                    doppler_factors_mirror = [doppler_factors_mirror(1); doppler_factors_mirror];
                else
                    doppler_factors_mirror = ones(num_points, 1);
                end
                
                segment_length = ceil(length(source_signal) / num_points);
                
                for seg = 1:num_points
                    seg_start = (seg-1) * segment_length + 1;
                    seg_end = min(seg * segment_length, length(source_signal));
                    
                    if seg_end >= seg_start
                        source_segment = source_signal(seg_start:seg_end);
                        source_segment_padded = [source_segment; zeros(NFFT - length(source_segment), 1)];
                        
                        source_freq = fft(source_segment_padded);
                        
                        doppler_factor = doppler_factors_mirror(seg);
                        delay = delays(seg);
                        attenuation = attenuations(seg);
                        
                        H_segment = attenuation * exp(-1j * 2 * pi * f * delay);
                        H_segment = [H_segment, conj(fliplr(H_segment(2:end-1)))];
                        
                        if is_moving && abs(doppler_factor - 1.0) > 0.01
                            source_freq_doppler = zeros(size(source_freq));
                            
                            for bin = 1:length(source_freq)
                                new_bin = round(bin * doppler_factor);
                                if new_bin >= 1 && new_bin <= length(source_freq)
                                    source_freq_doppler(new_bin) = source_freq_doppler(new_bin) + source_freq(bin);
                                end
                            end
                            
                            source_freq = source_freq_doppler;
                        end
                        
                        reflection_segment_freq = source_freq .* H_segment(:);
                        
                        time_position = ((seg-1) * segment_length) / length(source_signal);
                        phase_shift = exp(-1j * 2 * pi * (0:NFFT-1)' * time_position);
                        reflection_signal_freq = reflection_signal_freq + (reflection_segment_freq .* phase_shift);
                    end
                end
                
                if save_mappings_for_localization
                    reflection_path_data = struct();
                    reflection_path_data.mic_index = mic;
                    reflection_path_data.path_index = path_idx;
                    reflection_path_data.indices = current_path.indices;
                    reflection_path_data.source_times = sample_times;
                    reflection_path_data.receiver_times = sample_times + delays;
                    reflection_path_data.mirror_positions = mirror_positions;
                    reflection_path_data.distances = mirror_distances;
                    reflection_path_data.reflection_coef = current_path.coefficient;
                    reflection_path_data.doppler_factors = doppler_factors_mirror;
                    
                    localization_data.reflections(end+1) = reflection_path_data;
                end
            end
            
            reflection_signal = real(ifft(reflection_signal_freq));
            reflection_signal = reflection_signal(1:num_samples);
            reflection_signals(mic, :) = reflection_signal;
            
        else
            fprintf('  Using time-domain approach for reflections with accurate Doppler...\n');
            
            reflection_signal = zeros(num_samples, 1);
            
            for path_idx = 1:length(reflection_paths)
                current_path = reflection_paths(path_idx);
                
                if mod(path_idx, 5) == 0 || path_idx == length(reflection_paths)
                    fprintf('    Processing reflection path %d/%d\n', path_idx, length(reflection_paths));
                end
                
                if strcmp(room_type, 'rectangular') || strcmp(current_path.type, 'rectangular')
                    mirror_offset = current_path.offset;
                    mirror_trajectory = @(t) source_trajectory(t) + mirror_offset;
                else
                    switch current_path.type
                        case 'cylinder_wall'
                            normal = current_path.normal;
                            mirror_trajectory = @(t) computeCylindricalMirrorPos(source_trajectory(t), normal, tunnel_radius, tunnel_center);
                        case 'cylinder_floor'
                            mirror_trajectory = @(t) reflect_across_floor(source_trajectory(t));
                        case 'cylinder_entrance'
                            mirror_trajectory = @(t) reflect_across_entrance(source_trajectory(t));
                        case 'cylinder_exit'
                            mirror_trajectory = @(t) reflect_across_exit(source_trajectory(t), tunnel_length);
                        otherwise
                            continue;
                    end
                end
                
                mirror_positions = zeros(length(t_s), 3);
                for i = 1:length(t_s)
                    mirror_positions(i,:) = mirror_trajectory(t_s(i));
                end
                
                mirror_distances = sqrt(sum((mirror_positions - repmat(mic_pos, length(t_s), 1)).^2, 2));
                t_x_mirror = t_s + mirror_distances / c;
                
                if is_moving
                    mirror_dt_x = diff(t_x_mirror);
                    mirror_f_s_x = 1 ./ mirror_dt_x;
                    mirror_doppler_raw = fs ./ mirror_f_s_x;
                    mirror_doppler_raw = [mirror_doppler_raw(1); mirror_doppler_raw];
                    mirror_doppler_smooth = smoothdata(mirror_doppler_raw, 'movmean', window_size);
                else
                    mirror_doppler_smooth = ones(length(t_s), 1);
                end
                
                if save_mappings_for_localization
                    reflection_path_data = struct();
                    reflection_path_data.mic_index = mic;
                    reflection_path_data.path_index = path_idx;
                    reflection_path_data.indices = current_path.indices;
                    reflection_path_data.source_times = t_s;
                    reflection_path_data.receiver_times = t_x_mirror;
                    reflection_path_data.mirror_positions = mirror_positions;
                    reflection_path_data.distances = mirror_distances;
                    reflection_path_data.reflection_coef = current_path.coefficient;
                    reflection_path_data.doppler_factors = mirror_doppler_smooth;
                    
                    localization_data.reflections(end+1) = reflection_path_data;
                end
                
                chunk_size = 2000;
                t_min = min(t_x_mirror);
                t_max = max(t_x_mirror);
                
                mic_time_vector = (0:num_samples-1)'/fs;
                
                for chunk_start_idx = 1:chunk_size:num_samples
                    chunk_end_idx = min(chunk_start_idx + chunk_size - 1, num_samples);
                    chunk_times = mic_time_vector(chunk_start_idx:chunk_end_idx);
                    
                    for i = 1:length(chunk_times)
                        current_time = chunk_times(i);
                        idx = find(t_x_mirror <= current_time, 1, 'last');
                        
                        if ~isempty(idx) && idx <= length(source_signal)
                            distance = mirror_distances(idx);
                            attenuation = current_path.coefficient / (4 * pi * distance^2) * exp(-0.01 * distance);
                           
                            output_idx = chunk_start_idx + i - 1;
                            if output_idx <= num_samples
                                reflection_signal(output_idx) = reflection_signal(output_idx) + source_signal(idx) * attenuation;
                            end
                        end
                    end
                end
                
                reflection_signal = filtfilt(b_filter, a_filter, reflection_signal);
            end
            
            reflection_signals(mic, :) = reflection_signal;
        end
        
        fprintf('  Reflection processing completed in %.2f seconds\n', toc(start_time));
        
        %% Process Vehicle Body Reflections
        if use_vehicle_body
            fprintf('  Processing vehicle body reflections...\n');
            start_time = tic;
            
            vehicle_signal = zeros(num_samples, 1);
            
            for surface_idx = 1:length(vehicle_surfaces)
                current_surface = vehicle_surfaces(surface_idx);
                
                mirror_trajectory = @(t) calculateVehicleMirrorPosition(source_trajectory(t), vehicle_trajectory(t), current_surface);
                
                vehicle_mirror_positions = zeros(length(t_s), 3);
                for i = 1:length(t_s)
                    vehicle_mirror_positions(i,:) = mirror_trajectory(t_s(i));
                end
                
                vehicle_mirror_distances = sqrt(sum((vehicle_mirror_positions - repmat(mic_pos, length(t_s), 1)).^2, 2));
                t_x_vehicle_mirror = t_s + vehicle_mirror_distances / c;
                
                chunk_size = 2000;
                mic_time_vector = (0:num_samples-1)'/fs;
                
                for chunk_start_idx = 1:chunk_size:num_samples
                    chunk_end_idx = min(chunk_start_idx + chunk_size - 1, num_samples);
                    chunk_times = mic_time_vector(chunk_start_idx:chunk_end_idx);
                    
                    for i = 1:length(chunk_times)
                        current_time = chunk_times(i);
                        idx = find(t_x_vehicle_mirror <= current_time, 1, 'last');
                        
                        if ~isempty(idx) && idx <= length(source_signal)
                            distance = vehicle_mirror_distances(idx);
                            attenuation = current_surface.reflection / (4 * pi * distance^2) * exp(-0.01 * distance);
                             
                            output_idx = chunk_start_idx + i - 1;
                            if output_idx <= num_samples
                                vehicle_signal(output_idx) = vehicle_signal(output_idx) + source_signal(idx) * attenuation;
                            end
                        end
                    end
                end
            end
            
            vehicle_signal = filtfilt(b_filter, a_filter, vehicle_signal);
            vehicle_reflection_signals(mic, :) = vehicle_signal;
            
            fprintf('  Vehicle reflection processing completed in %.2f seconds\n', toc(start_time));
        end
        
        %% Combine All Components
        if use_vehicle_body
            output_signals(mic, :) = direct_signals(mic, :) + reflection_signals(mic, :) + vehicle_reflection_signals(mic, :);
        else
            output_signals(mic, :) = direct_signals(mic, :) + reflection_signals(mic, :);
        end
        
        fprintf('  Completed microphone %d/%d\n', mic, numMics);
    end

    total_processing_time = toc(total_start_time);
    fprintf('Processing completed in %.2f seconds\n', total_processing_time);

    %% Add Noise
    if enable_noise
        SNR_dB = 10;
        
        fprintf('Adding %s noise at SNR = %d dB...\n', noise_type, SNR_dB);
        
        noisy_signals = zeros(size(output_signals));
        
        switch noise_type
            case 'white'
                for mic = 1:numMics
                    signal_power = mean(output_signals(mic,:).^2);
                    noise_power = signal_power / (10^(SNR_dB/10));
                    noise = sqrt(noise_power) * randn(1, size(output_signals, 2));
                    noisy_signals(mic, :) = output_signals(mic, :) + noise;
                end
                
            case 'pink'
                for mic = 1:numMics
                    signal_power = mean(output_signals(mic,:).^2);
                    noise_power = signal_power / (10^(SNR_dB/10));
                    
                    white_noise = randn(1, size(output_signals, 2) + 1000);
                    
                    b = [0.049922035, -0.095993537, 0.050612699, -0.004408786];
                    a = [1, -2.494956002, 2.017265875, -0.522189400];
                    pink_noise = filter(b, a, white_noise);
                    
                    pink_noise = pink_noise(1001:end);
                    pink_noise = pink_noise / std(pink_noise) * sqrt(noise_power);
                    
                    noisy_signals(mic, :) = output_signals(mic, :) + pink_noise;
                end
                
            case 'ambient'
                nfft = size(output_signals, 2) + 1000;
                base_ambient = randn(1, nfft);
                
                b = [0.1, 0.2, 0.3, 0.2, 0.1];
                a = 1;
                base_ambient = filter(b, a, base_ambient);
                base_ambient = base_ambient(1001:end);
                
                for mic = 1:numMics
                    signal_power = mean(output_signals(mic,:).^2);
                    noise_power = signal_power / (10^(SNR_dB/10));
                    
                    mic_specific = randn(1, size(output_signals, 2));
                    mixed_noise = 0.7 * base_ambient + 0.3 * mic_specific;
                    mixed_noise = mixed_noise / std(mixed_noise) * sqrt(noise_power);
                    
                    noisy_signals(mic, :) = output_signals(mic, :) + mixed_noise;
                end
        end
        
        final_signals = noisy_signals;
        report_noise_info = struct('enabled', true, 'type', noise_type, 'snr_db', SNR_dB);
    else
        final_signals = output_signals;
        report_noise_info = struct('enabled', false, 'type', 'none', 'snr_db', inf);
    end

    %% Extract Impulse Response via Deconvolution
    fprintf('Extracting impulse responses via deconvolution...\n');
    
    simulated_impulse_responses = zeros(size(final_signals));
    
    for mic = 1:numMics
        ir_temp = deconvolveImpulseResponse(final_signals(mic, :), source_signal, fs, signal_choice, false, reg);
        ir_length = min(length(ir_temp), size(simulated_impulse_responses, 2));
        simulated_impulse_responses(mic, 1:ir_length) = ir_temp(1:ir_length);
        if ir_length < size(simulated_impulse_responses, 2)
            simulated_impulse_responses(mic, ir_length+1:end) = 0;
        end
    end
    
    fprintf('✓ Deconvolution completed for %d microphones\n', numMics);
    
    max_ir_length = round(fs * 2);
    if size(simulated_impulse_responses, 2) > max_ir_length
        simulated_impulse_responses = simulated_impulse_responses(:, 1:max_ir_length);
    end

    %% Visualization: Signal Results
    if enable_visualization
        if mod(numMics, 2) == 1
            center_mic = ceil(numMics/2);
        else
            center_mic = numMics/2;
        end
        
        figure(2);
        subplot(3,1,1);
        plot((0:length(direct_signals(center_mic,:))-1)/fs, direct_signals(center_mic,:));
        title('Direct Path Component at Center Microphone');
        xlabel('Time (s)'); ylabel('Amplitude'); grid on;
        
        subplot(3,1,2);
        plot((0:length(reflection_signals(center_mic,:))-1)/fs, reflection_signals(center_mic,:));
        title('Reflection Components at Center Microphone');
        xlabel('Time (s)'); ylabel('Amplitude'); grid on;
        
        subplot(3,1,3);
        plot((0:length(final_signals(center_mic,:))-1)/fs, final_signals(center_mic,:));
        title('Combined Signal at Center Microphone');
        xlabel('Time (s)'); ylabel('Amplitude'); grid on;
        
        figure(3);
        window_size = 1024;
        overlap = round(window_size*0.75);
        spectrogram(final_signals(center_mic,:), window_size, overlap, [], fs, 'yaxis');
        title(sprintf('Spectrogram of Signal at Center Microphone %d', center_mic));
        colorbar;
        
        figure(4);
        plot((0:length(doppler_factors(center_mic,:))-1)/fs, doppler_factors(center_mic,:));
        title('Doppler Factor at Center Microphone');
        xlabel('Time (s)'); ylabel('Doppler Factor'); grid on;
        ylim([0.8 1.2]);
        
        figure(5);
        hold on;
        
        num_mics_to_display = min(5, numMics);
        mic_indices = round(linspace(1, numMics, num_mics_to_display));
        
        for i = 1:length(mic_indices)
            mic_idx = mic_indices(i);
            plot((0:length(doppler_factors(mic_idx,:))-1)/fs, doppler_factors(mic_idx,:), ...
                'DisplayName', sprintf('Mic %d', mic_idx));
        end
        
        title('Doppler Factor Comparison Across Microphones');
        xlabel('Time (s)'); ylabel('Doppler Factor'); grid on;
        legend('Location', 'best'); ylim([0.8 1.2]);
        
        figure(6);
        hold on;
        
        for i = 1:length(mic_indices)
            mic_idx = mic_indices(i);
            window_len = 1000;
            signal = final_signals(mic_idx,:);
            rms_amplitude = zeros(1, length(signal));
            
            for j = 1:length(signal)
                start_idx = max(1, j - window_len/2);
                end_idx = min(length(signal), j + window_len/2);
                window_data = signal(start_idx:end_idx);
                rms_amplitude(j) = sqrt(mean(window_data.^2));
            end
            
            rms_amplitude = smoothdata(rms_amplitude, 'movmean', 500);
            
            plot((0:length(rms_amplitude)-1)/fs, rms_amplitude, ...
                'DisplayName', sprintf('Mic %d', mic_idx));
        end
        
        title('Signal Amplitude Comparison Across Microphones');
        xlabel('Time (s)'); ylabel('RMS Amplitude'); grid on;
        legend('Location', 'best');
        
        if save_mappings_for_localization
            figure(7);
            
            subplot(2,1,1);
            plot(localization_data.direct_path(center_mic).source_times, ...
                 localization_data.direct_path(center_mic).receiver_times, 'b.-');
            title('Source-Receiver Time Mapping (Direct Path)');
            xlabel('Source Time (s)'); ylabel('Receiver Time (s)'); grid on;
            
            subplot(2,1,2);
            hold on;
            
            mic_reflections = find([localization_data.reflections.mic_index] == center_mic);
            
            if ~isempty(mic_reflections)
                for i = 1:min(3, length(mic_reflections))
                    refl_idx = mic_reflections(i);
                    refl_data = localization_data.reflections(refl_idx);
                
                    color_map = 'rgb';
                    color_char = color_map(mod(i-1, 3) + 1);
                    line_spec = ['.-', color_char];
                
                    plot(refl_data.source_times, refl_data.receiver_times, ...
                        line_spec, ...
                        'DisplayName', sprintf('Reflection [%d,%d,%d]', ...
                        refl_data.indices(1), refl_data.indices(2), refl_data.indices(3)));
                end
            end
            
            title('Source-Receiver Time Mapping (Reflections)');
            xlabel('Source Time (s)'); ylabel('Receiver Time (s)'); grid on;
            legend('Location', 'northwest');
            
            figure(8);
            
            subplot(2,1,1);
            plot(localization_data.direct_path(center_mic).receiver_times, ...
                 localization_data.direct_path(center_mic).distances, 'b.-');
            title('Distance vs. Receiver Time (Direct Path)');
            xlabel('Receiver Time (s)'); ylabel('Distance (m)'); grid on;
            
            subplot(2,1,2);
            hold on;
            
            if ~isempty(mic_reflections)
                for i = 1:min(3, length(mic_reflections))
                    refl_idx = mic_reflections(i);
                    refl_data = localization_data.reflections(refl_idx);
                
                    color_map = 'rgb';
                    color_char = color_map(mod(i-1, 3) + 1);
                    line_spec = ['.-', color_char];
                
                    plot(refl_data.receiver_times, refl_data.distances, ...
                        line_spec, ...
                        'DisplayName', sprintf('Reflection [%d,%d,%d]', ...
                        refl_data.indices(1), refl_data.indices(2), refl_data.indices(3)));
                end
            end
            
            title('Distance vs. Receiver Time (Reflections)');
            xlabel('Receiver Time (s)'); ylabel('Distance (m)'); grid on;
            legend('Location', 'northwest');
        end
        
        for fig_num = 1:8
            if isgraphics(fig_num)
                fig_filename = sprintf('figure_%d_%s_%s.png', fig_num, room_type, array_type);
                saveas(figure(fig_num), fig_filename);
                fprintf('Saved figure to %s\n', fig_filename);
            end
        end
    end

    %% Save Results with Organized Structure
    fprintf('Saving results with organized file structure...\n');
    
    % Create organized output directories
    output_structure = createOutputDirectories(processed_params);
    
    % Generate consistent filename prefix
    config_prefix = sprintf('%s_%s_%s_mics%d', room_type, array_type, noise_type, numMics);
    
    % Save microphone signals to WAV files
    wav_files = {};
    for mic = 1:numMics
        signal_normalized = 0.9 * final_signals(mic,:) / max(abs(final_signals(mic,:)));
        filename = fullfile(output_structure.audio_dir, sprintf('%s_mic_%d_signal.wav', config_prefix, mic));
        audiowrite(filename, signal_normalized, fs);
        wav_files{end+1} = filename;
    end
    
    % Save simulation data
    simulation_file = fullfile(output_structure.data_dir, sprintf('%s_simulation_data.mat', config_prefix));
    save(simulation_file, 'final_signals', 'direct_signals', 'reflection_signals', ...
         'micPositions', 'doppler_factors', 'fs', 't_r', 'is_moving', ...
         'source_trajectory', 'roomDim', 'total_processing_time', 'room_type', ...
         'array_type', 'noise_type', 'F_parameter', 'report_noise_info');
    
    % Save impulse responses
    ir_file = fullfile(output_structure.ir_dir, sprintf('%s_impulse_responses.mat', config_prefix));
    save(ir_file, 'simulated_impulse_responses', 'fs', 'micPositions');
    
    % Save localization data
    if save_mappings_for_localization
        loc_file = fullfile(output_structure.localization_dir, sprintf('%s_localization_data.mat', config_prefix));
        save(loc_file, 'localization_data');
    end
    
    % Save parameters for reproducibility
    params_file = fullfile(output_structure.params_dir, sprintf('%s_parameters.mat', config_prefix));
    save(params_file, 'processed_params');
    
    %% Organize Results Structure
    results.success = true;
    results.execution_time = toc(start_time);
    
    % File paths
    results.files = struct();
    results.files.audio_files = wav_files;
    results.files.simulation_data = simulation_file;
    results.files.impulse_responses = ir_file;
    results.files.parameters = params_file;
    if save_mappings_for_localization
        results.files.localization_data = loc_file;
    end
    
    % Output directories
    results.directories = output_structure;
    
    % Core simulation data
    results.data.signals = struct();
    results.data.signals.final = final_signals;
    results.data.signals.direct = direct_signals;
    results.data.signals.reflections = reflection_signals;
    if use_vehicle_body
        results.data.signals.vehicle_reflections = vehicle_reflection_signals;
    end
    results.data.signals.impulse_responses = simulated_impulse_responses;
    
    % Array and source information
    results.data.array = struct();
    results.data.array.positions = micPositions;
    results.data.array.center = arrayCenter;
    results.data.array.type = array_type;
    results.data.array.orientation = array_orientation;
    results.data.array.num_mics = numMics;
    results.data.array.f_parameter = F_parameter;
    
    results.data.source = struct();
    results.data.source.start_position = start_pos;
    results.data.source.end_position = end_pos;
    results.data.source.is_moving = is_moving;
    if is_moving
        results.data.source.velocity = v_vec;
        results.data.source.velocity_magnitude = norm(v_vec);
    end
    results.data.source.signal_type = signal_choice;
    results.data.source.power = source_power;
    results.data.source.duration = duration;
    
    % Room information
    results.data.room = struct();
    results.data.room.type = room_type;
    results.data.room.dimensions = roomDim;
    if strcmp(room_type, 'cylindrical')
        results.data.room.tunnel_length = tunnel_length;
        results.data.room.tunnel_radius = tunnel_radius;
    end
    results.data.room.temperature = processed_params.temperature;
    results.data.room.humidity = processed_params.humidity;
    results.data.room.speed_of_sound = c;
    
    % Simulation settings
    results.data.simulation = struct();
    results.data.simulation.sampling_rate = fs;
    results.data.simulation.num_samples = num_samples;
    results.data.simulation.processing_time = total_processing_time;
    results.data.simulation.use_freq_domain_reflections = use_freq_domain_reflections;
    results.data.simulation.max_reflection_order = max_reflection_order;
    results.data.simulation.oversample_ratio = oversample_ratio;
    results.data.simulation.num_reflection_paths = length(reflection_paths);
    
    % Noise information
    results.data.noise = report_noise_info;
    
    % Analysis data
    results.data.analysis = struct();
    results.data.analysis.doppler_factors = doppler_factors;
    if save_mappings_for_localization
        results.data.analysis.localization_data = localization_data;
    end
    
    % Quality metrics
    results.data.metrics = struct();
    results.data.metrics.signal_to_noise_ratio = SNR_dB;
    results.data.metrics.array_quality = F_parameter;
    results.data.metrics.processing_efficiency = num_samples / total_processing_time; % samples per second
    
    % Create summary report
    summary_file = fullfile(output_structure.reports_dir, sprintf('%s_summary_report.txt', config_prefix));
    createSummaryReport(summary_file, results);
    results.files.summary_report = summary_file;
    
    fprintf('\n=== SIMULATION COMPLETED SUCCESSFULLY ===\n');
    fprintf('Total execution time: %.2f seconds\n', results.execution_time);
    fprintf('Output directory: %s\n', output_structure.base_dir);
    fprintf('Files created: %d\n', length(results.files_created));
    fprintf('Summary report: %s\n', summary_file);
    
catch ME
    results.success = false;
    results.error_message = ME.message;
    results.execution_time = toc(start_time);
    
    fprintf('\n=== SIMULATION FAILED ===\n');
    fprintf('Error: %s\n', ME.message);
    fprintf('Execution time before failure: %.2f seconds\n', results.execution_time);
    
    % Save error report
    if exist('output_structure', 'var') && isfield(output_structure, 'reports_dir')
        error_file = fullfile(output_structure.reports_dir, 'error_report.txt');
        fid = fopen(error_file, 'w');
        if fid > 0
            fprintf(fid, 'Simulation Error Report\n');
            fprintf(fid, 'Time: %s\n', datestr(now));
            fprintf(fid, 'Error: %s\n', ME.message);
            fprintf(fid, 'Stack:\n%s\n', ME.getReport());
            fclose(fid);
            % results.files.error_report = error_file;
        end
    end
    
    rethrow(ME);
end

end

%% ============= HELPER FUNCTIONS =============

function [processed_params, validation_result] = processInputParameters(params)
% Process and validate input parameters from GUI
validation_result = struct('valid', true, 'message', '');

try
    % Initialize processed parameters with defaults
    processed_params = getDefaultParameters();
    
    % Override with input parameters
    if isfield(params, 'room_type')
        processed_params.room_type = validateRoomType(params.room_type);
    end
    
    if isfield(params, 'roomDim') && ~isempty(params.roomDim)
        processed_params.roomDim = validateRoomDimensions(params.roomDim);
    elseif isfield(params, 'room_length') && isfield(params, 'room_width') && isfield(params, 'room_height')
        processed_params.roomDim = [params.room_length, params.room_width, params.room_height];
    end
    
    % Handle cylindrical room parameters
    if strcmp(processed_params.room_type, 'cylindrical')
        if isfield(params, 'tunnel_length')
            processed_params.tunnel_length = max(1, params.tunnel_length);
        else
            processed_params.tunnel_length = 20;
        end
        
        if isfield(params, 'tunnel_radius')
            processed_params.tunnel_radius = max(0.5, params.tunnel_radius);
        else
            processed_params.tunnel_radius = 5;
        end
        
        processed_params.roomDim = [processed_params.tunnel_length, 2*processed_params.tunnel_radius, 2*processed_params.tunnel_radius];
    end
    
    % Temperature and humidity
    if isfield(params, 'temperature')
        processed_params.temperature = params.temperature;
        processed_params.c = 331.3 * sqrt(1 + params.temperature/273.15);
    end
    
    if isfield(params, 'humidity')
        processed_params.humidity = params.humidity;
    end
    
    % Array configuration
    if isfield(params, 'array_type')
        processed_params.array_type = validateArrayType(params.array_type);
    end
    
    if isfield(params, 'array_orientation')
        processed_params.array_orientation = validateArrayOrientation(params.array_orientation);
    end
    
    if isfield(params, 'numMics')
        processed_params.numMics = max(1, round(params.numMics));
    end
    
    % Enhanced array center parsing (matching enhanced_3.m logic)
    if isfield(params, 'arrayCenter') && length(params.arrayCenter) == 3
        processed_params.arrayCenter = params.arrayCenter;
        fprintf('Debug: Using direct arrayCenter = [%.2f, %.2f, %.2f]\n', ...
            params.arrayCenter(1), params.arrayCenter(2), params.arrayCenter(3));
    elseif isfield(params, 'array_center_x') && isfield(params, 'array_center_y') && isfield(params, 'array_center_z')
        processed_params.arrayCenter = [params.array_center_x, params.array_center_y, params.array_center_z];
        fprintf('Debug: Parsed arrayCenter from individual fields = [%.2f, %.2f, %.2f]\n', ...
            params.array_center_x, params.array_center_y, params.array_center_z);
    else
        % Fallback to default if parsing fails
        processed_params.arrayCenter = [0.5, 0.5, 0.5];
        fprintf('Debug: Using default arrayCenter = [0.5, 0.5, 0.5] (no valid input found)\n');
    end
    
    % Verification print
    fprintf('Debug: Final processed_params.arrayCenter = [%.2f, %.2f, %.2f]\n', ...
        processed_params.arrayCenter(1), processed_params.arrayCenter(2), processed_params.arrayCenter(3));
    
    % Source configuration
    if isfield(params, 'signal_choice')
        processed_params.signal_choice = validateSignalChoice(params.signal_choice);
    end
    
    % Enhanced position parsing (matching enhanced_3.m logic)
    if isfield(params, 'start_pos') && length(params.start_pos) == 3
        processed_params.start_pos = params.start_pos;
    elseif isfield(params, 'start_x') && isfield(params, 'start_y') && isfield(params, 'start_z')
        processed_params.start_pos = [params.start_x, params.start_y, params.start_z];
    else
        % Use default if no valid start position found
        processed_params.start_pos = [2, 2, 2];
        fprintf('Debug: Using default start_pos = [2, 2, 2] (no valid input found)\n');
    end
    
    if isfield(params, 'end_pos') && length(params.end_pos) == 3
        processed_params.end_pos = params.end_pos;
    elseif isfield(params, 'end_x') && isfield(params, 'end_y') && isfield(params, 'end_z')
        processed_params.end_pos = [params.end_x, params.end_y, params.end_z];
    else
        processed_params.end_pos = processed_params.start_pos; % Stationary by default
    end
    
    % Motion detection debug
    is_moving_debug = ~isequal(processed_params.start_pos, processed_params.end_pos);
    
    
    if isfield(params, 'duration')
        processed_params.duration = max(0.1, params.duration);
    end
    
    if isfield(params, 'source_power')
        processed_params.source_power = max(0.01, params.source_power);
    end
    
    % Simulation options
    if isfield(params, 'enable_noise')
        processed_params.enable_noise = logical(params.enable_noise);
    end
    
    if isfield(params, 'noise_type')
        processed_params.noise_type = validateNoiseType(params.noise_type);
    end
    
    if isfield(params, 'fs')
        processed_params.fs = max(8000, params.fs);
    end
    
    if isfield(params, 'max_reflection_order')
        processed_params.max_reflection_order = max(0, round(params.max_reflection_order));
    end
    
    % Boolean options
    boolean_fields = {'use_vehicle_body', 'enable_visualization', 'use_freq_domain_reflections', ...
                     'use_batch_processing', 'save_mappings_for_localization'};
    
    for i = 1:length(boolean_fields)
        field = boolean_fields{i};
        if isfield(params, field)
            processed_params.(field) = logical(params.(field));
        end
    end
    
    % Numeric options
    if isfield(params, 'batch_size')
        processed_params.batch_size = max(100, round(params.batch_size));
    end
    
    if isfield(params, 'oversample_ratio')
        processed_params.oversample_ratio = max(2, round(params.oversample_ratio));
    end
    
    if isfield(params, 'reg')
        processed_params.reg = max(1e-6, params.reg);
    end
    
    % Validate parameter combinations
    validateParameterCombinations(processed_params);
    
catch ME
    validation_result.valid = false;
    validation_result.message = ME.message;
end

end

function defaults = getDefaultParameters()
% Get default parameter values
defaults = struct();

% Room configuration
defaults.room_type = 'rectangular';
defaults.roomDim = [4, 4, 4];
defaults.tunnel_length = 20;
defaults.tunnel_radius = 5;
defaults.temperature = 20;
defaults.humidity = 50;
defaults.c = 343;

% Array configuration
defaults.array_type = 'linear';
defaults.array_orientation = 'y-axis';
defaults.numMics = 2;
defaults.arrayCenter = [0.5, 0.5, 0.5];

% Source configuration
defaults.signal_choice = 'measurement_sweep';
defaults.start_pos = [2, 2, 2];
defaults.end_pos = [2, 2, 2];
defaults.duration = 1;
defaults.source_power = 0.5;

% Simulation options
defaults.enable_noise = true;
defaults.noise_type = 'white';
defaults.use_vehicle_body = false;
defaults.enable_visualization = true;
defaults.use_freq_domain_reflections = false;
defaults.use_batch_processing = true;
defaults.max_reflection_order = 1;
defaults.batch_size = 1000;
defaults.oversample_ratio = 8;
defaults.save_mappings_for_localization = true;
defaults.reg = 0.001;
defaults.fs = 16000;

end

function room_type = validateRoomType(input)
% Validate and normalize room type
valid_types = {'rectangular', 'cylindrical'};
input_lower = lower(input);
if any(strcmp(input_lower, valid_types))
    room_type = input_lower;
else
    error('Invalid room type: %s. Valid options: %s', input, strjoin(valid_types, ', '));
end
end

function array_type = validateArrayType(input)
% Validate and normalize array type
valid_types = {'linear', 'rectangular_grid', 'circular', 'spiral', 'spherical'};
input_processed = lower(strrep(input, ' ', '_'));
if any(strcmp(input_processed, valid_types))
    array_type = input_processed;
else
    error('Invalid array type: %s. Valid options: %s', input, strjoin(valid_types, ', '));
end
end

function orientation = validateArrayOrientation(input)
% Validate and normalize array orientation
valid_orientations = {'x-axis', 'y-axis', 'z-axis'};
input_lower = lower(input);
if any(strcmp(input_lower, valid_orientations))
    orientation = input_lower;
else
    error('Invalid array orientation: %s. Valid options: %s', input, strjoin(valid_orientations, ', '));
end
end

function signal_choice = validateSignalChoice(input)
% Validate and normalize signal choice
valid_signals = {'siren', 'chirp', 'measurement_sweep'};
input_processed = lower(strrep(input, ' ', '_'));
if any(strcmp(input_processed, valid_signals))
    signal_choice = input_processed;
else
    error('Invalid signal choice: %s. Valid options: %s', input, strjoin(valid_signals, ', '));
end
end

function noise_type = validateNoiseType(input)
% Validate and normalize noise type
valid_types = {'white', 'pink', 'ambient'};
input_lower = lower(input);
if any(strcmp(input_lower, valid_types))
    noise_type = input_lower;
else
    error('Invalid noise type: %s. Valid options: %s', input, strjoin(valid_types, ', '));
end
end

function dimensions = validateRoomDimensions(input)
% Validate room dimensions
if length(input) ~= 3
    error('Room dimensions must be a 3-element vector [length, width, height]');
end

if any(input <= 0)
    error('Room dimensions must be positive');
end

if any(input > 100)
    error('Room dimensions seem unrealistic (>100m). Please check input.');
end

dimensions = input;
end

function validateParameterCombinations(params)
% Validate parameter combinations make sense

% Check source positions are within room bounds
if strcmp(params.room_type, 'rectangular')
    if any(params.start_pos < 0) || any(params.start_pos > params.roomDim)
        warning('Start position may be outside room boundaries');
    end
    if any(params.end_pos < 0) || any(params.end_pos > params.roomDim)
        warning('End position may be outside room boundaries');
    end
end

% Check array center is reasonable
if strcmp(params.room_type, 'rectangular')
    if any(params.arrayCenter < 0) || any(params.arrayCenter > params.roomDim)
        warning('Array center may be outside room boundaries');
    end
end

% Check sampling rate vs signal type
if strcmp(params.signal_choice, 'measurement_sweep') && params.fs < 16000
    warning('Low sampling rate for measurement sweep may affect quality');
end

end

function output_structure = createOutputDirectories(params)
% Create organized output directory structure

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
config_name = sprintf('%s_%s_mics%d', params.room_type, params.array_type, params.numMics);
base_dir = sprintf('simulation_results_%s_%s', config_name, timestamp);

% Create directory structure
output_structure = struct();
output_structure.base_dir = base_dir;
output_structure.audio_dir = fullfile(base_dir, 'audio_signals');
output_structure.data_dir = fullfile(base_dir, 'simulation_data');
output_structure.ir_dir = fullfile(base_dir, 'impulse_responses');
output_structure.localization_dir = fullfile(base_dir, 'localization_data');
output_structure.params_dir = fullfile(base_dir, 'parameters');
output_structure.reports_dir = fullfile(base_dir, 'reports');
output_structure.figures_dir = fullfile(base_dir, 'figures');

% Create directories
dirs = struct2cell(output_structure);
for i = 1:length(dirs)
    if ~exist(dirs{i}, 'dir')
        mkdir(dirs{i});
    end
end

fprintf('Created output directory structure: %s\n', base_dir);

end

function createSummaryReport(filename, results)
% Create a comprehensive summary report

fid = fopen(filename, 'w');
if fid < 0
    warning('Could not create summary report file: %s', filename);
    return;
end

try
    fprintf(fid, '=== ACOUSTIC SIMULATION SUMMARY REPORT ===\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));
    
    fprintf(fid, 'EXECUTION STATUS:\n');
    fprintf(fid, '  Success: %s\n', string(results.success));
    fprintf(fid, '  Execution Time: %.2f seconds\n\n', results.execution_time);
    
    fprintf(fid, 'ROOM CONFIGURATION:\n');
    fprintf(fid, '  Type: %s\n', results.data.room.type);
    fprintf(fid, '  Dimensions: [%.1f, %.1f, %.1f] m\n', results.data.room.dimensions(1), results.data.room.dimensions(2), results.data.room.dimensions(3));
    fprintf(fid, '  Temperature: %.1f°C\n', results.data.room.temperature);
    fprintf(fid, '  Humidity: %.1f%%\n', results.data.room.humidity);
    fprintf(fid, '  Speed of Sound: %.1f m/s\n\n', results.data.room.speed_of_sound);
    
    fprintf(fid, 'MICROPHONE ARRAY:\n');
    fprintf(fid, '  Type: %s\n', results.data.array.type);
    fprintf(fid, '  Orientation: %s\n', results.data.array.orientation);
    fprintf(fid, '  Number of Microphones: %d\n', results.data.array.num_mics);
    fprintf(fid, '  Center Position: [%.2f, %.2f, %.2f] m\n', results.data.array.center(1), results.data.array.center(2), results.data.array.center(3));
    fprintf(fid, '  F-parameter (Quality): %.4f\n\n', results.data.array.f_parameter);
    
    fprintf(fid, 'SOURCE CONFIGURATION:\n');
    fprintf(fid, '  Signal Type: %s\n', results.data.source.signal_type);
    fprintf(fid, '  Power: %.3f W\n', results.data.source.power);
    fprintf(fid, '  Duration: %.2f s\n', results.data.source.duration);
    fprintf(fid, '  Start Position: [%.2f, %.2f, %.2f] m\n', results.data.source.start_position(1), results.data.source.start_position(2), results.data.source.start_position(3));
    fprintf(fid, '  End Position: [%.2f, %.2f, %.2f] m\n', results.data.source.end_position(1), results.data.source.end_position(2), results.data.source.end_position(3));
    fprintf(fid, '  Moving: %s\n', string(results.data.source.is_moving));
    if results.data.source.is_moving
        fprintf(fid, '  Velocity: [%.2f, %.2f, %.2f] m/s (%.2f m/s magnitude)\n', results.data.source.velocity(1), results.data.source.velocity(2), results.data.source.velocity(3), results.data.source.velocity_magnitude);
    end
    fprintf(fid, '\n');
    
    fprintf(fid, 'SIMULATION SETTINGS:\n');
    fprintf(fid, '  Sampling Rate: %d Hz\n', results.data.simulation.sampling_rate);
    fprintf(fid, '  Number of Samples: %d\n', results.data.simulation.num_samples);
    fprintf(fid, '  Frequency Domain Reflections: %s\n', string(results.data.simulation.use_freq_domain_reflections));
    fprintf(fid, '  Maximum Reflection Order: %d\n', results.data.simulation.max_reflection_order);
    fprintf(fid, '  Number of Reflection Paths: %d\n', results.data.simulation.num_reflection_paths);
    fprintf(fid, '  Oversample Ratio: %d\n', results.data.simulation.oversample_ratio);
    fprintf(fid, '  Processing Time: %.2f s\n\n', results.data.simulation.processing_time);
    
    fprintf(fid, 'NOISE CONFIGURATION:\n');
    fprintf(fid, '  Enabled: %s\n', string(results.data.noise.enabled));
    if results.data.noise.enabled
        fprintf(fid, '  Type: %s\n', results.data.noise.type);
        fprintf(fid, '  SNR: %.1f dB\n', results.data.noise.snr_db);
    end
    fprintf(fid, '\n');
    
    fprintf(fid, 'QUALITY METRICS:\n');
    fprintf(fid, '  Array Quality (F-parameter): %.4f\n', results.data.metrics.array_quality);
    fprintf(fid, '  Signal-to-Noise Ratio: %.1f dB\n', results.data.metrics.signal_to_noise_ratio);
    fprintf(fid, '  Processing Efficiency: %.0f samples/second\n\n', results.data.metrics.processing_efficiency);
    
    fprintf(fid, 'OUTPUT FILES:\n');
    fprintf(fid, '  Audio Signals: %d files\n', length(results.files.audio_files));
    fprintf(fid, '  Simulation Data: %s\n', results.files.simulation_data);
    fprintf(fid, '  Impulse Responses: %s\n', results.files.impulse_responses);
    fprintf(fid, '  Parameters: %s\n', results.files.parameters);
    if isfield(results.files, 'localization_data')
        fprintf(fid, '  Localization Data: %s\n', results.files.localization_data);
    end
    fprintf(fid, '\n');
    
    fprintf(fid, 'DIRECTORIES:\n');
    fprintf(fid, '  Base: %s\n', results.directories.base_dir);
    fprintf(fid, '  Audio: %s\n', results.directories.audio_dir);
    fprintf(fid, '  Data: %s\n', results.directories.data_dir);
    fprintf(fid, '  Reports: %s\n', results.directories.reports_dir);
    
    fclose(fid);
    fprintf('Summary report created: %s\n', filename);
    
catch ME
    fclose(fid);
     warning('%s: %s', ME.identifier, ME.message);
end

end

%% ============= SIGNAL PROCESSING HELPER FUNCTIONS =============

function sweepSignal = generateExpSweep(fs, duration, fStart, fEnd)
% Generate exponential sine sweep
N = round(duration * fs);
t = (0:N-1)' / fs;

K = (duration * fStart) / log(fEnd/fStart);
L = duration / log(fEnd/fStart);

sweepSignal = sin(2 * pi * K * (exp(t/L) - 1));

fadeLen = round(0.01 * fs);
fadeIn = (0:fadeLen-1)' / fadeLen;
fadeOut = flipud(fadeIn);

sweepSignal(1:fadeLen) = sweepSignal(1:fadeLen) .* fadeIn;
sweepSignal(end-fadeLen+1:end) = sweepSignal(end-fadeLen+1:end) .* fadeOut;

sweepSignal = sweepSignal / max(abs(sweepSignal));
end

function ir = deconvolveImpulseResponse(recordedSignal, referenceSignal, fs, ~, ~, reg_value)
% Perform FFT-based deconvolution to extract impulse response

recordedSignal = recordedSignal(:);
referenceSignal = referenceSignal(:);

N = length(recordedSignal) + length(referenceSignal) - 1;
N = 2^nextpow2(N);

X = fft(referenceSignal, N);
Y = fft(recordedSignal, N);

reg = reg_value;
H = Y .* conj(X) ./ (X .* conj(X) + reg * max(X .* conj(X)));

ir = real(ifft(H));

irLength = min(length(ir), round(fs * 2));
ir = ir(1:irLength);

fadeLen = round(0.05 * fs);
if fadeLen < length(ir)
    fadeOut = ones(size(ir));
    fadeOut(end-fadeLen+1:end) = (fadeLen-1:-1:0)' / fadeLen;
    ir = ir .* fadeOut;
end

ir = ir / max(abs(ir));
end

%% ============= GEOMETRIC HELPER FUNCTIONS =============

function mirror_pos = computeCylindricalMirrorPos(source_pos, ~, tunnel_radius, tunnel_center)
% Compute mirror position for cylindrical wall reflections
source_yz = source_pos(2:3);
center_yz = tunnel_center(2:3);

vec_to_source = source_yz - center_yz;
dir_to_source = vec_to_source / norm(vec_to_source);

point_on_cylinder = center_yz + tunnel_radius * dir_to_source;
reflection_point_yz = 2 * point_on_cylinder - source_yz;

mirror_pos = [source_pos(1), reflection_point_yz(1), reflection_point_yz(2)];
end

function reflected = reflect_across_floor(pos)
% Mirror position across z=0 plane
reflected = [pos(1), pos(2), -pos(3)];
end

function reflected = reflect_across_entrance(pos)
% Mirror position across x=0 plane
reflected = [-pos(1), pos(2), pos(3)];
end

function reflected = reflect_across_exit(pos, tunnel_length)
% Mirror position across x=tunnel_length plane
reflected = [2*tunnel_length - pos(1), pos(2), pos(3)];
end

function mirror_pos = calculateVehicleMirrorPosition(source_pos, vehicle_pos, surface)
% Calculate mirror position across vehicle surface
normal = surface.normal;
p0 = vehicle_pos + normal * surface.offset;
dist = dot(normal, source_pos - p0);
mirror_pos = source_pos - 2 * dist * normal;
end