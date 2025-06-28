%% Enhanced Fast Hybrid Doppler-Reverb Simulator with Acoustic Beamforming
% This script implements a hybrid approach that combines:
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

%%
% all the parts needed to be added to the code will be comments starting
% with _ 
%%

%%
% tests conductedand geneerate .wav and graphs for the following;
% isolate the ecoes and lsiten to them only and issue 
% isolate the doppler 
% isolate the noises 
% different mic geometry
% different room geometry
% different number od mics 
% call the localization mecahnizm
% call the beamforming and the ITD and ILD 
% we need to add a graph that shows the mirrors 
% the spiral shape is not a spiral so it needs to be adjusted
% the chapes needs to be handeled in a more generic way , so with different
% number of mics it can position itself correctly

%% Check for external call and prepare
if exist('EXTERNAL_CALL_FLAG', 'var') && EXTERNAL_CALL_FLAG
    fprintf('Simulator: Received external call, cleaning workspace...\n');
    clear;
    close all;
    clc;
    % Send "ready" signal by throwing a specific error
    error('Simulator_check_completed'); % This is our "handshake" signal
end

%% Normal operation - check if external variables exist
if exist('external_start_pos', 'var') 
    % External call (second time) - external variables are set, don't clear
    fprintf('Simulator is Ready for normal external usage\n');
else
    % Independent run - clear workspace
    clear;
    close all;
    clc;
    fprintf('Simulator is Ready for internal usage\n');
end



%% Enhanced Configuration Settings with Interactive Input

% Determine if using predefined configuration or running interactively
if ~exist('use_custom_config', 'var')
    use_custom_config = false;
end


% ask user only if not in use_custom_config mode
if ~use_custom_config 
    fprintf('\n=== CONFIGURATION SETUP ===\n');
    user_input = strtrim(lower(input('Do you want to use custom configuration? (y/n) [default: n]: ', 's')));
    use_custom_config = any(strcmp(user_input, {'y', 'yes'}));
end


 
%% Default Settings
default_values = struct( ...
    'room_type', 'rectangular', ...
    'rect_roomDim', [4, 4, 4], ...
    'tunnel_length', 20, ...
    'tunnel_radius', 5, ...
    'c', 343, ...
    'fs', 16000, ...
    'duration', 10, ...
    'start_pos', [2, 2, 2], ...
    'end_pos', [2, 2, 2], ...
    'signal_choice', 'measurement_sweep', ...
    'source_power', 0.5, ...
    'array_type', 'linear', ...
    'array_orientation', 'y-axis', ...
    'noise_type', 'white', ...
    'enable_noise', false, ...
    'use_vehicle_body', false, ...
    'enable_visualization', true, ...
    'use_freq_domain_reflections', false, ...
    'optimized_rev_time_method', false, ...
    'use_batch_processing', true, ...
    'max_reflection_order', 3, ...
    'batch_size', 1000, ...
    'oversample_ratio', 10, ...
    'save_mappings_for_localization', true, ...
    'reg',0.001 ...
);

%% Check for external override of default values
% If external_use is true, override defaults with external values
if exist('external_use', 'var') && external_use
    % Override start and end positions with external values if they exist
    if exist('external_start_pos', 'var') && ~isempty(external_start_pos)
        default_values.start_pos = external_start_pos;
        fprintf('Overriding default start position with external: [%.2f, %.2f, %.2f]\n', ...
            external_start_pos(1), external_start_pos(2), external_start_pos(3));
    end
    
    if exist('external_end_pos', 'var') && ~isempty(external_end_pos)
        default_values.end_pos = external_end_pos;
        fprintf('Overriding default end position with external: [%.2f, %.2f, %.2f]\n', ...
            external_end_pos(1), external_end_pos(2), external_end_pos(3));
    end
    
    % Can also override other external parameters here if needed 
    % for example:
    if exist('external_numMics', 'var') && ~isempty(external_numMics)
        numMics_input = external_numMics;
        fprintf('Overriding default microphone count with external: %d\n', external_numMics);
    end
    
    if exist('external_arrayCenter', 'var') && ~isempty(external_arrayCenter)
        arrayCenter_input = external_arrayCenter;
        fprintf('Overriding default array center with external: [%.2f, %.2f, %.2f]\n', ...
            external_arrayCenter(1), external_arrayCenter(2), external_arrayCenter(3));
    end

    % Override noise settings with external values if they exist
    if exist('external_enable_noise', 'var') && ~isempty(external_enable_noise)
        default_values.enable_noise = external_enable_noise;
        fprintf('Overriding default noise enable with external: %s\n', string(external_enable_noise));
    end
    
    if exist('external_noise_type', 'var') && ~isempty(external_noise_type)
        default_values.noise_type = external_noise_type;
        fprintf('Overriding default noise type with external: %s\n', external_noise_type);
    end
end

%% Interactive Input - Single Configuration Block
if use_custom_config
    fprintf('\n--- Room Configuration ---\n');
    disp('Room type options: rectangular, cylindrical');
    room_type = get_user_input('Room type', default_values.room_type);

    fprintf('\n--- Array Configuration ---\n');
    disp('Array type options: rectangular_grid, linear, circular, spiral, spherical');
    array_type = get_user_input('Array type', default_values.array_type);

    disp('Array orientation options: x-axis, y-axis, z-axis');
    array_orientation = get_user_input('Array orientation', default_values.array_orientation);

    % Get number of microphones
    numMics_input = get_user_number('Number of microphones', 2);
    
    % Get array center position
    fprintf('Array center position:\n');
    user_input = input('Array center [x, y, z] in meters [default: [0.5, 0.5, 0.5]]: ', 's');
    if ~isempty(user_input)
        arrayCenter_input = str2num(user_input);
        if length(arrayCenter_input) ~= 3
            fprintf('Invalid input. Using default array center.\n');
            arrayCenter_input = [0.5, 0.5, 0.5];
        end
    else
        arrayCenter_input = [0.5, 0.5, 0.5];
    end

    fprintf('\n--- Source Signal Configuration ---\n');
    disp('Signal type options: siren, chirp, measurement_sweep');
    signal_choice = get_user_input('Signal type', default_values.signal_choice);
    
    source_power = get_user_number('Source acoustic power (W)', default_values.source_power);

    fprintf('\n--- Vehicle Body Configuration ---\n');
    use_vehicle_body = get_user_boolean('Use vehicle body reflections', default_values.use_vehicle_body);

    fprintf('\n--- Optimization Settings ---\n');
    enable_visualization = get_user_boolean('Enable visualization', default_values.enable_visualization);
    use_freq_domain_reflections = get_user_boolean('Use frequency domain reflections', default_values.use_freq_domain_reflections);
    use_batch_processing = get_user_boolean('Use batch processing', default_values.use_batch_processing);

    max_reflection_order = get_user_number('Maximum reflection order', default_values.max_reflection_order);
    batch_size = get_user_number('Batch size', default_values.batch_size);
    oversample_ratio = get_user_number('Oversample ratio', default_values.oversample_ratio);
    save_mappings_for_localization = get_user_boolean('Save mappings for localization', default_values.save_mappings_for_localization);
    reg = default_values.reg;  
    fprintf('\n--- Room Dimensions and Acoustic Properties ---\n');
    
    if strcmp(room_type, 'rectangular')
        fprintf('Current room type: Rectangular\n');
        user_input = input(sprintf('Room dimensions [L, W, H] in meters [default: [%.1f, %.1f, %.1f]]: ', ...
            default_values.rect_roomDim(1), default_values.rect_roomDim(2), default_values.rect_roomDim(3)), 's');
        if ~isempty(user_input)
            roomDim = str2num(user_input);
            if length(roomDim) ~= 3
                fprintf('Invalid input. Using default room dimensions.\n');
                roomDim = default_values.rect_roomDim;
            end
        else
            roomDim = default_values.rect_roomDim;
        end
        
    elseif strcmp(room_type, 'cylindrical')
        fprintf('Current room type: Cylindrical Tunnel\n');
        tunnel_length = get_user_number('Tunnel length in meters', default_values.tunnel_length);
        tunnel_radius = get_user_number('Tunnel radius in meters', default_values.tunnel_radius);
        roomDim = [tunnel_length, 2*tunnel_radius, 2*tunnel_radius]; % Equivalent bounding box
    end
    
    fprintf('\n--- Physical Constants ---\n');
    c = get_user_number('Speed of sound (m/s)', default_values.c);
    fs = get_user_number('Sampling frequency (Hz)', default_values.fs);
    
    fprintf('\n--- Source Trajectory Parameters ---\n');
    duration = get_user_number('Simulation duration (s)', default_values.duration);


    %% special case here to know whether custom internal or if custom external availble 
    
    %% Unified Position Handling - handles both internal and external usage
    % Check if we should use interactive mode or external values
    if exist('external_use', 'var') && external_use
        if exist('external_start_pos', 'var') && exist('external_end_pos', 'var') && ~isempty(external_start_pos) && ~isempty(external_end_pos)
            % External mode - use the values already set in default_values (overridden above)
            start_pos = external_start_pos;
            end_pos = external_end_pos;
            fprintf('Using external positions: Start=[%.2f, %.2f, %.2f], End=[%.2f, %.2f, %.2f]\n', ...
                start_pos(1), start_pos(2), start_pos(3), end_pos(1), end_pos(2), end_pos(3));
        end

        if exist('external_enable_noise', 'var') && exist('external_noise_type', 'var') && ~isempty(external_enable_noise) && ~isempty(external_noise_type)
            % External mode - use the values already set in default_values (overridden above)
            enable_noise = external_enable_noise;
            noise_type = external_noise_type;
            fprintf('Using external noise enabling: enable_noise=[%s], noise_type=[%s]\n', ...
                enable_noise, noise_type);
        end


    else
        % Interactive mode - ask user for positions
        user_input = input(sprintf('Starting position [x, y, z] in meters [default: [%.1f, %.1f, %.1f]]: ', ...
            default_values.start_pos(1), default_values.start_pos(2), default_values.start_pos(3)), 's');
        if ~isempty(user_input)
            start_pos = str2num(user_input);
            if length(start_pos) ~= 3
                fprintf('Invalid input. Using default starting position.\n');
                start_pos = default_values.start_pos;
            end
        else
            start_pos = default_values.start_pos;
        end
        
        user_input = input(sprintf('Ending position [x, y, z] in meters [default: [%.1f, %.1f, %.1f]]: ', ...
            default_values.end_pos(1), default_values.end_pos(2), default_values.end_pos(3)), 's');
        if ~isempty(user_input)
            end_pos = str2num(user_input);
            if length(end_pos) ~= 3
                fprintf('Invalid input. Using default ending position.\n');
                end_pos = default_values.end_pos;
            end
        else
            end_pos = default_values.end_pos;
        end

        user_input = input(sprintf('Enable noise (true/false) [default: %s]: ', string(default_values.enable_noise)), 's');

        if ~isempty(user_input)
            % Convert string input to logical
            if strcmpi(user_input, 'true')
                enable_noise = true;
            elseif strcmpi(user_input, 'false')
                enable_noise = false;
            else
                error('Invalid input. Please enter "true" or "false".');
            end
        else
            enable_noise = default_values.enable_noise;
        end

        
        if enable_noise
            user_input = input(sprintf('Noise type (white - pink - ambient) [default: %s]: ', default_values.noise_type), 's');
            if ~isempty(user_input)
                noise_type = user_input;
            else
                noise_type = default_values.noise_type;
            end
        else
            noise_type = default_values.noise_type;
        end
    end

    %% end special case here to know whether custom internal or if custom external availble 

    fprintf('\n=== CONFIGURATION COMPLETE ===\n');

else
    % Use all default values
    room_type = default_values.room_type;
    signal_choice = default_values.signal_choice;
    source_power = default_values.source_power;
    array_type = default_values.array_type;
    array_orientation = default_values.array_orientation;
    noise_type = default_values.noise_type;
    enable_noise  = default_values.enable_noise ;
    use_vehicle_body = default_values.use_vehicle_body;
    enable_visualization = default_values.enable_visualization;
    use_freq_domain_reflections = default_values.use_freq_domain_reflections;
    optimized_rev_time_method = default_values.optimized_rev_time_method;
    use_batch_processing = default_values.use_batch_processing;
    max_reflection_order = default_values.max_reflection_order;
    batch_size = default_values.batch_size;
    oversample_ratio = default_values.oversample_ratio;
    save_mappings_for_localization = default_values.save_mappings_for_localization;
    reg = default_values.reg;

    % Use default values for room dimensions and physical constants
    if strcmp(room_type, 'rectangular')
        roomDim = default_values.rect_roomDim;
    elseif strcmp(room_type, 'cylindrical')
        tunnel_length = default_values.tunnel_length;
        tunnel_radius = default_values.tunnel_radius;
        roomDim = [tunnel_length, 2*tunnel_radius, 2*tunnel_radius];
    end
    
    c = default_values.c;
    fs = default_values.fs;
    duration = default_values.duration;
    start_pos = default_values.start_pos;
    end_pos = default_values.end_pos;
end

%% Display Final Configuration
fprintf('\nCurrent Configuration:\n');
fprintf('  Room type: %s\n', room_type);
fprintf('  Array type: %s\n', array_type);
fprintf('  Array orientation: %s\n', array_orientation);
fprintf('  enable_noise: %s\n', string(enable_noise));
if enable_noise
    fprintf('  Noise type: %s\n', noise_type);
end
fprintf('  Vehicle body: %s\n', string(use_vehicle_body));
fprintf('  Visualization: %s\n', string(enable_visualization));
fprintf('  Freq domain reflections: %s\n', string(use_freq_domain_reflections));
fprintf('  Optimized reflection time method: %s\n', string(optimized_rev_time_method));
fprintf('  Batch processing: %s\n', string(use_batch_processing));
fprintf('  Max reflection order: %d\n', max_reflection_order);
fprintf('  Batch size: %d\n', batch_size);
fprintf('  Oversample ratio: %d\n', oversample_ratio);
fprintf('  Save localization mappings: %s\n', string(save_mappings_for_localization));
fprintf('\n');

%% input validation check

% Validate room boundaries for source positions
if strcmp(room_type, 'rectangular')
    % Check if positions are within room bounds
    if any(start_pos < 0) || any(start_pos > roomDim)
        fprintf('Warning: Starting position is outside room boundaries. Adjusting...\n');
        start_pos = max(0.1, min(start_pos, roomDim - 0.1));
    end
    if any(end_pos < 0) || any(end_pos > roomDim)
        fprintf('Warning: Ending position is outside room boundaries. Adjusting...\n');
        end_pos = max(0.1, min(end_pos, roomDim - 0.1));
    end
elseif strcmp(room_type, 'cylindrical')
    % Check if positions are within cylindrical tunnel bounds
    for pos_idx = 1:2
        if pos_idx == 1
            pos = start_pos;
        else
            pos = end_pos;
        end
        
        if pos(1) < 0 || pos(1) > tunnel_length
            fprintf('Warning: Position x-coordinate is outside tunnel length. Adjusting...\n');
            if pos_idx == 1
                start_pos(1) = max(0.1, min(start_pos(1), tunnel_length - 0.1));
            else
                end_pos(1) = max(0.1, min(end_pos(1), tunnel_length - 0.1));
            end
        end
        
        % Check if y,z coordinates are within circular cross-section
        dist_from_center = sqrt((pos(2) - tunnel_radius)^2 + (pos(3) - tunnel_radius)^2);
        if dist_from_center > tunnel_radius - 0.1
            fprintf('Warning: Position is outside tunnel radius. Adjusting...\n');
            % Move position towards center
            direction = [pos(2) - tunnel_radius, pos(3) - tunnel_radius] / dist_from_center;
            new_pos_yz = [tunnel_radius, tunnel_radius] + direction * (tunnel_radius - 0.1);
            if pos_idx == 1
                start_pos(2:3) = new_pos_yz;
            else
                end_pos(2:3) = new_pos_yz;
            end
        end
    end
end

% Display current parameters
fprintf('\nCurrent Parameters:\n');
if strcmp(room_type, 'rectangular')
    fprintf('  Room dimensions [L×W×H]: [%.1f, %.1f, %.1f] m\n', roomDim(1), roomDim(2), roomDim(3));
elseif strcmp(room_type, 'cylindrical')
    fprintf('  Tunnel length: %.1f m\n', tunnel_length);
    fprintf('  Tunnel radius: %.1f m\n', tunnel_radius);
    fprintf('  Equivalent bounding box: [%.1f, %.1f, %.1f] m\n', roomDim(1), roomDim(2), roomDim(3));
end
fprintf('  Speed of sound: %d m/s\n', c);
fprintf('  Sampling frequency: %d Hz\n', fs);
fprintf('  Simulation duration: %.1f s\n', duration);
fprintf('  Starting position: [%.1f, %.1f, %.1f] m\n', start_pos(1), start_pos(2), start_pos(3));
fprintf('  Ending position: [%.1f, %.1f, %.1f] m\n', end_pos(1), end_pos(2), end_pos(3));
if isequal(start_pos, end_pos)
    fprintf('  Source motion: Stationary\n');
else
    velocity = norm(end_pos - start_pos) / duration;
    fprintf('  Source motion: Moving (velocity: %.2f m/s)\n', velocity);
end
fprintf('\n');

is_moving = ~isequal(start_pos, end_pos); % Flag to check if source is moving

if is_moving
    v_vec = (end_pos - start_pos) / duration; % Velocity vector
    source_trajectory = @(t) start_pos + v_vec * t; % Linear motion function
    fprintf('Source is moving from [%.2f, %.2f, %.2f] to [%.2f, %.2f, %.2f]\n', ...
        start_pos(1), start_pos(2), start_pos(3), end_pos(1), end_pos(2), end_pos(3));
    fprintf('Source velocity: [%.2f, %.2f, %.2f] m/s\n', v_vec(1), v_vec(2), v_vec(3));
else
    source_trajectory = @(t) start_pos; % Stationary source
    fprintf('Source is stationary at position [%.2f, %.2f, %.2f]\n', ...
        start_pos(1), start_pos(2), start_pos(3));
end

% Vehicle body parameters (if enabled)
if use_vehicle_body
    vehicle_length = 4.5;  % Length of vehicle in meters
    vehicle_width = 1.8;   % Width of vehicle in meters
    vehicle_height = 1.5;  % Height of vehicle in meters
    
    % Vehicle body is centered at the source position
    vehicle_trajectory = @(t) source_trajectory(t); % Vehicle follows the source
end

%% Define Sensor Array

array_size = 0.3;      % Array size (m x m)
array_elements = 2;  % Number of elements in each direction (2x2 = 4 total)
element_spacing = array_size/(array_elements-1); % Spacing between array elements (m)

% Use input array center or default
if exist('arrayCenter_input', 'var') && ~isempty(arrayCenter_input)
    arrayCenter = arrayCenter_input;
else
    arrayCenter = [0.5, 0.5, 0.5]; % Fixed position in room coordinates [x, y, z]
end


% Number of microphones based on array type
% Use input values, with fallback to defaults based on array type
% Determine number of microphones - prioritize external input
if exist('external_numMics', 'var') && ~isempty(external_numMics)
    numMics = external_numMics;
    fprintf('Using external microphone count: %d\n', numMics);
elseif exist('numMics_input', 'var') && ~isempty(numMics_input)
    numMics = numMics_input;
    fprintf('Using input microphone count: %d\n', numMics);
else
    % Only use array type defaults if no external input provided
    switch array_type
        case 'single'  % NEW: Explicit single microphone array type
            numMics = 1;
        case 'rectangular_grid'
            numMics = array_elements^2;
        case 'linear'
            numMics = 4;
        case 'circular'
            numMics = 6;
        case 'spiral'
            numMics = 6;
        case 'spherical'
            numMics = 6;
            
        otherwise
            numMics = array_elements^2;
    end
    fprintf('Using default microphone count for %s: %d\n', array_type, numMics);
end

% Initialize microphone positions
micPositions = zeros(numMics, 3);

% Set microphone positions based on array type
switch array_type
    case 'single'  % NEW: Single microphone case
        micPositions(1, :) = arrayCenter;

    case 'rectangular_grid'
    % Standard rectangular grid with configurable orientation
    [x_grid, z_grid] = meshgrid(-(array_size/2):element_spacing:(array_size/2), ...
                               -(array_size/2):element_spacing:(array_size/2));
    
        switch array_orientation
            case 'x-axis'  % Array faces along X-axis (grid in YZ plane)
                array_x = repmat(arrayCenter(1), size(x_grid));
                array_y = arrayCenter(2) + x_grid;
                array_z = arrayCenter(3) + z_grid;
            case 'y-axis'  % Array faces along Y-axis (grid in XZ plane)
                array_x = arrayCenter(1) + x_grid;
                array_y = repmat(arrayCenter(2), size(x_grid));
                array_z = arrayCenter(3) + z_grid;
            case 'z-axis'  % Array faces along Z-axis (grid in XY plane)
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
        
    case 'linear'
        if numMics == 1
            % Single microphone at array center
            micPositions(1, :) = arrayCenter;
        else
            % Existing linear array code for multiple microphones
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
            % Existing circular array code...
            radius = array_size/2;
            angles = linspace(0, 2*pi, numMics+1);
            angles = angles(1:end-1); % Remove duplicate endpoint
            
            for i = 1:numMics
                switch array_orientation
                    case 'x-axis'  % Circle in YZ plane
                        micPositions(i, :) = [arrayCenter(1), ...
                                             arrayCenter(2) + radius*cos(angles(i)), ...
                                             arrayCenter(3) + radius*sin(angles(i))];
                    case 'y-axis'  % Circle in XZ plane
                        micPositions(i, :) = [arrayCenter(1) + radius*cos(angles(i)), ...
                                             arrayCenter(2), ...
                                             arrayCenter(3) + radius*sin(angles(i))];
                    case 'z-axis'  % Circle in XY plane
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
            % Spiral array in plane perpendicular to the selected axis
            radius = array_size/2;
            num_arms = 4; % Number of spiral arms
            points_per_arm = ceil(numMics / num_arms);
            actual_mics = num_arms * points_per_arm;
            
            mic_idx = 1;
            for arm = 1:num_arms
                phi0 = 2*pi*(arm-1)/num_arms; % Initial angle for this arm
                for pt = 1:points_per_arm
                    if mic_idx <= numMics
                        r = radius * sqrt(pt/points_per_arm);
                        phi = phi0 + 0.1*r; % Small offset for better distribution
                        
                        switch array_orientation
                            case 'x-axis'  % Spiral in YZ plane
                                micPositions(mic_idx, :) = [arrayCenter(1), ...
                                                          arrayCenter(2) + r*cos(phi), ...
                                                          arrayCenter(3) + r*sin(phi)];
                            case 'y-axis'  % Spiral in XZ plane
                                micPositions(mic_idx, :) = [arrayCenter(1) + r*cos(phi), ...
                                                          arrayCenter(2), ...
                                                          arrayCenter(3) + r*sin(phi)];
                            case 'z-axis'  % Spiral in XY plane
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
            % Spherical array (approximately uniform distribution on sphere)
            radius = array_size/2;
            
            % Generate points using Golden Spiral method
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

fprintf('Number of microphones: %d\n', numMics);
if numMics == 1
    F_parameter = 1.0; % Single microphone: perfect unicity (no redundancy)
    fprintf('Single microphone: F-parameter = 1.0 (perfect unicity)\n');
else
    % Calculate F-parameter (unicity metric) for the array
    % - Create co-array matrix of vector spacings
    co_array = zeros(numMics^2, 3);
    counter = 1;
    for i = 1:numMics
        for j = 1:numMics
            co_array(counter, :) = micPositions(i, :) - micPositions(j, :);
            counter = counter + 1;
        end
    end
    
    % - Find unique vector spacings 
    % Using a small tolerance to account for floating-point errors
    unique_vectors = 0;
    tolerance = 1e-10;
    is_unique = false(size(co_array, 1), 1);
    
    for i = 1:size(co_array, 1)
        if ~any(is_unique)
            is_unique(i) = true;
            unique_vectors = unique_vectors + 1;
        else
            % Check if this vector is unique compared to previously found unique vectors
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
    
    % - Calculate the maximum possible unique vectors
    max_unique_vectors = numMics^2 - (numMics - 1);
    
    % - Calculate F-parameter (unicity ratio)
    F_parameter = unique_vectors / max_unique_vectors;
    fprintf('Array F-parameter (unicity): %.4f\n', F_parameter);
end




%% FIXED: Proper Frequency-Dependent Reflection Processing
% Replace the section from line ~580 to ~950 in your original code

%% Define Acoustic Properties of Walls/Surfaces
% Frequencies for acoustic properties
freq = [125, 250, 500, 1000, 2000, 4000]; % Frequencies in Hz

if strcmp(room_type, 'rectangular')
    % Define absorption coefficients for each wall in rectangular room
    % [front, left, right, ceiling, ground, back]
    absorption_coef = zeros(6, length(freq));
    
    % Front and left side
    absorption_coef(1:2, :) = repmat([0.18, 0.06, 0.04, 0.03, 0.02, 0.02], 2, 1);
    
    % Right side and ceiling
    absorption_coef(3:4, :) = repmat([0.1, 0.05, 0.06, 0.07, 0.09, 0.08], 2, 1);
    
    % Ground
    absorption_coef(5, :) = [0.01, 0.01, 0.01, 0.01, 0.02, 0.02];
    
    % Back side
    absorption_coef(6, :) = [0.35, 0.25, 0.18, 0.12, 0.07, 0.04];
else  % Cylindrical tunnel
    % Define absorption coefficients for tunnel surfaces
    % [wall, entrance, exit, N/A, floor, N/A]
    absorption_coef = zeros(6, length(freq));
    
    % Brick tunnel walls - approximation
    absorption_coef(1, :) = [0.03, 0.03, 0.03, 0.04, 0.05, 0.07]; 
    
    % Tunnel entrance and exit (mostly reflective)
    absorption_coef(2:3, :) = repmat([0.01, 0.01, 0.01, 0.01, 0.02, 0.02], 2, 1);
    
    % Floor (asphalt-like)
    absorption_coef(5, :) = [0.02, 0.02, 0.03, 0.03, 0.05, 0.05];
end

% Calculate reflection coefficients using formula β = √(1-γ)
reflection_coef = sqrt(1 - absorption_coef);

% Analysis frequency (for demonstration we'll use 1000 Hz)
freq_analysis = 1000; % Hz
analysis_frequencies = [125, 250, 500, 1000, 2000, 4000];
num_freq_bands = length(analysis_frequencies);
[~, freq_idx] = min(abs(freq - freq_analysis));

% Get reflection coefficients at ALL frequencies (not just one)
if strcmp(room_type, 'rectangular')
    % FIXED: Create frequency-dependent beta coefficients structure
    beta_coefficients = struct();
    beta_coefficients.x1 = reflection_coef(1, :); % Front wall - ALL frequencies
    beta_coefficients.x2 = reflection_coef(3, :); % Right wall - ALL frequencies  
    beta_coefficients.y1 = reflection_coef(6, :); % Back wall - ALL frequencies
    beta_coefficients.y2 = reflection_coef(2, :); % Left wall - ALL frequencies
    beta_coefficients.z1 = reflection_coef(5, :); % Ground - ALL frequencies
    beta_coefficients.z2 = reflection_coef(4, :); % Ceiling - ALL frequencies
    
    fprintf('Reflection coefficients loaded for %d frequency bands\n', num_freq_bands);
    fprintf('Frequencies: %s Hz\n', mat2str(analysis_frequencies));
    
else  % Cylindrical tunnel
    % FIXED: Create frequency-dependent coefficients for cylindrical room
    beta_tunnel_wall = reflection_coef(1, :); % ALL frequencies
    beta_entrance = reflection_coef(2, :);    % ALL frequencies
    beta_exit = reflection_coef(3, :);        % ALL frequencies
    beta_floor = reflection_coef(5, :);       % ALL frequencies
    
    fprintf('Cylindrical tunnel reflection coefficients loaded for %d frequency bands\n', num_freq_bands);
end


%% Load or Generate Source Signal
% signal_choice = 'measurement_sweep';  % Options: 'siren', 'chirp', 'measurement_sweep'

switch signal_choice
    case 'siren'
        % Try to load siren sound
        try
            [source_signal, src_fs] = audioread('siren_sound/original.wav');
            
            % If stereo, convert to mono
            if size(source_signal, 2) > 1
                source_signal = mean(source_signal, 2);
            end
            
            % Resample if necessary
            if src_fs ~= fs
                source_signal = resample(source_signal, fs, src_fs);
            end
            
            fprintf('Loaded siren sound, length: %.2f seconds\n', length(source_signal)/fs);
        catch
            error('Siren sound file not found!');
        end
        
    case 'chirp'
        % Generate a chirp
        fprintf('Generating chirp signal...\n');
        Tsig = 0.5; % Duration of signal in seconds
        t_sig = 0:1/fs:Tsig-1/fs;
        source_signal = chirp(t_sig, 500, Tsig, 2000)';
        
    case 'measurement_sweep'
        % Generate same exponential sweep as measurements
        fprintf('Generating exponential sweep matching measurements...\n');
        signalDuration = duration;  % Same as measurement
        freqStart = 20;       % Same as measurement
        freqEnd = 20000;      % Same as measurement
        source_signal = generateExpSweep(fs, signalDuration, freqStart, freqEnd);
        
    otherwise
        error('Unknown signal choice: %s', signal_choice);
end

% Scale the source signal to a reasonable amplitude
% source_power = 0.5;  % Acoustic power in W
reference_distance = 1.0; % meters
scaling_factor = sqrt(source_power/(4*pi*reference_distance^2));
source_signal = source_signal * scaling_factor;

% % Standalone visualization figure
% figure('Position', [100, 100, 800, 600]);
% 
% % Time vector for plotting
% t_plot = (0:length(source_signal)-1) / fs;
% 
% % Plot time domain
% subplot(2,1,1);
% plot(t_plot, source_signal);
% title(sprintf('Source Signal at 1m Reference Distance (%s)', signal_choice));
% xlabel('Time (s)');
% ylabel('Amplitude (Pa)');
% grid on;
% 
% % Plot frequency domain
% subplot(2,1,2);
% [Pxx, f] = pwelch(source_signal, [], [], [], fs);
% semilogx(f, 10*log10(Pxx));
% title('Source Signal Spectrum');
% xlabel('Frequency (Hz)');
% ylabel('Power (dB)');
% grid on;
% xlim([20 fs/2]);
% 
% fprintf('Source signal peak amplitude at 1m: %.4f Pa\n', max(abs(source_signal)));
% fprintf('Source signal RMS at 1m: %.4f Pa\n', rms(source_signal));
% 
% drawnow;

%% simulator timing handeling 

% Source time vector (SOURCE-TIME DOMINANT APPROACH)
source_duration = length(source_signal)/fs;  
t_s = (0:length(source_signal)-1)'/fs;   % |__ts__|__ts__|__ts__ sampling time of source singal

% Determine output duration and create receiver time vector
% maximum of the duration specified by the user and the signal source 
% if duration > source_duration then the dration - source_duration will be
% the reflection and noise 

output_duration = max(duration, source_duration);   
t_r = (0:1/fs:output_duration-1/fs)';   % 0 second|__tr__|__tr__|__tr__|one second
num_samples = length(t_r);

%% Enhanced Filter Parameters

fs_over = fs * oversample_ratio;  % Oversampled frequency

%% Calculate Reflection Paths (Pre-computation) - COMPLETELY REWRITTEN
fprintf('Pre-computing reflection paths with frequency-dependent coefficients...\n');
tic;

% Initialize structure to store reflection paths
reflection_paths = struct();
path_counter = 1;

if strcmp(room_type, 'rectangular')
    % Generate reflection paths for rectangular room
    for nx = -max_reflection_order:max_reflection_order
        for ny = -max_reflection_order:max_reflection_order
            for nz = -max_reflection_order:max_reflection_order
                if nx == 0 && ny == 0 && nz == 0, continue; end
                if sum(abs([nx ny nz])) > max_reflection_order, continue; end
                
                % FIXED: Calculate reflection coefficients for ALL frequencies
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
                
                % FIXED: Store comprehensive frequency-dependent information
                reflection_paths(path_counter).indices = [nx, ny, nz];
                reflection_paths(path_counter).coefficients_freq = total_reflection_all_freq; % Array of coefficients for all frequencies
                reflection_paths(path_counter).frequencies = analysis_frequencies; % Frequency vector
                reflection_paths(path_counter).offset = [2*nx*roomDim(1), 2*ny*roomDim(2), 2*nz*roomDim(3)];
                reflection_paths(path_counter).type = 'rectangular';
                
                % FIXED: Also store single coefficient for backwards compatibility (use 1kHz)
                [~, ref_freq_idx] = min(abs(analysis_frequencies - 1000));
                reflection_paths(path_counter).coefficient = total_reflection_all_freq(ref_freq_idx);
                
                path_counter = path_counter + 1;
            end
        end
    end
    
else  % Cylindrical tunnel
    % FIXED: Frequency-dependent cylindrical reflections
    num_segments = 16;
    
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

fprintf('Pre-computed %d reflection paths with %d frequency bands in %.2f seconds\n', ...
    length(reflection_paths), num_freq_bands, toc);

% VERIFICATION: Print first few reflection paths to verify frequency data
fprintf('\nVerification - First 3 reflection paths:\n');
for i = 1:min(3, length(reflection_paths))
    fprintf('Path %d: indices=[%d,%d,%d], type=%s\n', i, ...
        reflection_paths(i).indices(1), reflection_paths(i).indices(2), reflection_paths(i).indices(3), ...
        reflection_paths(i).type);
    fprintf('  Coefficients by frequency: %s\n', mat2str(reflection_paths(i).coefficients_freq, 3));
    fprintf('  Single coeff (1kHz): %.3f\n', reflection_paths(i).coefficient);
end

%% Vehicle Body Reflection Setup (if enabled)
if use_vehicle_body
    fprintf('Setting up vehicle body reflections...\n');
    
    % Define vehicle surfaces as planes with normals
    vehicle_surfaces = struct();
    
    % Top surface
    vehicle_surfaces(1).normal = [0, 0, 1];  % Points up
    vehicle_surfaces(1).offset = vehicle_height/2;
    vehicle_surfaces(1).reflection = 0.8;    % Reflection coefficient
    
    % Bottom surface
    vehicle_surfaces(2).normal = [0, 0, -1]; % Points down
    vehicle_surfaces(2).offset = vehicle_height/2;
    vehicle_surfaces(2).reflection = 0.7;
    
    % Front surface
    vehicle_surfaces(3).normal = [1, 0, 0];  % Points forward
    vehicle_surfaces(3).offset = vehicle_length/2;
    vehicle_surfaces(3).reflection = 0.8;
    
    % Back surface
    vehicle_surfaces(4).normal = [-1, 0, 0]; % Points backward
    vehicle_surfaces(4).offset = vehicle_length/2;
    vehicle_surfaces(4).reflection = 0.8;
    
    % Left side
    vehicle_surfaces(5).normal = [0, 1, 0];  % Points left
    vehicle_surfaces(5).offset = vehicle_width/2;
    vehicle_surfaces(5).reflection = 0.8;
    
    % Right side
    vehicle_surfaces(6).normal = [0, -1, 0]; % Points right
    vehicle_surfaces(6).offset = vehicle_width/2;
    vehicle_surfaces(6).reflection = 0.8;
    
    % We'll check for vehicle reflections during processing
end

%% NOTAP Source-Time Dominant Approach
fprintf('Initializing Source-Time Dominant processing...\n');

% Initialize Output Buffer for Each Microphone
output_signals = zeros(numMics, num_samples);

% For visualization purposes
direct_signals = zeros(numMics, num_samples);
reflection_signals = zeros(numMics, num_samples);
vehicle_reflection_signals = zeros(numMics, num_samples);
doppler_factors = zeros(numMics, num_samples);

% NEW: Data structures for localization
if save_mappings_for_localization
    % For direct path
    localization_data = struct();
    localization_data.direct_path = struct('mic_index', {}, 'source_times', {}, ...
        'receiver_times', {}, 'source_positions', {}, 'distances', {}, 'doppler_factors', {});
    
    % For reflections - UPDATED to include frequency-dependent data
    localization_data.reflections = struct('mic_index', {}, 'path_index', {}, 'indices', {}, ...
    'source_times', {}, 'receiver_times', {}, 'mirror_positions', {}, ...
    'distances', {}, 'reflection_coef', {}, 'doppler_factors', {}, ...
    'reflection_coef_freq', {}, 'frequencies', {});
end

%% Visualization: Room Setup
if enable_visualization
    fprintf('Setting up visualizations...\n');
    figure(1);
    hold on;
    
    % Plot room
    if strcmp(room_type, 'rectangular')
        % Plot rectangular room
        x = [0 0 0 0 roomDim(1) roomDim(1) roomDim(1) roomDim(1)];
        y = [0 0 roomDim(2) roomDim(2) 0 0 roomDim(2) roomDim(2)];
        z = [0 roomDim(3) 0 roomDim(3) 0 roomDim(3) 0 roomDim(3)];
        faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
        patch('Vertices', [x(:) y(:) z(:)], 'Faces', faces, ...
            'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    else  % Cylindrical tunnel
        % Plot cylindrical tunnel
        tunnel_center = [tunnel_length/2, tunnel_radius, tunnel_radius];
        
        % Draw the cylindrical part
        [X, Y, Z] = cylinder(tunnel_radius, 36);
        X = X * tunnel_length;
        Y = Y + tunnel_radius;
        Z = Z + tunnel_radius;
        surf(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineStyle', ':');
        
        % Draw the circular ends
        theta = linspace(0, 2*pi, 37);
        x_circle = tunnel_radius * cos(theta) + tunnel_radius;
        y_circle = tunnel_radius * sin(theta) + tunnel_radius;
        
        % Entrance circle
        plot3(zeros(size(theta)), x_circle, y_circle, 'k-', 'LineWidth', 1.5);
        
        % Exit circle
        plot3(repmat(tunnel_length, size(theta)), x_circle, y_circle, 'k-', 'LineWidth', 1.5);
        
        % Draw the floor
        floor_x = [0, tunnel_length, tunnel_length, 0];
        floor_y = [0, 0, 2*tunnel_radius, 2*tunnel_radius];
        floor_z = zeros(1, 4);
        patch('Vertices', [floor_x(:) floor_y(:) floor_z(:)], 'Faces', [1 2 3 4], ...
            'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    end
    
    % Plot source trajectory
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
    
    % Plot vehicle body if enabled
    if use_vehicle_body
        % Draw simple vehicle body at starting position
        vehicle_pos = source_trajectory(0);
        vehicle_corners = [
            vehicle_pos(1) - vehicle_length/2, vehicle_pos(2) - vehicle_width/2, vehicle_pos(3) - vehicle_height/2;  % Bottom-left-back
            vehicle_pos(1) - vehicle_length/2, vehicle_pos(2) + vehicle_width/2, vehicle_pos(3) - vehicle_height/2;  % Bottom-right-back
            vehicle_pos(1) + vehicle_length/2, vehicle_pos(2) + vehicle_width/2, vehicle_pos(3) - vehicle_height/2;  % Bottom-right-front
            vehicle_pos(1) + vehicle_length/2, vehicle_pos(2) - vehicle_width/2, vehicle_pos(3) - vehicle_height/2;  % Bottom-left-front
            vehicle_pos(1) - vehicle_length/2, vehicle_pos(2) - vehicle_width/2, vehicle_pos(3) + vehicle_height/2;  % Top-left-back
            vehicle_pos(1) - vehicle_length/2, vehicle_pos(2) + vehicle_width/2, vehicle_pos(3) + vehicle_height/2;  % Top-right-back
            vehicle_pos(1) + vehicle_length/2, vehicle_pos(2) + vehicle_width/2, vehicle_pos(3) + vehicle_height/2;  % Top-right-front
            vehicle_pos(1) + vehicle_length/2, vehicle_pos(2) - vehicle_width/2, vehicle_pos(3) + vehicle_height/2   % Top-left-front
        ];
        
        % Define faces based on corners
        vehicle_faces = [
            1 2 6 5;  % Back face
            3 4 8 7;  % Front face
            1 4 8 5;  % Left face
            2 3 7 6;  % Right face
            5 6 7 8;  % Top face
            1 2 3 4   % Bottom face
        ];
        
        % Plot the vehicle
        patch('Vertices', vehicle_corners, 'Faces', vehicle_faces, ...
            'FaceColor', [0.7 0.7 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    end
    
    % Plot microphone array
    scatter3(micPositions(:,1), micPositions(:,2), micPositions(:,3), 25, 'b', 'filled');
    text(arrayCenter(1), arrayCenter(2)-0.2, arrayCenter(3), 'Microphone Array', 'Color', 'b');
    
    % Set plot properties
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
disp('Starting Fast Hybrid Doppler-Reverb simulation with localization mapping...');
total_start_time = tic;

% For each microphone
for mic = 1:numMics
    mic_pos = micPositions(mic, :);
    fprintf('Processing microphone %d/%d...\n', mic, numMics);
    
    % SOURCE-TIME DOMINANT APPROACH
    % Step 1: Calculate transmission times from source to receiver
    % This is forward propagation, not requiring iterative solutions
    fprintf('  Calculating transmission times using source-time dominant approach...\n');
    start_time = tic;
    
    % For each source sample, calculate when it arrives at the receiver
    source_positions = zeros(length(t_s), 3);
    for i = 1:length(t_s)
        source_positions(i,:) = source_trajectory(t_s(i));
    end
    
    % Vectorized distance calculation
    distances = sqrt(sum((source_positions - repmat(mic_pos, length(t_s), 1)).^2, 2));
    
    % Calculate transmission times
    t_x_direct = t_s + distances / c;
    
    % Calculate distance changes for Doppler factor
    if is_moving
        % dr = diff(distances);
        dt_x = diff(t_x_direct);
        
        % Calculate instantaneous sampling frequency ratio (Doppler factor)
        f_s_x = 1 ./ dt_x;
        doppler_raw = fs ./ f_s_x;
        
        % Pad first value for consistency
        doppler_raw = [doppler_raw(1); doppler_raw];
        
        % Simple smoothing for Doppler factor
        window_size = 11;
        doppler_smooth = smoothdata(doppler_raw, 'movmean', window_size);
    else
        % Stationary source - no Doppler shift
        doppler_smooth = ones(length(t_s), 1);
    end
    
    % Store for output (this line stays the same)
    doppler_factors(mic, 1:length(doppler_smooth)) = doppler_smooth(1:min(length(doppler_smooth), size(doppler_factors, 2)));
    
    % NEW: Store mapping data for localization
    if save_mappings_for_localization
        direct_path_data = struct();
        direct_path_data.mic_index = mic;
        direct_path_data.source_times = t_s;
        direct_path_data.receiver_times = t_x_direct;
        direct_path_data.source_positions = source_positions;
        direct_path_data.distances = distances;
        direct_path_data.doppler_factors = doppler_smooth;
        
        % Add to localization data
        localization_data.direct_path(mic) = direct_path_data;
    end
    
    % Step 2: Apply NOTAP approach for direct path
    fprintf('  Processing direct path using NOTAP approach...\n');
    
    % Create anti-aliasing filter (needed for both moving and stationary)
    cutoff_freq = 0.8 * fs/2;  % Cutoff just below Nyquist
    [b_filter, a_filter] = butter(6, cutoff_freq/(fs_over/2));
    
    if is_moving
        % Create oversampled grid for interpolation (NOTAP key step)
        t_over = min(t_x_direct):1/fs_over:max(t_x_direct);
        
        % Zero-order hold interpolation (as described in NOTAP method)
        p_over = zeros(size(t_over));
        
        for i = 1:length(t_over)
            % Find the last source sample that would have reached the receiver by this time
            idx = find(t_x_direct <= t_over(i), 1, 'last');
            
            if ~isempty(idx) && idx <= length(source_signal)
                % Apply distance-based amplitude scaling
                distance = distances(idx);
                % attenuation = 1 / (1 + attenuation_beta * distance + attenuation_gamma * distance^2); % _this is not physically accurate 
                % Physically accurate: inverse square law with atmospheric absorption
                attenuation = 1 / (4 * pi * distance^2) * exp(-0.01 * distance);
                p_over(i) = source_signal(idx) * attenuation;
            end
        end
        
        % Apply anti-aliasing filter
        p_filtered = filtfilt(b_filter, a_filter, p_over);
        
        % Downsample to receiver rate
        p_direct = p_filtered(1:oversample_ratio:end);
        
    else
        % Stationary source - use simple impulse response approach
        distance = distances(1); % Constant distance for stationary source
        delay_samples = round(distance/c * fs); % use this t_x_direct instead and no need for this redundant part  
        
        % Initialize direct path signal
        p_direct = zeros(num_samples, 1);
        
        % Apply direct sound with proper delay and amplitude scaling
        if delay_samples < num_samples
            % attenuation = 1 / (1 + attenuation_beta * distance + attenuation_gamma * distance^2); %_not phyisically accurate
            % Physically accurate: inverse square law with atmospheric absorption
            attenuation = 1 / (4 * pi * distance^2) * exp(-0.01 * distance);
            
            % For stationary source, simply place the source signal at the correct delay
            signal_end = min(delay_samples + length(source_signal), num_samples);
            signal_length = signal_end - delay_samples;
            
            if signal_length > 0
                p_direct(delay_samples+1:signal_end) = source_signal(1:signal_length) * attenuation;
            end
        end
    end
    
    % Trim to match output length (this section stays the same)
    if length(p_direct) > num_samples
        p_direct = p_direct(1:num_samples);
    else
        % Pad with zeros if shorter
        p_direct(end+1:num_samples) = 0;
    end
    
    % Store direct path signal
    direct_signals(mic, :) = p_direct;
    
    fprintf('  Direct path processing completed in %.2f seconds\n', toc(start_time));
    % Generate WAV file for direct path only (first microphone)
    if mic == 1
        signal_normalized = 0.9 * p_direct / max(abs(p_direct));
        audiowrite('direct_path_mic1_only.wav', signal_normalized, fs);
        fprintf('  Direct path WAV saved for microphone 1\n');
    end
    
    % Step 3: Process reflections
    fprintf('  Processing reflections...\n');
    start_time = tic;

    %% FIXED: Reflection Processing with Proper Frequency-Dependent Coefficients
    % This replaces the reflection processing section in your frequency domain approach
    
    if use_freq_domain_reflections
        % FREQUENCY-DOMAIN APPROACH - FIXED for frequency-dependent reflections
        fprintf('  Using efficient frequency-domain approach with frequency-dependent reflections...\n');
        
        % Prepare frequency domain
        NFFT = 2^nextpow2(max(length(source_signal), num_samples));
        f = fs/2 * linspace(0, 1, NFFT/2+1);
        
        % FIXED: Process reflections with proper frequency-dependent coefficients
        reflection_signal_freq = zeros(NFFT, 1);
        
        for path_idx = 1:length(reflection_paths)
            current_path = reflection_paths(path_idx);
            
            if mod(path_idx, 5) == 0 || path_idx == length(reflection_paths)
                fprintf('    Processing reflection path %d/%d\n', path_idx, length(reflection_paths));
            end
            
            % Handle different room types
            if strcmp(room_type, 'rectangular') || strcmp(current_path.type, 'rectangular')
                mirror_offset = current_path.offset;
                mirror_trajectory = @(t) source_trajectory(t) + mirror_offset;
            else  % Cylindrical tunnel reflections
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
            
            % Calculate time-varying delay and attenuation
            num_points = 20;
            sample_times = linspace(0, min(duration, source_duration), num_points);
            
            % Precalculate mirror positions for this path
            mirror_positions = zeros(num_points, 3);
            for i = 1:num_points
                mirror_positions(i,:) = mirror_trajectory(sample_times(i));
            end
            
            % Calculate distances and delays
            mirror_distances = sqrt(sum((mirror_positions - repmat(mic_pos, num_points, 1)).^2, 2));
            delays = mirror_distances / c;
            
            % FIXED: Use frequency-dependent reflection coefficients
            % Get the reflection coefficients for all frequencies
            refl_coeffs_all_freq = current_path.coefficients_freq;
            refl_frequencies = current_path.frequencies;
            
            % Calculate Doppler factors
            if is_moving
                velocity_term = diff(mirror_distances) ./ diff(sample_times');
                doppler_factors_mirror = c ./ (c - velocity_term);
                doppler_factors_mirror = [doppler_factors_mirror(1); doppler_factors_mirror];
            else
                doppler_factors_mirror = ones(num_points, 1);
            end
            
            % Process each time segment with frequency-dependent coefficients
            segment_length = ceil(length(source_signal) / num_points);
            
            for seg = 1:num_points
                seg_start = (seg-1) * segment_length + 1;
                seg_end = min(seg * segment_length, length(source_signal));
                
                if seg_end >= seg_start
                    % Get source segment
                    source_segment = source_signal(seg_start:seg_end);
                    source_segment_padded = [source_segment; zeros(NFFT - length(source_segment), 1)];
                    
                    % Transform to frequency domain
                    source_freq = fft(source_segment_padded);
                    
                    % Apply Doppler shift
                    doppler_factor = doppler_factors_mirror(seg);
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
                    
                    % FIXED: Create frequency-dependent transfer function
                    delay = delays(seg);
                    
                    % Create full frequency vector for transfer function
                    f_full = fs * (0:NFFT-1)' / NFFT;
                    
                    % Interpolate reflection coefficients to match frequency grid
                    if length(refl_coeffs_all_freq) > 1
                        % Interpolate reflection coefficients to the full frequency grid
                        refl_coeff_interp = interp1(refl_frequencies, refl_coeffs_all_freq, f_full(1:NFFT/2+1), 'linear', 'extrap');
                        % Mirror for negative frequencies
                        refl_coeff_full = [refl_coeff_interp; flipud(refl_coeff_interp(2:end-1))];
                    else
                        % Single frequency case (fallback)
                        refl_coeff_full = repmat(refl_coeffs_all_freq, NFFT, 1);
                    end
                    
                    % Apply distance attenuation
                    distance = mirror_distances(seg);
                    distance_attenuation = 1 / (4 * pi * distance^2) * exp(-0.01 * distance);
                    
                    % Create frequency-dependent transfer function
                    H_segment = refl_coeff_full .* distance_attenuation .* exp(-1j * 2 * pi * f_full * delay);
                    
                    % Apply reflection transfer function
                    reflection_segment_freq = source_freq .* H_segment;
                    
                    % Accumulate with time position offset
                    time_position = ((seg-1) * segment_length) / length(source_signal);
                    phase_shift = exp(-1j * 2 * pi * (0:NFFT-1)' * time_position);
                    reflection_signal_freq = reflection_signal_freq + (reflection_segment_freq .* phase_shift);
                end
            end
            
            % Store reflection path data for localization (unchanged)
            if save_mappings_for_localization
                reflection_path_data = struct();
                reflection_path_data.mic_index = mic;
                reflection_path_data.path_index = path_idx;
                reflection_path_data.indices = current_path.indices;
                reflection_path_data.source_times = sample_times;
                reflection_path_data.receiver_times = sample_times + delays;
                reflection_path_data.mirror_positions = mirror_positions;
                reflection_path_data.distances = mirror_distances;
                reflection_path_data.reflection_coef = current_path.coefficient; % Use reference coefficient
                reflection_path_data.reflection_coef_freq = refl_coeffs_all_freq; % NEW: Store full frequency data
                reflection_path_data.frequencies = refl_frequencies; % NEW: Store frequency vector
                reflection_path_data.doppler_factors = doppler_factors_mirror;
                
                localization_data.reflections(end+1) = reflection_path_data;
            end
        end
        
        % Convert back to time domain
        reflection_signal = real(ifft(reflection_signal_freq));
        reflection_signal = reflection_signal(1:num_samples);
        reflection_signals(mic, :) = reflection_signal;
        
        else
            % TIME-DOMAIN APPROACH - WITH OPTIMIZED METHOD SELECTION - more
            % accurate
            if optimized_rev_time_method
                % OPTIMIZED: Use same NOTAP source-time approach as direct path
                fprintf('  Using optimized time-domain approach (NOTAP-style) with frequency-dependent reflections...\n');
                
                reflection_signal = zeros(num_samples, 1);
                
                for path_idx = 1:length(reflection_paths)
                    current_path = reflection_paths(path_idx);
                    
                    if mod(path_idx, 5) == 0 || path_idx == length(reflection_paths)
                        fprintf('    Processing reflection path %d/%d\n', path_idx, length(reflection_paths));
                    end
                    
                    % Handle different room types
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
                    
                    % NOTAP APPROACH: Calculate mirror source positions and transmission times
                    mirror_positions = zeros(length(t_s), 3);
                    for i = 1:length(t_s)
                        mirror_positions(i,:) = mirror_trajectory(t_s(i));
                    end
                    
                    % Vectorized distance calculation
                    mirror_distances = sqrt(sum((mirror_positions - repmat(mic_pos, length(t_s), 1)).^2, 2));
                    
                    % Calculate transmission times (source-time dominant)
                    t_x_mirror = t_s + mirror_distances / c;
                    
                    % Get reflection coefficient (use 1kHz reference for time-domain)
                    if length(current_path.coefficients_freq) > 1
                        refl_coeff = current_path.coefficient; % 1kHz reference
                    else
                        refl_coeff = current_path.coefficients_freq(1);
                    end
                    
                    % Calculate Doppler factors for this reflection path
                    if is_moving
                        mirror_dt_x = diff(t_x_mirror);
                        mirror_f_s_x = 1 ./ mirror_dt_x;
                        mirror_doppler_raw = fs ./ mirror_f_s_x;
                        mirror_doppler_raw = [mirror_doppler_raw(1); mirror_doppler_raw];
                        mirror_doppler_smooth = smoothdata(mirror_doppler_raw, 'movmean', window_size);
                    else
                        mirror_doppler_smooth = ones(length(t_s), 1);
                    end
                    
                    % NOTAP METHOD: Create oversampled grid and interpolate
                    if is_moving
                        % Create oversampled grid for interpolation
                        t_over = min(t_x_mirror):1/fs_over:max(t_x_mirror);
                        
                        % Zero-order hold interpolation (NOTAP approach)
                        p_over = zeros(size(t_over));
                        
                        for i = 1:length(t_over)
                            % Find the last source sample that would have reached the receiver by this time
                            idx = find(t_x_mirror <= t_over(i), 1, 'last');
                            
                            if ~isempty(idx) && idx <= length(source_signal)
                                % Apply distance-based amplitude scaling and reflection coefficient
                                distance = mirror_distances(idx);
                                attenuation = refl_coeff / (4 * pi * distance^2) * exp(-0.01 * distance);
                                p_over(i) = source_signal(idx) * attenuation;
                            end
                        end
                        
                        % Apply anti-aliasing filter
                        p_filtered = filtfilt(b_filter, a_filter, p_over);
                        
                        % Downsample to receiver rate
                        p_reflection = p_filtered(1:oversample_ratio:end);
                        
                    else
                        % Stationary source - use simple impulse response approach
                        distance = mirror_distances(1); % Constant distance
                        delay_samples = round(distance/c * fs);
                        
                        % Initialize reflection signal for this path
                        p_reflection = zeros(num_samples, 1);
                        
                        % Apply reflection with proper delay and amplitude scaling
                        if delay_samples < num_samples
                            attenuation = refl_coeff / (4 * pi * distance^2) * exp(-0.01 * distance);
                            
                            signal_end = min(delay_samples + length(source_signal), num_samples);
                            signal_length = signal_end - delay_samples;
                            
                            if signal_length > 0
                                p_reflection(delay_samples+1:signal_end) = source_signal(1:signal_length) * attenuation;
                            end
                        end
                    end
                    
                    % Trim to match output length
                    if length(p_reflection) > num_samples
                        p_reflection = p_reflection(1:num_samples);
                    else
                        p_reflection(end+1:num_samples) = 0;
                    end
                    
                    % Accumulate this reflection path
                    reflection_signal = reflection_signal + p_reflection;
                    
                    % Store reflection path data for localization
                    if save_mappings_for_localization
                        reflection_path_data = struct();
                        reflection_path_data.mic_index = mic;
                        reflection_path_data.path_index = path_idx;
                        reflection_path_data.indices = current_path.indices;
                        reflection_path_data.source_times = t_s;
                        reflection_path_data.receiver_times = t_x_mirror;
                        reflection_path_data.mirror_positions = mirror_positions;
                        reflection_path_data.distances = mirror_distances;
                        reflection_path_data.reflection_coef = refl_coeff;
                        reflection_path_data.reflection_coef_freq = current_path.coefficients_freq;
                        reflection_path_data.frequencies = current_path.frequencies;
                        reflection_path_data.doppler_factors = mirror_doppler_smooth;
                        
                        localization_data.reflections(end+1) = reflection_path_data;
                    end
                end
                
                % Store reflection signal
                reflection_signals(mic, :) = reflection_signal;
                
            else
                % CURRENT: Use the existing chunked time-domain approach  - More memory-efficient but slower
                fprintf('  Using current chunked time-domain approach with frequency-dependent reflections...\n');
                
                reflection_signal = zeros(num_samples, 1);
                
                for path_idx = 1:length(reflection_paths)
                    current_path = reflection_paths(path_idx);
                    
                    if mod(path_idx, 5) == 0 || path_idx == length(reflection_paths)
                        fprintf('    Processing reflection path %d/%d\n', path_idx, length(reflection_paths));
                    end
                    
                    % Handle different room types (same as before)
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
                    
                    % Calculate mirror positions and distances
                    mirror_positions = zeros(length(t_s), 3);
                    for i = 1:length(t_s)
                        mirror_positions(i,:) = mirror_trajectory(t_s(i));
                    end
                    
                    mirror_distances = sqrt(sum((mirror_positions - repmat(mic_pos, length(t_s), 1)).^2, 2));
                    t_x_mirror = t_s + mirror_distances / c;
                    
                    % Calculate Doppler factors
                    if is_moving
                        mirror_dt_x = diff(t_x_mirror);
                        mirror_f_s_x = 1 ./ mirror_dt_x;
                        mirror_doppler_raw = fs ./ mirror_f_s_x;
                        mirror_doppler_raw = [mirror_doppler_raw(1); mirror_doppler_raw];
                        mirror_doppler_smooth = smoothdata(mirror_doppler_raw, 'movmean', window_size);
                    else
                        mirror_doppler_smooth = ones(length(t_s), 1);
                    end
                    
                    % Get reflection coefficient
                    if length(current_path.coefficients_freq) > 1
                        refl_coeff = current_path.coefficient; % 1kHz reference
                    else
                        refl_coeff = current_path.coefficients_freq(1);
                    end
                    
                    % Process in chunks (existing approach)
                    chunk_size = 2000;
                    mic_time_vector = (0:num_samples-1)'/fs;
                    
                    for chunk_start_idx = 1:chunk_size:num_samples
                        chunk_end_idx = min(chunk_start_idx + chunk_size - 1, num_samples);
                        chunk_times = mic_time_vector(chunk_start_idx:chunk_end_idx);
                        
                        for i = 1:length(chunk_times)
                            current_time = chunk_times(i);
                            idx = find(t_x_mirror <= current_time, 1, 'last');
                            
                            if ~isempty(idx) && idx <= length(source_signal)
                                distance = mirror_distances(idx);
                                attenuation = refl_coeff / (4 * pi * distance^2) * exp(-0.01 * distance);
                                
                                output_idx = chunk_start_idx + i - 1;
                                if output_idx <= num_samples
                                    reflection_signal(output_idx) = reflection_signal(output_idx) + source_signal(idx) * attenuation;
                                end
                            end
                        end
                    end
                    
                    % Store reflection path data for localization
                    if save_mappings_for_localization
                        reflection_path_data = struct();
                        reflection_path_data.mic_index = mic;
                        reflection_path_data.path_index = path_idx;
                        reflection_path_data.indices = current_path.indices;
                        reflection_path_data.source_times = t_s;
                        reflection_path_data.receiver_times = t_x_mirror;
                        reflection_path_data.mirror_positions = mirror_positions;
                        reflection_path_data.distances = mirror_distances;
                        reflection_path_data.reflection_coef = refl_coeff;
                        reflection_path_data.reflection_coef_freq = current_path.coefficients_freq;
                        reflection_path_data.frequencies = current_path.frequencies;
                        reflection_path_data.doppler_factors = mirror_doppler_smooth;
                        
                        localization_data.reflections(end+1) = reflection_path_data;
                    end
                end
                
                % Apply filter and store
                reflection_signal = filtfilt(b_filter, a_filter, reflection_signal);
                reflection_signals(mic, :) = reflection_signal;
            end
        end
    
    fprintf('  Reflection processing completed in %.2f seconds\n', toc(start_time));
    
    % Step 4: Process vehicle body reflections (if enabled)
    if use_vehicle_body
        fprintf('  Processing vehicle body reflections...\n');
        start_time = tic;
        
        % Process vehicle reflections in time domain
        vehicle_signal = zeros(num_samples, 1);
        
        % Process each vehicle surface
        for surface_idx = 1:length(vehicle_surfaces)
            current_surface = vehicle_surfaces(surface_idx);
            
            % Create a function to calculate the mirror position for this surface
            mirror_trajectory = @(t) calculateVehicleMirrorPosition(source_trajectory(t), vehicle_trajectory(t), current_surface);
            
            % Calculate mirror positions, distances and transmission times
            vehicle_mirror_positions = zeros(length(t_s), 3);
            for i = 1:length(t_s)
                vehicle_mirror_positions(i,:) = mirror_trajectory(t_s(i));
            end
            
            % Vectorized distance calculation
            vehicle_mirror_distances = sqrt(sum((vehicle_mirror_positions - repmat(mic_pos, length(t_s), 1)).^2, 2));
            
            % Calculate transmission times
            t_x_vehicle_mirror = t_s + vehicle_mirror_distances / c;
            
            % Use the same chunking approach as for room reflections
            chunk_size = 2000;
            mic_time_vector = (0:num_samples-1)'/fs;
            
            % Process in time chunks
            for chunk_start_idx = 1:chunk_size:num_samples
                chunk_end_idx = min(chunk_start_idx + chunk_size - 1, num_samples);
                chunk_times = mic_time_vector(chunk_start_idx:chunk_end_idx);
                
                for i = 1:length(chunk_times)
                    current_time = chunk_times(i);
                    idx = find(t_x_vehicle_mirror <= current_time, 1, 'last');
                    
                    if ~isempty(idx) && idx <= length(source_signal)
                        % Apply distance and reflection attenuation
                        distance = vehicle_mirror_distances(idx);
                        attenuation = current_surface.reflection / (1 + attenuation_beta * distance + attenuation_gamma * distance^2);
                        
                        output_idx = chunk_start_idx + i - 1;
                        if output_idx <= num_samples
                            vehicle_signal(output_idx) = vehicle_signal(output_idx) + source_signal(idx) * attenuation;
                        end
                    end
                end
            end
        end
        
        % Apply filter to vehicle reflections
        vehicle_signal = filtfilt(b_filter, a_filter, vehicle_signal);
        
        % Store vehicle reflection signal
        vehicle_reflection_signals(mic, :) = vehicle_signal;
        
        fprintf('  Vehicle reflection processing completed in %.2f seconds\n', toc(start_time));
    end
    
    % Combine direct and reflection components
    if use_vehicle_body
        output_signals(mic, :) = direct_signals(mic, :) + reflection_signals(mic, :) + vehicle_reflection_signals(mic, :);
    else
        output_signals(mic, :) = direct_signals(mic, :) + reflection_signals(mic, :);
    end
    
    fprintf('  Completed microphone %d/%d\n', mic, numMics);
end

% Processing time
total_processing_time = toc(total_start_time);
fprintf('Processing completed in %.2f seconds\n', total_processing_time);


% CREATE IMPULSE RESPONSE OUTPUT DIRECTORY
ir_output_dir = 'impulse_response_data';
if ~exist(ir_output_dir, 'dir')
    mkdir(ir_output_dir);
end

% Create a filename prefix for IR files
ir_config_prefix = sprintf('%s_%s_%s_mics%d', room_type, array_type, noise_type, numMics);



%% Add Noise
if enable_noise
    SNR_dB = 10; % Desired signal-to-noise ratio in dB
    
    % Add noise to signals
    disp('Adding noise to signals...');
    
    % Create noisy version of signals
    noisy_signals = zeros(size(output_signals));
    
    % Generate noise based on selected type
    switch noise_type
        case 'white'
            fprintf('Using white noise model (default)...\n');
            for mic = 1:numMics
                % Calculate signal power
                signal_power = mean(output_signals(mic,:).^2);
                
                % Calculate required noise power for the target SNR
                noise_power = signal_power / (10^(SNR_dB/10));
                
                % Generate white noise with the correct power
                noise = sqrt(noise_power) * randn(1, size(output_signals, 2));
                
                % Add noise to the signal
                noisy_signals(mic, :) = output_signals(mic, :) + noise;
            end
            
        case 'pink'
            fprintf('Using pink noise model...\n');
            for mic = 1:numMics
                % Calculate signal power
                signal_power = mean(output_signals(mic,:).^2);
                
                % Calculate required noise power for the target SNR
                noise_power = signal_power / (10^(SNR_dB/10));
                
                % Generate white noise first
                white_noise = randn(1, size(output_signals, 2) + 1000); % Extra samples for filter warmup
                
                % Apply 1/f filter to create pink noise (simple approximation)
                b = [0.049922035, -0.095993537, 0.050612699, -0.004408786];
                a = [1, -2.494956002, 2.017265875, -0.522189400];
                pink_noise = filter(b, a, white_noise);
                
                % Remove extra samples, normalize and scale
                pink_noise = pink_noise(1001:end);
                pink_noise = pink_noise / std(pink_noise) * sqrt(noise_power);
                
                % Add noise to the signal
                noisy_signals(mic, :) = output_signals(mic, :) + pink_noise;
            end
            
        case 'ambient'
            fprintf('Using spatially correlated ambient noise model...\n');
            % Generate a base ambient noise field
            nfft = size(output_signals, 2) + 1000; % Extra samples for filter warmup
            base_ambient = randn(1, nfft);
            
            % Apply a filter to shape the spectrum to be more like ambient noise
            b = [0.1, 0.2, 0.3, 0.2, 0.1]; % Simple low-pass filter
            a = 1;
            base_ambient = filter(b, a, base_ambient);
            base_ambient = base_ambient(1001:end);
            
            for mic = 1:numMics
                % Calculate signal power
                signal_power = mean(output_signals(mic,:).^2);
                
                % Calculate required noise power for the target SNR
                noise_power = signal_power / (10^(SNR_dB/10));
                
                % Create a microphone-specific variation of the ambient noise
                % that's correlated with the base noise but has some unique components
                mic_specific = randn(1, size(output_signals, 2));
                
                % Mix base ambient (70%) with mic-specific noise (30%)
                mixed_noise = 0.7 * base_ambient + 0.3 * mic_specific;
                mixed_noise = mixed_noise / std(mixed_noise) * sqrt(noise_power);
                
                % Add noise to the signal
                noisy_signals(mic, :) = output_signals(mic, :) + mixed_noise;
            end
    end
    
    % For visualization and saving, we'll use the noisy signals
    final_signals = noisy_signals;
    fprintf('Noise added at SNR = %d dB\n', SNR_dB);

else
    % Skip noise - use clean signals
    fprintf('Noise disabled - using clean signals\n');
    final_signals = output_signals;
end

report_noise_info = struct();
report_noise_info.enabled = enable_noise;
report_noise_info.type = noise_type;
if enable_noise
    report_noise_info.snr_db = SNR_dB;
else
    report_noise_info.snr_db = -1;
end


%% Extract Impulse Response via Deconvolution
fprintf('Extracting impulse responses via deconvolution...\n');

% Only deconvolve if using measurement signal
% Applies fade-out optimized for sweep artifacts
% Assumes the source has good spectral coverage
% if strcmp(signal_choice, 'measurement_sweep')
% Deconvolve to get actual impulse responses
simulated_impulse_responses = zeros(size(final_signals));

for mic = 1:numMics
    ir_temp = deconvolveImpulseResponse(final_signals(mic, :), source_signal, fs, signal_choice, use_custom_config, reg);
    % Handle size mismatch - take only what fits
    ir_length = min(length(ir_temp), size(simulated_impulse_responses, 2));
    simulated_impulse_responses(mic, 1:ir_length) = ir_temp(1:ir_length);
    if ir_length < size(simulated_impulse_responses, 2)
        simulated_impulse_responses(mic, ir_length+1:end) = 0; % Pad with zeros
    end
end

fprintf('✓ Deconvolution completed for %d microphones\n', numMics);
% else
%     % For other signals, use output directly
%     simulated_impulse_responses = final_signals;
%     fprintf('✓ Using final signals directly (no deconvolution)\n');
% end

% Trim to 2 seconds if needed for comparison with measurements
max_ir_length = round(fs * 2); % 2 seconds
if size(simulated_impulse_responses, 2) > max_ir_length
    simulated_impulse_responses = simulated_impulse_responses(:, 1:max_ir_length);
end

% Save for comparison - using same name as before
ir_filename = fullfile(ir_output_dir, sprintf('%s_impulse_responses.mat', ir_config_prefix));
save(ir_filename, 'simulated_impulse_responses', 'fs', 'micPositions');
fprintf('Simulated impulse responses saved to: %s\n', ir_filename);




%% Visualization: Signal Results
if enable_visualization
    % Select a representative microphone 
    if numMics == 1
        center_mic = 1; % Only one microphone available
    else
        % If there's an odd number of mics, choose the middle one, otherwise choose one of the central ones
        if mod(numMics, 2) == 1
            center_mic = ceil(numMics/2);
        else
            % For even number of mics, choose one of the central ones
            center_mic = numMics/2;
        end
    end
    
    % Plot time-domain signals
    figure(2);
    subplot(3,1,1);
    plot((0:length(direct_signals(center_mic,:))-1)/fs, direct_signals(center_mic,:));
    title('Direct Path Component at Center Microphone');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    subplot(3,1,2);
    plot((0:length(reflection_signals(center_mic,:))-1)/fs, reflection_signals(center_mic,:));
    title('Reflection Components at Center Microphone');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    subplot(3,1,3);
    plot((0:length(final_signals(center_mic,:))-1)/fs, final_signals(center_mic,:));
    title('Combined Signal at Center Microphone');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    % Plot spectrogram for center microphone
    figure(3);
    window_size = 1024;
    overlap = round(window_size*0.75);
    spectrogram(final_signals(center_mic,:), window_size, overlap, [], fs, 'yaxis');
    title(sprintf('Spectrogram of Signal at Center Microphone %d', center_mic));
    colorbar;
    
    % Plot Doppler factor
    figure(4);
    plot((0:length(doppler_factors(center_mic,:))-1)/fs, doppler_factors(center_mic,:));
    title('Doppler Factor at Center Microphone');
    xlabel('Time (s)');
    ylabel('Doppler Factor');
    grid on;
    ylim([0.8 1.2]);  % Adjust as needed for better visualization
    
    % NEW: Compare Doppler factors across microphones (only if multiple mics)
    if numMics > 1
        figure(5);
        hold on;
        
        % Choose a few microphones to display (not too many to keep the plot readable)
        num_mics_to_display = min(5, numMics);
        mic_indices = round(linspace(1, numMics, num_mics_to_display));
        
        for i = 1:length(mic_indices)
            mic_idx = mic_indices(i);
            plot((0:length(doppler_factors(mic_idx,:))-1)/fs, doppler_factors(mic_idx,:), ...
                'DisplayName', sprintf('Mic %d', mic_idx));
        end
        
        title('Doppler Factor Comparison Across Microphones');
        xlabel('Time (s)');
        ylabel('Doppler Factor');
        grid on;
        legend('Location', 'best');
        ylim([0.8 1.2]);
        
        % NEW: Compare signal amplitudes across microphones
        figure(6);
        hold on;
        
        for i = 1:length(mic_indices)
            mic_idx = mic_indices(i);
            % Plot RMS amplitude over time using a moving window
            window_len = 1000; % samples
            signal = final_signals(mic_idx,:);
            rms_amplitude = zeros(1, length(signal));
            
            for j = 1:length(signal)
                start_idx = max(1, j - window_len/2);
                end_idx = min(length(signal), j + window_len/2);
                window_data = signal(start_idx:end_idx);
                rms_amplitude(j) = sqrt(mean(window_data.^2));
            end
            
            % Smooth the RMS curve
            rms_amplitude = smoothdata(rms_amplitude, 'movmean', 500);
            
            % Plot
            plot((0:length(rms_amplitude)-1)/fs, rms_amplitude, ...
                'DisplayName', sprintf('Mic %d', mic_idx));
        end
        
        title('Signal Amplitude Comparison Across Microphones');
        xlabel('Time (s)');
        ylabel('RMS Amplitude');
        grid on;
        legend('Location', 'best');
    else
        fprintf('Single microphone: Skipping multi-microphone comparison plots (Figures 5-6)\n');
    end
    
    % Visualize source-receiver relationship for localization
    if save_mappings_for_localization
        % Plot timing relationships
        figure(7);
        
        % Direct path timing relationship
        subplot(2,1,1);
        plot(localization_data.direct_path(center_mic).source_times, ...
             localization_data.direct_path(center_mic).receiver_times, 'b.-');
        title('Source-Receiver Time Mapping (Direct Path)');
        xlabel('Source Time (s)');
        ylabel('Receiver Time (s)');
        grid on;
        
        % Reflection path timing relationships
        subplot(2,1,2);
        hold on;
        
        % Find reflections for center mic
        mic_reflections = find([localization_data.reflections.mic_index] == center_mic);
        
        if ~isempty(mic_reflections)
            % Plot first 3 reflections at most (to avoid clutter)
            for i = 1:min(3, length(mic_reflections))
                refl_idx = mic_reflections(i);
                refl_data = localization_data.reflections(refl_idx);
            
                % Define color map
                color_map = 'rgb';
                color_char = color_map(mod(i-1, 3) + 1);
            
                % Construct line spec (e.g., '.-r')
                line_spec = ['.-', color_char];
            
                % Plot with correct syntax
                plot(refl_data.source_times, refl_data.receiver_times, ...
                    line_spec, ...
                    'DisplayName', sprintf('Reflection [%d,%d,%d]', ...
                    refl_data.indices(1), refl_data.indices(2), refl_data.indices(3)));
            end
        end
        
        title('Source-Receiver Time Mapping (Reflections)');
        xlabel('Source Time (s)');
        ylabel('Receiver Time (s)');
        grid on;
        legend('Location', 'northwest');
        
        % Plot distance vs. time
        figure(8);
        
        % Direct path distance
        subplot(2,1,1);
        plot(localization_data.direct_path(center_mic).receiver_times, ...
             localization_data.direct_path(center_mic).distances, 'b.-');
        title('Distance vs. Receiver Time (Direct Path)');
        xlabel('Receiver Time (s)');
        ylabel('Distance (m)');
        grid on;
        
        % Reflection path distances
        subplot(2,1,2);
        hold on;
        
        if ~isempty(mic_reflections)
            % Plot first 3 reflections at most
            for i = 1:min(3, length(mic_reflections))
                refl_idx = mic_reflections(i);
                refl_data = localization_data.reflections(refl_idx);
            
                % Define color map
                color_map = 'rgb';
                color_char = color_map(mod(i-1, 3) + 1);
            
                % Construct line spec (e.g., '.-r')
                line_spec = ['.-', color_char];
            
                % Plot with correct syntax
                plot(refl_data.receiver_times, refl_data.distances, ...
                    line_spec, ...
                    'DisplayName', sprintf('Reflection [%d,%d,%d]', ...
                    refl_data.indices(1), refl_data.indices(2), refl_data.indices(3)));
            end
        end
        
        title('Distance vs. Receiver Time (Reflections)');
        xlabel('Receiver Time (s)');
        ylabel('Distance (m)');
        grid on;
        legend('Location', 'northwest');
    end
    
    % Save figures to files
    for fig_num = 1:8
        if isgraphics(fig_num)
            fig_filename = sprintf('figure_%d_%s_%s.png', fig_num, room_type, array_type);
            saveas(figure(fig_num), fig_filename);
            fprintf('Saved figure to %s\n', fig_filename);
        end
    end
end

%% Save Results
% Create output directory if it doesn't exist
output_dir = 'output_data';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Create a filename prefix based on configuration
config_prefix = sprintf('%s_%s_%s_mics%d', room_type, array_type, noise_type, numMics);

% Save microphone signals to WAV file
for mic = 1:numMics
    % Normalize signal to range [-1, 1] for WAV output
    signal_normalized = 0.9 * final_signals(mic,:) / max(abs(final_signals(mic,:)));
    % Generate filename
    filename = fullfile(output_dir, sprintf('%s_mic_%d_signal.wav', config_prefix, mic));
    % Write WAV file
    audiowrite(filename, signal_normalized, fs);
end

% Save all variables to MAT file
mat_filename = fullfile(output_dir, sprintf('%s_simulation_data.mat', config_prefix));
save(mat_filename, 'final_signals', 'direct_signals', ...
     'reflection_signals', 'micPositions', 'doppler_factors', 'fs', 't_r', 'is_moving', ...
     'source_trajectory', 'roomDim', 'total_processing_time', 'room_type', 'array_type', ...
     'noise_type', 'F_parameter');

% Save localization data separately
if save_mappings_for_localization
    loc_filename = fullfile(output_dir, sprintf('%s_localization_data.mat', config_prefix));
    save(loc_filename, 'localization_data');
    fprintf('Localization mapping data saved for source tracking analysis.\n');
end

fprintf('Enhanced Fast Hybrid Simulation with Accurate Doppler complete.\n');

% Helper functions for cylindrical room reflections
function mirror_pos = computeCylindricalMirrorPos(source_pos, ~, tunnel_radius, tunnel_center)
    % For cylindrical wall reflections, compute the reflection point
    % First, project to 2D (y-z plane)
    source_yz = source_pos(2:3);
    center_yz = tunnel_center(2:3);
    
    % Vector from center to source
    vec_to_source = source_yz - center_yz;
    
    % Normalize to get direction
    dir_to_source = vec_to_source / norm(vec_to_source);
    
    % Find the point on the cylinder in this direction
    point_on_cylinder = center_yz + tunnel_radius * dir_to_source;
    
    % Reflect source across this point
    reflection_point_yz = 2 * point_on_cylinder - source_yz;
    
    % Combine with original x coordinate
    mirror_pos = [source_pos(1), reflection_point_yz(1), reflection_point_yz(2)];
end

function reflected = reflect_across_floor(pos)
    % Mirror a position across the z=0 plane
    reflected = [pos(1), pos(2), -pos(3)];
end

function reflected = reflect_across_entrance(pos)
    % Mirror a position across the x=0 plane
    reflected = [-pos(1), pos(2), pos(3)];
end

function reflected = reflect_across_exit(pos, tunnel_length)
    % Mirror a position across the x=tunnel_length plane
    reflected = [2*tunnel_length - pos(1), pos(2), pos(3)];
end

function mirror_pos = calculateVehicleMirrorPosition(source_pos, vehicle_pos, surface)
    % Calculate the mirror position of the source across a vehicle surface
    
    % Adjust the normal vector based on vehicle orientation
    % For simplicity, we assume the vehicle doesn't rotate
    normal = surface.normal;
    
    % Calculate the plane equation: normal·(x - p0) = 0
    % where p0 is a point on the plane
    p0 = vehicle_pos + normal * surface.offset;
    
    % Distance from source to plane
    dist = dot(normal, source_pos - p0);
    
    % Calculate mirror position by reflecting across the plane
    mirror_pos = source_pos - 2 * dist * normal;
end

function sweepSignal = generateExpSweep(fs, duration, fStart, fEnd)
   % Generate exponential sine sweep
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
   
   % if use_custom_config
   %      fprintf('\nRegularization options:\n');
   %      reg_choice = input('Use default regularization or advanced? (d/a) [default: d]: ', 's');
   % 
   %      if strcmpi(reg_choice, 'a') || strcmpi(reg_choice, 'advanced')
   %          reg = calculateOptimalReg(X,signal_choice);
   %      else
   %          reg = reg_value; % Use the reg value passed from default_values
   %      end
   % 
   %  else
   %      reg = reg_value; % Use default when not in custom config mode
   %  end
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

%% Helper Functions
function val = get_user_input(prompt, default)
    user_input = strtrim(input(sprintf('%s [default: %s]: ', prompt, default), 's'));
    if isempty(user_input)
        val = default;
    else
        val = user_input;
    end
end

function val = get_user_boolean(prompt, default)
    user_input = strtrim(lower(input(sprintf('%s (true/false) [default: %s]: ', prompt, string(default)), 's')));
    if isempty(user_input)
        val = default;
    else
        val = any(strcmp(user_input, {'true', '1', 'yes', 'y'}));
    end
end

function val = get_user_number(prompt, default)
    user_input = strtrim(input(sprintf('%s [default: %d]: ', prompt, default), 's'));
    if isempty(user_input)
        val = default;
    else
        num_val = str2double(user_input);
        if isnan(num_val)
            warning('Invalid input. Using default value: %d\n', default);
            val = default;
        else
            val = num_val;
        end
    end
end

function reg = calculateAdaptiveReg(X, signal_type)
    % Calculate signal energy statistics
    X_energy = abs(X).^2;
    max_energy = max(X_energy);
    mean_energy = mean(X_energy);
    min_energy = min(X_energy(X_energy > 0)); % Exclude zeros
    
    % Energy ratio indicates spectral uniformity
    energy_ratio = min_energy / max_energy;
    
    switch signal_type
        case 'measurement_sweep'
            % Sweeps have good spectral coverage
            reg = 0.001 * max_energy;
        case 'chirp'
            % Similar to sweep but might be less uniform
            reg = 0.002 * max_energy;
        case 'siren'
            % Poor spectral coverage, need stronger regularization
            if energy_ratio < 0.01  % Very uneven spectrum
                reg = 0.01 * max_energy;  % 10x stronger
            else
                reg = 0.005 * max_energy; % 5x stronger
            end
        otherwise
            % Conservative default for unknown signals
            reg = 0.01 * max_energy;
    end
end


function reg = calculateSpectralAdaptiveReg(X)
    % Analyze spectral characteristics
    X_mag = abs(X);
    max_mag = max(X_mag);
    
    % Count frequency bins with significant energy (>1% of max)
    significant_bins = sum(X_mag > 0.01 * max_mag);
    total_bins = length(X);
    coverage_ratio = significant_bins / total_bins;
    
    % Calculate spectral flatness (measure of whiteness)
    % Closer to 1 = more uniform (white-like)
    % Closer to 0 = more tonal (spike-like)
    geo_mean = exp(mean(log(X_mag(X_mag > 0))));
    arith_mean = mean(X_mag);
    spectral_flatness = geo_mean / arith_mean;
    
    % Adaptive regularization
    base_reg = 0.001 * max(X_mag.^2);
    
    % Adjust based on coverage
    coverage_factor = 1 / (coverage_ratio + 0.1); % More coverage = less reg
    
    % Adjust based on flatness  
    flatness_factor = 1 / (spectral_flatness + 0.1); % More flat = less reg
    
    reg = base_reg * coverage_factor * flatness_factor;
    
    % Clamp to reasonable bounds
    reg = max(0.0001 * max(X_mag.^2), min(reg, 0.1 * max(X_mag.^2)));
end

function reg = calculateOptimalReg(X,signal_choice)
    % Combined approach for optimal regularization
    
    % Method 1: Signal type based baseline
    reg1 = calculateAdaptiveReg(X, signal_choice);
    
    % Method 2: Spectral analysis based
    reg2 = calculateSpectralAdaptiveReg(X);
    
    % Combine methods (geometric mean for balance)
    reg = sqrt(reg1 * reg2);
    
    % Safety bounds relative to signal energy
    max_energy = max(abs(X).^2);
    reg = max(0.00001, min(reg, 0.1 * max_energy));
    
    % Provide feedback
    fprintf('  Signal type factor: %.6f\n', reg1/max_energy);
    fprintf('  Spectral factor: %.6f\n', reg2/max_energy);
    fprintf('  Final regularization: %.6f (%.3f%% of max energy)\n', ...
        reg, 100*reg/max_energy);
end

function visualizeBlackmanFilter(FIR_length, fs, Q)
    % Visualize the effect of Blackman windowing on filter design
    % 
    % Inputs:
    %   FIR_length - Length of FIR filter (e.g., 31)
    %   fs - Sampling frequency (e.g., 16000)
    %   Q - Quality factor (e.g., 0.95)
    
    if nargin < 3
        Q = 0.95;
    end
    if nargin < 2
        fs = 16000;
    end
    if nargin < 1
        FIR_length = 31;
    end
    
    % Calculate cutoff frequency
    cutoff_freq = Q * fs/2;
    
    % Design ideal low-pass filter (sinc function)
    K = (FIR_length - 1) / 2;
    n = -(FIR_length-1)/2 : (FIR_length-1)/2;  % Symmetric indices
    
    % Ideal sinc filter (infinite length truncated)
    wc = 2*pi*cutoff_freq/fs;  % Normalized cutoff frequency
    h_ideal = sin(wc * n) ./ (pi * n);
    h_ideal(K+1) = wc/pi;  % Handle n=0 case (sinc(0) = 1)
    
    % Apply different windows
    win_rect = ones(1, FIR_length);                    % Rectangular (no window)
    win_hann = hann(FIR_length)';                      % Hann window
    win_hamming = hamming(FIR_length)';                % Hamming window
    win_blackman = blackman(FIR_length)';              % Blackman window
    
    % Normalize windows
    win_rect = win_rect / sum(win_rect);
    win_hann = win_hann / sum(win_hann);
    win_hamming = win_hamming / sum(win_hamming);
    win_blackman = win_blackman / sum(win_blackman);
    
    % Apply windows to ideal filter
    h_rect = h_ideal .* win_rect;
    h_hann = h_ideal .* win_hann;
    h_hamming = h_ideal .* win_hamming;
    h_blackman = h_ideal .* win_blackman;
    
    % Calculate frequency responses
    NFFT = 2048;
    [H_rect, f] = freqz(h_rect, 1, NFFT, fs);
    [H_hann, ~] = freqz(h_hann, 1, NFFT, fs);
    [H_hamming, ~] = freqz(h_hamming, 1, NFFT, fs);
    [H_blackman, ~] = freqz(h_blackman, 1, NFFT, fs);
    
    % Create comprehensive visualization
    figure('Position', [100, 100, 1200, 800]);
    
    % 1. Time domain - Window functions
    subplot(2,3,1);
    plot(1:FIR_length, win_rect, 'k-', 'LineWidth', 1.5); hold on;
    plot(1:FIR_length, win_hann, 'b-', 'LineWidth', 1.5);
    plot(1:FIR_length, win_hamming, 'g-', 'LineWidth', 1.5);
    plot(1:FIR_length, win_blackman, 'r-', 'LineWidth', 2);
    xlabel('Sample Index'); ylabel('Amplitude');
    title('Window Functions (Time Domain)');
    legend('Rectangular', 'Hann', 'Hamming', 'Blackman', 'Location', 'best');
    grid on;
    
    % 2. Time domain - Filter impulse responses
    subplot(2,3,2);
    plot(1:FIR_length, h_rect, 'k-', 'LineWidth', 1.5); hold on;
    plot(1:FIR_length, h_hann, 'b-', 'LineWidth', 1.5);
    plot(1:FIR_length, h_hamming, 'g-', 'LineWidth', 1.5);
    plot(1:FIR_length, h_blackman, 'r-', 'LineWidth', 2);
    xlabel('Sample Index'); ylabel('Amplitude');
    title('Filter Impulse Responses');
    legend('Rectangular', 'Hann', 'Hamming', 'Blackman', 'Location', 'best');
    grid on;
    
    % 3. Frequency response - Magnitude (Linear)
    subplot(2,3,3);
    plot(f/1000, abs(H_rect), 'k-', 'LineWidth', 1.5); hold on;
    plot(f/1000, abs(H_hann), 'b-', 'LineWidth', 1.5);
    plot(f/1000, abs(H_hamming), 'g-', 'LineWidth', 1.5);
    plot(f/1000, abs(H_blackman), 'r-', 'LineWidth', 2);
    % Mark cutoff frequency
    xline(cutoff_freq/1000, '--', 'Cutoff', 'LineWidth', 1);
    xlabel('Frequency (kHz)'); ylabel('Magnitude');
    title('Frequency Response (Linear Scale)');
    legend('Rectangular', 'Hann', 'Hamming', 'Blackman', 'Location', 'best');
    grid on; xlim([0 fs/2000]);
    
    % 4. Frequency response - Magnitude (dB)
    subplot(2,3,4);
    plot(f/1000, 20*log10(abs(H_rect)), 'k-', 'LineWidth', 1.5); hold on;
    plot(f/1000, 20*log10(abs(H_hann)), 'b-', 'LineWidth', 1.5);
    plot(f/1000, 20*log10(abs(H_hamming)), 'g-', 'LineWidth', 1.5);
    plot(f/1000, 20*log10(abs(H_blackman)), 'r-', 'LineWidth', 2);
    % Mark cutoff frequency
    xline(cutoff_freq/1000, '--', 'Cutoff', 'LineWidth', 1);
    xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
    title('Frequency Response (dB Scale)');
    legend('Rectangular', 'Hann', 'Hamming', 'Blackman', 'Location', 'best');
    grid on; xlim([0 fs/2000]); ylim([-80 5]);
    
    % 5. Phase response
    subplot(2,3,5);
    plot(f/1000, unwrap(angle(H_rect)), 'k-', 'LineWidth', 1.5); hold on;
    plot(f/1000, unwrap(angle(H_hann)), 'b-', 'LineWidth', 1.5);
    plot(f/1000, unwrap(angle(H_hamming)), 'g-', 'LineWidth', 1.5);
    plot(f/1000, unwrap(angle(H_blackman)), 'r-', 'LineWidth', 2);
    xlabel('Frequency (kHz)'); ylabel('Phase (radians)');
    title('Phase Response');
    legend('Rectangular', 'Hann', 'Hamming', 'Blackman', 'Location', 'best');
    grid on; xlim([0 fs/2000]);
    
    % 6. Stopband detail (zoomed view)
    subplot(2,3,6);
    stopband_start = cutoff_freq * 1.2;  % Start viewing 20% above cutoff
    idx_start = find(f >= stopband_start, 1);
    plot(f(idx_start:end)/1000, 20*log10(abs(H_rect(idx_start:end))), 'k-', 'LineWidth', 1.5); hold on;
    plot(f(idx_start:end)/1000, 20*log10(abs(H_hann(idx_start:end))), 'b-', 'LineWidth', 1.5);
    plot(f(idx_start:end)/1000, 20*log10(abs(H_hamming(idx_start:end))), 'g-', 'LineWidth', 1.5);
    plot(f(idx_start:end)/1000, 20*log10(abs(H_blackman(idx_start:end))), 'r-', 'LineWidth', 2);
    xlabel('Frequency (kHz)'); ylabel('Magnitude (dB)');
    title('Stopband Detail (Side Lobe Suppression)');
    legend('Rectangular', 'Hann', 'Hamming', 'Blackman', 'Location', 'best');
    grid on; ylim([-80 -10]);
    
    % Print filter characteristics
    fprintf('\n=== FILTER CHARACTERISTICS ===\n');
    fprintf('Filter Length: %d samples\n', FIR_length);
    fprintf('Sampling Frequency: %d Hz\n', fs);
    fprintf('Cutoff Frequency: %.1f Hz (Q=%.2f)\n', cutoff_freq, Q);
    
    % Calculate and display stopband attenuation
    stopband_freq = cutoff_freq * 1.5;  % 50% above cutoff
    [~, idx] = min(abs(f - stopband_freq));
    
    fprintf('\nStopband Attenuation at %.1f Hz:\n', stopband_freq);
    fprintf('  Rectangular: %.1f dB\n', 20*log10(abs(H_rect(idx))));
    fprintf('  Hann:        %.1f dB\n', 20*log10(abs(H_hann(idx))));
    fprintf('  Hamming:     %.1f dB\n', 20*log10(abs(H_hamming(idx))));
    fprintf('  Blackman:    %.1f dB\n', 20*log10(abs(H_blackman(idx))));
    
    % Calculate transition bandwidth (3dB to -40dB)
    mag_blackman_db = 20*log10(abs(H_blackman));
    idx_3db = find(mag_blackman_db <= -3, 1);
    idx_40db = find(mag_blackman_db <= -40, 1);
    if ~isempty(idx_3db) && ~isempty(idx_40db)
        transition_bw = f(idx_40db) - f(idx_3db);
        fprintf('\nBlackman Transition Bandwidth (3dB to 40dB): %.1f Hz\n', transition_bw);
    end
    
    sgtitle(sprintf('Blackman vs Other Windows - FIR Filter Comparison (Length=%d)', FIR_length));
end
