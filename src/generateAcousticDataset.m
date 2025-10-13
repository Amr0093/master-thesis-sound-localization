function dataset_results = generateAcousticDataset(params,positions,mic_config)
%GENERATEACOUSTICDATASET - Create WAV dataset for DNN training (GUI-Ready)
% 
% This function generates WAV files for multiple source positions using the
% integrated FastHybridDopplerReverbSimulator_GUI_Ready for acoustic source
% localization training datasets.
%
% INPUT:
%   params - Complete parameter structure containing:
%            - dataset_positions: Nx3 array of [x,y,z] source positions
%            - dataset_mic_config: struct with microphone selection settings
%            - dataset_noise_config: struct with noise settings
%            - All other simulation parameters from GUI
%
% OUTPUT:
%   dataset_results - struct with file paths, statistics, and metadata

% fprintf('\n=======================================================\n');
% fprintf('ACOUSTIC DATASET GENERATOR - GUI Ready Version\n');
% fprintf('=======================================================\n');
% 
% %% Extract dataset-specific parameters from params
% if isfield(params, 'dataset_positions') && ~isempty(params.dataset_positions)
%     positions = params.dataset_positions;
%     fprintf('Using passed positions: %d positions\n', size(positions, 1));
% else
%     % Default grid positions (safe for 4x4x4 room)
%     [x_grid, y_grid, z_grid] = meshgrid([1.5, 2.5, 3.5], [1.5, 2.5], [1.5, 2.5]);
%     positions = [x_grid(:), y_grid(:), z_grid(:)];
%     fprintf('Using default position grid: %d positions\n', size(positions, 1));
% end
% 
% if isfield(params, 'dataset_mic_config') && ~isempty(params.dataset_mic_config)
%     mic_config = params.dataset_mic_config;
% else
%     mic_config = struct();
%     mic_config.selection_mode = 'range';
%     mic_config.range = [1, 4];
%     mic_config.custom_indices = [1, 2, 3, 4];
%     mic_config.num_random = 4;
% end
% 
% if isfield(params, 'dataset_noise_config') && ~isempty(params.dataset_noise_config)
%     noise_config = params.dataset_noise_config;
% else
%     noise_config = struct();
%     noise_config.enable_noise = true;
%     noise_config.noise_type = 'white';
% end
% 
% %% Setup Output Directory
% outputDir = 'acoustic_wavs';
% if ~exist(outputDir, 'dir')
%     mkdir(outputDir);
%     fprintf('Created output directory: %s\n', outputDir);
% end
% 
% runExternalSimulatorComplete(params);
% 
% % Step 2: Clear the flag for normal operation
% clear EXTERNAL_CALL_FLAG;

%% Initialize Processing Variables
%% Setup Output Directory (MISSING - needs to be added)
outputDir = 'acoustic_wavs';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% Initialize Processing Variables
num_positions = size(positions, 1);
original_positions = positions; % Store a copy to protect from modification
wav_file_length = 8001;

% Initialize metadata arrays for CSV
metadata = cell(1, 12);
metadata{1,1} = 'file_name';
metadata{1,2} = 'source_idx';
metadata{1,3} = 'source_x';
metadata{1,4} = 'source_y';
metadata{1,5} = 'source_z';
metadata{1,6} = 'mic_idx';
metadata{1,7} = 'mic_x';
metadata{1,8} = 'mic_y';
metadata{1,9} = 'mic_z';
metadata{1,10} = 'noise_enabled';
metadata{1,11} = 'noise_type';
metadata{1,12} = 'noise_snr_db';

row_idx = 2; % Start from row 2 (after headers)

%% Main Processing Loop
fprintf('\nStarting dataset generation for %d positions...\n', num_positions);
total_files_created = 0;
processing_errors = 0;

%% Main Processing Loop
for pos_idx = 1:num_positions
    sourcePos = original_positions(pos_idx, :);
    
    try
        %% Prepare Parameters for GUI Simulator (create copy to avoid modifying original)
        sim_params = params; % Create a copy of the original params
        
        % REMOVE dataset-specific fields that simulator doesn't understand
        if isfield(sim_params, 'dataset_positions')
            sim_params = rmfield(sim_params, 'dataset_positions');
        end
        if isfield(sim_params, 'dataset_mic_config')
            sim_params = rmfield(sim_params, 'dataset_mic_config');
        end
        if isfield(sim_params, 'dataset_noise_config')
            sim_params = rmfield(sim_params, 'dataset_noise_config');
        end
        
        % Only set the varying source position - everything else comes from passed params
        sim_params.start_pos = sourcePos;
        sim_params.end_pos = sourcePos; % Stationary source
        
        % Override specific settings for batch processing
        sim_params.enable_visualization = false; % Disable for batch processing
        
        %% Call GUI Simulator (now with clean params identical to direct GUI call)
        results = FastHybridDopplerReverbSimulator_GUI_Ready(sim_params);
        
        if ~results.success
            error('Simulator failed: %s', results.error_message);
        end
        
        % Extract results from the returned structure
        final_signals = results.data.signals.final;
        micPositions = results.data.array.positions;
        fs = results.data.simulation.sampling_rate;
        
        % Get noise information
        if isfield(results.data, 'noise')
            noise_info = results.data.noise;
        else
            noise_info = struct('enabled', false, 'type', 'unknown', 'snr_db', NaN);
        end
        
        fprintf('  Simulator completed successfully\n');
        
        %% Determine Microphone Selection
        if strcmp(mic_config.selection_mode, 'all')
            mic_indices = 1:size(micPositions, 1);
        elseif strcmp(mic_config.selection_mode, 'range')
            max_mic = size(micPositions, 1);
            mic_indices = max(1, mic_config.range(1)):min(max_mic, mic_config.range(2));
        elseif strcmp(mic_config.selection_mode, 'custom')
            mic_indices = mic_config.custom_indices;
            mic_indices = mic_indices(mic_indices <= size(micPositions, 1));
        elseif strcmp(mic_config.selection_mode, 'random')
            num_available_mics = size(micPositions, 1);
            num_to_select = min(mic_config.num_random, num_available_mics);
            mic_indices = randperm(num_available_mics, num_to_select);
        else
            mic_indices = 1:min(4, size(micPositions, 1)); % Default fallback
        end
        
        % Validate mic_indices
        if isempty(mic_indices)
            warning('No valid microphones selected for position %d, using first microphone', pos_idx);
            mic_indices = 1;
        end
        
        fprintf('  Using %d microphones: %s\n', length(mic_indices), mat2str(mic_indices));
        
        %% Save WAV Files
        position_files_created = 0;
        
        for m = 1:length(mic_indices)
            mic_idx = mic_indices(m);
            mic_pos = micPositions(mic_idx, :);
            
            % Create filename
            source_str = sprintf('x%.2f_y%.2f_z%.2f', sourcePos(1), sourcePos(2), sourcePos(3));
            source_str = strrep(source_str, '.', 'p'); % Replace decimals with 'p'
            
            if noise_info.enabled
                signal_filename = sprintf('src%02d_%s_mic%d_%s_snr%ddB.wav', ...
                    pos_idx, source_str, mic_idx, noise_info.type, round(noise_info.snr_db));
            else
                signal_filename = sprintf('src%02d_%s_mic%d_clean.wav', pos_idx, source_str, mic_idx);
            end
            
            signal_filepath = fullfile(outputDir, signal_filename);
            
            % Extract and process signal
            mic_signal = final_signals(mic_idx, :);
            
            if length(mic_signal) > 0 && max(abs(mic_signal)) > 0
                % Control WAV file length
                if length(mic_signal) > wav_file_length
                    mic_signal = mic_signal(1:wav_file_length);
                elseif length(mic_signal) < wav_file_length
                    mic_signal(end+1:wav_file_length) = 0;
                end
                
                % Normalize to prevent clipping
                mic_signal = 0.9 * mic_signal / max(abs(mic_signal));
                audiowrite(signal_filepath, mic_signal, fs);
                
                % Add to metadata
                metadata{row_idx,1} = signal_filename;
                metadata{row_idx,2} = pos_idx;
                metadata{row_idx,3} = sourcePos(1);
                metadata{row_idx,4} = sourcePos(2);
                metadata{row_idx,5} = sourcePos(3);
                metadata{row_idx,6} = mic_idx;
                metadata{row_idx,7} = mic_pos(1);
                metadata{row_idx,8} = mic_pos(2);
                metadata{row_idx,9} = mic_pos(3);
                metadata{row_idx,10} = noise_info.enabled;
                metadata{row_idx,11} = noise_info.type;
                metadata{row_idx,12} = noise_info.snr_db;
                row_idx = row_idx + 1;
                
                position_files_created = position_files_created + 1;
                total_files_created = total_files_created + 1;
                
            else
                warning('Empty or zero signal for mic %d at position %d, skipping', mic_idx, pos_idx);
            end
        end
        
        fprintf('  Created %d WAV files for position %d\n', position_files_created, pos_idx);
        
    catch ME
        fprintf('  ❌ Error processing position %d: %s\n', pos_idx, ME.message);
        fprintf('  Skipping this position and continuing...\n');
        processing_errors = processing_errors + 1;
        continue;
    end
    
    % Clear large variables to manage memory
    clear final_signals results;
end

%% Save Metadata Files
fprintf('\n=======================================================\n');
fprintf('Saving metadata and summary files...\n');

% Write metadata to CSV file
metadata_filename = fullfile(outputDir, 'sound_source_metadata.csv');
fid = fopen(metadata_filename, 'w');
if fid == -1
    error('Could not create metadata file');
end

% Write headers and data
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
    metadata{1,1}, metadata{1,2}, metadata{1,3}, metadata{1,4}, ...
    metadata{1,5}, metadata{1,6}, metadata{1,7}, metadata{1,8}, metadata{1,9}, ...
    metadata{1,10}, metadata{1,11}, metadata{1,12});

for i = 2:(row_idx-1)
    fprintf(fid, '%s,%d,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%s,%s,%.1f\n', ...
        metadata{i,1}, metadata{i,2}, metadata{i,3}, metadata{i,4}, ...
        metadata{i,5}, metadata{i,6}, metadata{i,7}, metadata{i,8}, metadata{i,9}, ...
        string(metadata{i,10}), metadata{i,11}, metadata{i,12});
end
fclose(fid);

% Create positional lookup table
pos_lookup_filename = fullfile(outputDir, 'position_lookup.csv');
fid = fopen(pos_lookup_filename, 'w');
if fid == -1
    error('Could not create position lookup file');
end

fprintf(fid, 'position_idx,source_x,source_y,source_z\n');
for i = 1:num_positions
    fprintf(fid, '%d,%.6f,%.6f,%.6f\n', i, original_positions(i,1), ...
        original_positions(i,2), original_positions(i,3));
end
fclose(fid);

% Save microphone and configuration info
% Extract noise config from params (since it's not passed as separate parameter)
if isfield(params, 'enable_noise') && isfield(params, 'noise_type')
    noise_config = struct();
    noise_config.enable_noise = params.enable_noise;
    noise_config.noise_type = params.noise_type;
else
    noise_config = struct('enable_noise', true, 'noise_type', 'white');
end

% Save microphone and configuration info
config_filename = fullfile(outputDir, 'dataset_config.mat');
dataset_config = struct();
dataset_config.mic_config = mic_config;
dataset_config.noise_config = noise_config;
dataset_config.positions = original_positions;
dataset_config.wav_file_length = wav_file_length;
dataset_config.generation_timestamp = datestr(now);
save(config_filename, 'dataset_config');

%% Prepare Results Structure
dataset_results = struct();
dataset_results.success = true;
dataset_results.output_directory = outputDir;
dataset_results.total_positions = num_positions;
dataset_results.successful_positions = num_positions - processing_errors;
dataset_results.processing_errors = processing_errors;
dataset_results.total_wav_files = total_files_created;
dataset_results.files = struct();
dataset_results.files.metadata_csv = metadata_filename;
dataset_results.files.position_lookup_csv = pos_lookup_filename;
dataset_results.files.config_mat = config_filename;
dataset_results.config = dataset_config;

%% Summary Report
fprintf('\n=======================================================\n');
fprintf('DATASET GENERATION COMPLETE!\n');
fprintf('=======================================================\n');
fprintf('✓ Processed %d source positions\n', num_positions);
fprintf('✓ Successfully generated %d positions\n', dataset_results.successful_positions);
if processing_errors > 0
    fprintf('⚠ Processing errors: %d positions\n', processing_errors);
end
fprintf('✓ Total WAV files created: %d\n', total_files_created);
fprintf('✓ Microphone selection mode: %s\n', mic_config.selection_mode);
fprintf('✓ Noise configuration: %s (enabled: %s)\n', noise_config.noise_type, string(noise_config.enable_noise));
fprintf('\n📁 Output files:\n');
fprintf('   • WAV files: %s/\n', outputDir);
fprintf('   • Metadata: %s\n', metadata_filename);
fprintf('   • Position lookup: %s\n', pos_lookup_filename);
fprintf('   • Configuration: %s\n', config_filename);
fprintf('=======================================================\n');

end
function runExternalSimulatorComplete(params)
    fprintf('Debug: === runExternalSimulatorComplete STARTED ===\n');
    fprintf('Debug: Received %d parameters\n', length(fieldnames(params)));
    
    try
        fprintf('Debug: About to call FastHybridDopplerReverbSimulator_GUI_Ready\n');
        
        % Call the simulator directly
        results = FastHybridDopplerReverbSimulator_GUI_Ready(params);
        
        fprintf('Debug: Simulator returned successfully\n');
        fprintf('Debug: results.success = %s\n', string(results.success));
        
        % Store essential results in base workspace for other modules
        assignin('base', 'simulation_results', results);
        assignin('base', 'final_signals', results.data.signals.final);
        assignin('base', 'micPositions', results.data.array.positions);
        assignin('base', 'fs', results.data.simulation.sampling_rate);
        assignin('base', 'sourcePos', results.data.source.start_position);
        
        % Additional variables for compatibility with other modules
        if isfield(results.data, 'room')
            assignin('base', 'roomDim', results.data.room.dimensions);
        end
        if isfield(results.data, 'noise')
            assignin('base', 'SNR_dB', results.data.noise.snr_db);
        end
        
        fprintf('Debug: Results assigned to base workspace\n');
        
        % Check if simulation completed successfully
        if results.success
            fprintf('✅ Simulation completed successfully!\n');
            fprintf('Execution time: %.2f seconds\n', results.execution_time);
            fprintf('Output directory: %s\n', results.directories.base_dir);
        else
            error('Simulation failed: %s', results.error_message);
        end
        
    catch ME
        fprintf('Debug: ERROR in runExternalSimulatorComplete: %s\n', ME.message);
        fprintf('Debug: Error identifier: %s\n', ME.identifier);
        fprintf('Debug: Error stack:\n');
        for i = 1:length(ME.stack)
            fprintf('  %s (line %d) in %s\n', ME.stack(i).name, ME.stack(i).line, ME.stack(i).file);
        end
        rethrow(ME);
    end
end