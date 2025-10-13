%% Combined 3D Localization DNN Script
% Compatible with main workflow data format:
% - position_mic_signals: [total_samples x 8001] - each row is one mic signal
% - Outputs: azimuth, distance, elevation simultaneously

%% Check if running in training/prediction mode or just defining architecture
if ~exist('localization_data', 'var')
    % Architecture definition mode - use default size for network creation
    fprintf('Defining combined 3D localization network architecture...\n');
    
    %% Define Network Architecture for Single Microphone Signals
    layers = [
        % Input Layer (single microphone signal: 1 x 8001 x 1)
        imageInputLayer([1 8001 1], 'Name', 'input')
        
        %% 1. Shared Frequency Analysis Layer
        convolution2dLayer([1 256], 128, 'Name', 'freq_conv', 'Padding', 'same')
        batchNormalizationLayer('Name', 'freq_bn')
        reluLayer('Name', 'freq_relu')
        maxPooling2dLayer([1 2], 'Name', 'freq_pool', 'Stride', [1 2])
        
        %% 2. Multi-scale Feature Extraction Layers
        % Layer 1: Large-scale temporal patterns (for distance/reverberation)
        convolution2dLayer([1 64], 128, 'Name', 'temp_conv1', 'Padding', 'same')
        batchNormalizationLayer('Name', 'temp_bn1')
        reluLayer('Name', 'temp_relu1')
        maxPooling2dLayer([1 4], 'Name', 'temp_pool1', 'Stride', [1 4])
        
        % Layer 2: Medium-scale patterns (for azimuth/phase differences)
        convolution2dLayer([1 32], 128, 'Name', 'mid_conv2', 'Padding', 'same')
        batchNormalizationLayer('Name', 'mid_bn2')
        reluLayer('Name', 'mid_relu2')
        maxPooling2dLayer([1 4], 'Name', 'mid_pool2', 'Stride', [1 4])
        
        % Layer 3: Fine-scale patterns (for elevation/spectral notches)
        convolution2dLayer([1 16], 96, 'Name', 'fine_conv3', 'Padding', 'same')
        batchNormalizationLayer('Name', 'fine_bn3')
        reluLayer('Name', 'fine_relu3')
        maxPooling2dLayer([1 2], 'Name', 'fine_pool3', 'Stride', [1 2])
        
        %% 3. Shared Fully Connected Layers
        flattenLayer('Name', 'flatten')
        
        % First shared dense layer
        fullyConnectedLayer(1024, 'Name', 'shared_fc1')
        batchNormalizationLayer('Name', 'shared_bn1')
        reluLayer('Name', 'shared_relu1')
        dropoutLayer(0.5, 'Name', 'shared_dropout1')
        
        % Second shared dense layer
        fullyConnectedLayer(512, 'Name', 'shared_fc2')
        batchNormalizationLayer('Name', 'shared_bn2')
        reluLayer('Name', 'shared_relu2')
        dropoutLayer(0.3, 'Name', 'shared_dropout2')
        
        % Third shared dense layer
        fullyConnectedLayer(256, 'Name', 'shared_fc3')
        reluLayer('Name', 'shared_relu3')
        dropoutLayer(0.2, 'Name', 'shared_dropout3')
        
        %% 4. Output layer - Multi-output regression for [azimuth, distance, elevation]
        fullyConnectedLayer(3, 'Name', 'fc_out')
        regressionLayer('Name', 'output')
    ];
    
    %% Create Network
    trainedNetwork_combined = layerGraph(layers);
    
    fprintf('Combined 3D localization network architecture created.\n');
    fprintf('Network expects input size: [1 x 8001 x 1] (single microphone signals)\n');
    fprintf('Network outputs: [azimuth, distance, elevation] simultaneously\n');
    return; % Exit if just defining architecture
end

%% Training Mode - Data Available
fprintf('Starting combined 3D localization training with actual data...\n');

% Validate input data format
fprintf('Input data validation:\n');
fprintf('  localization_data size: [%d x %d]\n', size(localization_data, 1), size(localization_data, 2));
fprintf('  localization_labels.azimuth length: %d\n', length(localization_labels.azimuth));
fprintf('  localization_labels.distance length: %d\n', length(localization_labels.distance));
fprintf('  localization_labels.elevation length: %d\n', length(localization_labels.elevation));

% Ensure data consistency
num_samples = size(localization_data, 1);
if length(localization_labels.azimuth) ~= num_samples || ...
   length(localization_labels.distance) ~= num_samples || ...
   length(localization_labels.elevation) ~= num_samples
    error('Mismatch in sample counts between signals and labels');
end

%% Get GUI Parameters (with fallbacks) - ONLY ONCE!
if exist('gui_epochs', 'var') && ~isempty(gui_epochs)
    epochs = gui_epochs;
else
    epochs = 100; % Default
end

if exist('gui_batch_size', 'var') && ~isempty(gui_batch_size)
    batch_size = gui_batch_size;
else
    batch_size = 32; % Default
end

if exist('gui_learning_rate', 'var') && ~isempty(gui_learning_rate)
    learning_rate = gui_learning_rate;
else
    learning_rate = 3e-4; % Default slightly lower for combined
end

if exist('gui_validation_split', 'var') && ~isempty(gui_validation_split)
    validation_split = gui_validation_split;
else
    validation_split = 0.2; % Default
end

if exist('gui_data_augmentation', 'var') && ~isempty(gui_data_augmentation)
    use_data_augmentation = gui_data_augmentation;
else
    use_data_augmentation = true; % Default
end

if exist('gui_gpu_training', 'var') && ~isempty(gui_gpu_training)
    gpu_training = gui_gpu_training;
else
    gpu_training = true; % Default
end

fprintf('Using GUI parameters: epochs=%d, batch_size=%d, lr=%.4f, validation_split=%.2f\n', ...
    epochs, batch_size, learning_rate, validation_split);

%% Prepare Training Data
fprintf('Preparing training data from %d microphone signals...\n', num_samples);

% Reshape data for network: [1 x 8001 x 1 x num_samples]
X_train = zeros(1, 8001, 1, num_samples);

for i = 1:num_samples
    % Each row of localization_data is one microphone signal [1 x 8001]
    mic_signal = localization_data(i, :);
    
    % Normalize signal to [-1, 1] range
    if max(abs(mic_signal)) > 0
        mic_signal = mic_signal / max(abs(mic_signal));
    end
    
    % Store in training format [1 x 8001 x 1]
    X_train(1, :, 1, i) = mic_signal;
end

% Combine labels into multi-output format [num_samples x 3]
% Format: [azimuth, distance, elevation] per row
Y_train = [localization_labels.azimuth(:), ...
           localization_labels.distance(:), ...
           localization_labels.elevation(:)];

fprintf('Training data prepared: %d samples\n', num_samples);
fprintf('Azimuth range: [%.1f, %.1f] degrees\n', min(Y_train(:,1)), max(Y_train(:,1)));
fprintf('Distance range: [%.2f, %.2f] meters\n', min(Y_train(:,2)), max(Y_train(:,2)));
fprintf('Elevation range: [%.1f, %.1f] degrees\n', min(Y_train(:,3)), max(Y_train(:,3)));

%% Normalize outputs for better training
fprintf('Normalizing outputs for stable training...\n');

% Store normalization parameters
azimuth_mean = mean(Y_train(:,1));
azimuth_std = std(Y_train(:,1));
distance_mean = mean(Y_train(:,2));
distance_std = std(Y_train(:,2));
elevation_mean = mean(Y_train(:,3));
elevation_std = std(Y_train(:,3));

% Normalize each output
Y_train_norm = Y_train;
Y_train_norm(:,1) = (Y_train(:,1) - azimuth_mean) / azimuth_std;
Y_train_norm(:,2) = (Y_train(:,2) - distance_mean) / distance_std;
Y_train_norm(:,3) = (Y_train(:,3) - elevation_mean) / elevation_std;

%% Data Augmentation (controlled by GUI)
if use_data_augmentation
    fprintf('Applying data augmentation (enabled in GUI)...\n');
    
    % Add some noise for robustness
    noise_factor = 0.01;
    X_augmented = X_train;
    Y_augmented = Y_train_norm;
    
    for i = 1:num_samples
        % Add Gaussian noise
        noisy_signal = X_train(1, :, 1, i) + noise_factor * randn(1, 8001);
        
        % Append to dataset
        X_augmented = cat(4, X_augmented, reshape(noisy_signal, [1, 8001, 1, 1]));
        Y_augmented = [Y_augmented; Y_train_norm(i, :)];
    end
    
    % Use augmented data
    X_train = X_augmented;
    Y_train_norm = Y_augmented;
    num_samples = size(Y_train_norm, 1);
    
    fprintf('After augmentation: %d samples\n', num_samples);
    
else
    fprintf('Skipping data augmentation (disabled in GUI)...\n');
    % Keep original data unchanged
    fprintf('Using original dataset: %d samples\n', num_samples);
end

%% Split into training and validation sets (handle both small and large datasets)
fprintf('Splitting data into training and validation sets...\n');

if num_samples < 10
    % Too few samples for validation split - use all for training
    fprintf('Warning: Only %d samples available. Using all for training (no validation).\n', num_samples);
    XTrain = X_train;
    YTrain = Y_train_norm;
    XVal = [];
    YVal = [];
    use_validation = false;
    fprintf('Training with %d samples (no validation)\n', num_samples);
else
    % Normal validation split for larger datasets
    cv = cvpartition(num_samples, 'HoldOut', validation_split);
    XTrain = X_train(:,:,:,cv.training);
    YTrain = Y_train_norm(cv.training, :);
    XVal = X_train(:,:,:,cv.test);
    YVal = Y_train_norm(cv.test, :);
    use_validation = true;
    fprintf('Data split: %d training, %d validation samples\n', ...
        size(XTrain, 4), size(YVal, 1));
end

%% Update Network Architecture for Training
layers = [
    % Input Layer
    imageInputLayer([1 8001 1], 'Name', 'input')
    
    %% 1. Shared Frequency Analysis Layer
    convolution2dLayer([1 256], 128, 'Name', 'freq_conv', 'Padding', 'same')
    batchNormalizationLayer('Name', 'freq_bn')
    reluLayer('Name', 'freq_relu')
    maxPooling2dLayer([1 2], 'Name', 'freq_pool', 'Stride', [1 2])
    
    %% 2. Multi-scale Feature Extraction Layers
    convolution2dLayer([1 64], 128, 'Name', 'temp_conv1', 'Padding', 'same')
    batchNormalizationLayer('Name', 'temp_bn1')
    reluLayer('Name', 'temp_relu1')
    maxPooling2dLayer([1 4], 'Name', 'temp_pool1', 'Stride', [1 4])
    
    convolution2dLayer([1 32], 128, 'Name', 'mid_conv2', 'Padding', 'same')
    batchNormalizationLayer('Name', 'mid_bn2')
    reluLayer('Name', 'mid_relu2')
    maxPooling2dLayer([1 4], 'Name', 'mid_pool2', 'Stride', [1 4])
    
    convolution2dLayer([1 16], 96, 'Name', 'fine_conv3', 'Padding', 'same')
    batchNormalizationLayer('Name', 'fine_bn3')
    reluLayer('Name', 'fine_relu3')
    maxPooling2dLayer([1 2], 'Name', 'fine_pool3', 'Stride', [1 2])
    
    %% 3. Shared Fully Connected Layers
    flattenLayer('Name', 'flatten')
    
    fullyConnectedLayer(1024, 'Name', 'shared_fc1')
    batchNormalizationLayer('Name', 'shared_bn1')
    reluLayer('Name', 'shared_relu1')
    dropoutLayer(0.5, 'Name', 'shared_dropout1')
    
    fullyConnectedLayer(512, 'Name', 'shared_fc2')
    batchNormalizationLayer('Name', 'shared_bn2')
    reluLayer('Name', 'shared_relu2')
    dropoutLayer(0.3, 'Name', 'shared_dropout2')
    
    fullyConnectedLayer(256, 'Name', 'shared_fc3')
    reluLayer('Name', 'shared_relu3')
    dropoutLayer(0.2, 'Name', 'shared_dropout3')
    
    %% 4. Multi-output regression
    fullyConnectedLayer(3, 'Name', 'fc_out')
    regressionLayer('Name', 'output')
];

trainedNetwork_combined = layerGraph(layers);

%% Determine execution environment
if gpu_training
    execEnv = 'auto';
else
    execEnv = 'cpu';
end

%% Define Training Options - SINGLE DEFINITION ONLY
if use_validation
    options = trainingOptions('adam', ...
        'InitialLearnRate', learning_rate, ...           % 🎛️ GUI
        'MaxEpochs', epochs, ...                         % 🎛️ GUI  
        'MiniBatchSize', batch_size, ...                 % 🎛️ GUI
        'Shuffle', 'every-epoch', ...                    % 🔒 HARDCODED
        'ValidationData', {XVal, YVal}, ...              % 🔒 HARDCODED
        'ValidationFrequency', 30, ...                   % 🔒 HARDCODED
        'ValidationPatience', 15, ...                    % 🔒 HARDCODED
        'LearnRateSchedule', 'piecewise', ...            % 🔒 HARDCODED
        'LearnRateDropFactor', 0.5, ...                  % 🔒 HARDCODED
        'LearnRateDropPeriod', 25, ...                   % 🔒 HARDCODED
        'Verbose', true, ...                             % 🔒 HARDCODED
        'Plots', 'training-progress', ...                % 🔒 HARDCODED
        'ExecutionEnvironment', execEnv);                % 🎛️ GUI
else
    options = trainingOptions('adam', ...
        'InitialLearnRate', learning_rate * 0.5, ...    % 🎛️ GUI (conservative)
        'MaxEpochs', round(epochs * 0.8), ...           % 🎛️ GUI (reduced)
        'MiniBatchSize', min(batch_size/2, size(XTrain, 4)), ... % 🎛️ GUI (adaptive)
        'Shuffle', 'every-epoch', ...                    % 🔒 HARDCODED
        'Verbose', true, ...                             % 🔒 HARDCODED
        'ExecutionEnvironment', execEnv);                % 🎛️ GUI
end

%% Train Network
fprintf('Starting combined 3D localization network training...\n');
tic;
[trainedNetwork_combined, info] = trainNetwork(XTrain, YTrain, trainedNetwork_combined, options);
training_time = toc;
fprintf('Training completed in %.2f seconds\n', training_time);

%% Evaluate on validation set (if available)
if use_validation
    fprintf('Evaluating combined 3D localization on validation set...\n');
    predictions_norm = predict(trainedNetwork_combined, XVal);
    eval_data_norm = YVal;
    eval_type = 'Validation';
else
    fprintf('No validation set - evaluating on training data...\n');
    predictions_norm = predict(trainedNetwork_combined, XTrain);
    eval_data_norm = YTrain;
    eval_type = 'Training';
end

%% Denormalize predictions and ground truth for evaluation
predictions = predictions_norm;
predictions(:,1) = predictions_norm(:,1) * azimuth_std + azimuth_mean;
predictions(:,2) = predictions_norm(:,2) * distance_std + distance_mean;
predictions(:,3) = predictions_norm(:,3) * elevation_std + elevation_mean;

eval_data = eval_data_norm;
eval_data(:,1) = eval_data_norm(:,1) * azimuth_std + azimuth_mean;
eval_data(:,2) = eval_data_norm(:,2) * distance_std + distance_mean;
eval_data(:,3) = eval_data_norm(:,3) * elevation_std + elevation_mean;

%% Calculate metrics for each output
% Azimuth metrics
az_rmse = sqrt(mean((predictions(:,1) - eval_data(:,1)).^2));
az_mae = mean(abs(predictions(:,1) - eval_data(:,1)));

% Distance metrics
dist_rmse = sqrt(mean((predictions(:,2) - eval_data(:,2)).^2));
dist_mae = mean(abs(predictions(:,2) - eval_data(:,2)));

% Elevation metrics
elev_rmse = sqrt(mean((predictions(:,3) - eval_data(:,3)).^2));
elev_mae = mean(abs(predictions(:,3) - eval_data(:,3)));

fprintf('\n%s Results:\n', eval_type);
fprintf('  Azimuth   - RMSE: %.2f°, MAE: %.2f°\n', az_rmse, az_mae);
fprintf('  Distance  - RMSE: %.3f m, MAE: %.3f m\n', dist_rmse, dist_mae);
fprintf('  Elevation - RMSE: %.2f°, MAE: %.2f°\n', elev_rmse, elev_mae);

%% Visualize Results
figure('Position', [100, 100, 1500, 400]);

% Subplot 1: Azimuth
subplot(1, 3, 1);
scatter(eval_data(:,1), predictions(:,1), 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot([min(eval_data(:,1)) max(eval_data(:,1))], [min(eval_data(:,1)) max(eval_data(:,1))], 'r--', 'LineWidth', 2);
xlabel('True Azimuth (degrees)');
ylabel('Predicted Azimuth (degrees)');
title(sprintf('Azimuth\n(RMSE: %.1f°)', az_rmse));
grid on;
axis equal;

% Subplot 2: Distance
subplot(1, 3, 2);
scatter(eval_data(:,2), predictions(:,2), 50, 'g', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot([min(eval_data(:,2)) max(eval_data(:,2))], [min(eval_data(:,2)) max(eval_data(:,2))], 'r--', 'LineWidth', 2);
xlabel('True Distance (meters)');
ylabel('Predicted Distance (meters)');
title(sprintf('Distance\n(RMSE: %.3f m)', dist_rmse));
grid on;
axis equal;

% Subplot 3: Elevation
subplot(1, 3, 3);
scatter(eval_data(:,3), predictions(:,3), 50, 'm', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot([min(eval_data(:,3)) max(eval_data(:,3))], [min(eval_data(:,3)) max(eval_data(:,3))], 'r--', 'LineWidth', 2);
xlabel('True Elevation (degrees)');
ylabel('Predicted Elevation (degrees)');
title(sprintf('Elevation\n(RMSE: %.1f°)', elev_rmse));
grid on;
axis equal;

%% Save Trained Network and Results
results = struct();
results.azimuth_rmse = az_rmse;
results.azimuth_mae = az_mae;
results.distance_rmse = dist_rmse;
results.distance_mae = dist_mae;
results.elevation_rmse = elev_rmse;
results.elevation_mae = elev_mae;
results.training_time = training_time;
results.training_info = info;
results.validation_predictions = predictions;
results.validation_truth = eval_data;
results.dataset_type = eval_type;

% Save normalization parameters for future use
results.normalization = struct();
results.normalization.azimuth_mean = azimuth_mean;
results.normalization.azimuth_std = azimuth_std;
results.normalization.distance_mean = distance_mean;
results.normalization.distance_std = distance_std;
results.normalization.elevation_mean = elevation_mean;
results.normalization.elevation_std = elevation_std;

save('combined_3D_model.mat', 'trainedNetwork_combined', 'results');
fprintf('\nTrained combined 3D localization model saved to: combined_3D_model.mat\n');

%% Display Final Summary
fprintf('\n=======================================================\n');
fprintf('COMBINED 3D LOCALIZATION TRAINING COMPLETE\n');
fprintf('=======================================================\n');
fprintf('Training samples: %d\n', size(XTrain, 4));
if use_validation
    fprintf('Validation samples: %d\n', size(XVal, 4));
else
    fprintf('Validation samples: 0 (too few total samples)\n');
end
fprintf('Training time: %.1f seconds\n', training_time);
fprintf('Final Results:\n');
fprintf('  Azimuth   - RMSE: %.2f°, MAE: %.2f°\n', az_rmse, az_mae);
fprintf('  Distance  - RMSE: %.3f m, MAE: %.3f m\n', dist_rmse, dist_mae);
fprintf('  Elevation - RMSE: %.2f°, MAE: %.2f°\n', elev_rmse, elev_mae);
fprintf('Model saved to: combined_3D_model.mat\n');
fprintf('=======================================================\n');