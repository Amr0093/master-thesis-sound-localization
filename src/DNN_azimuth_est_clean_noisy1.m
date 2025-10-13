%% Azimuth Estimation DNN Script - Using Distance Estimation Techniques
% Compatible with main workflow data format:
% - position_mic_signals: [total_samples x 8001] - each row is one mic signal
% - source_to_mic_azimuths: [1 x total_samples] - one azimuth per signal

%% Check if running in training/prediction mode or just defining architecture
if ~exist('localization_data', 'var')
    % Architecture definition mode - use default size for network creation
    fprintf('Defining azimuth network architecture (distance techniques)...\n');
    
    %% Define Network Architecture for Single Microphone Signals
    layers = [
        % Input Layer (single microphone signal: 1 x 8001 x 1)
        imageInputLayer([1 8001 1], 'Name', 'input')
        
        %% 1. Frequency Analysis Layer
        convolution2dLayer([1 256], 64, 'Name', 'freq_conv', 'Padding', 'same')
        batchNormalizationLayer('Name', 'freq_bn')
        reluLayer('Name', 'freq_relu')
        maxPooling2dLayer([1 2], 'Name', 'freq_pool', 'Stride', [1 2])
        
        %% 2. Azimuth Feature Extraction Layers - Enhanced with distance techniques
        % Layer 1: Capture phase and directional patterns (using distance filter sizes)
        convolution2dLayer([1 64], 128, 'Name', 'azi_conv1', 'Padding', 'same')
        batchNormalizationLayer('Name', 'azi_bn1')
        reluLayer('Name', 'azi_relu1')
        maxPooling2dLayer([1 4], 'Name', 'azi_pool1', 'Stride', [1 4])
        
        % Layer 2: Deeper directional features
        convolution2dLayer([1 32], 128, 'Name', 'azi_conv2', 'Padding', 'same')
        batchNormalizationLayer('Name', 'azi_bn2')
        reluLayer('Name', 'azi_relu2')
        maxPooling2dLayer([1 4], 'Name', 'azi_pool2', 'Stride', [1 4])
        
        % Layer 3: Fine directional features
        convolution2dLayer([1 16], 96, 'Name', 'azi_conv3', 'Padding', 'same')
        batchNormalizationLayer('Name', 'azi_bn3')
        reluLayer('Name', 'azi_relu3')
        maxPooling2dLayer([1 2], 'Name', 'azi_pool3', 'Stride', [1 2])
        
        %% 3. Fully Connected Layers
        flattenLayer('Name', 'flatten')
        
        % First dense layer
        fullyConnectedLayer(512, 'Name', 'fc1')
        batchNormalizationLayer('Name', 'fc1_bn')
        reluLayer('Name', 'fc1_relu')
        dropoutLayer(0.5, 'Name', 'fc1_dropout')
        
        % Second dense layer
        fullyConnectedLayer(256, 'Name', 'fc2')
        batchNormalizationLayer('Name', 'fc2_bn')
        reluLayer('Name', 'fc2_relu')
        dropoutLayer(0.3, 'Name', 'fc2_dropout')
        
        % Third dense layer
        fullyConnectedLayer(128, 'Name', 'fc3')
        reluLayer('Name', 'fc3_relu')
        dropoutLayer(0.2, 'Name', 'fc3_dropout')
        
        % Output layer - regression for azimuth angle prediction
        fullyConnectedLayer(1, 'Name', 'fc_out')
        regressionLayer('Name', 'output')
    ];
    
    %% Create Network
    trainedNetwork_azimuth = layerGraph(layers);
    
    fprintf('Azimuth estimation network architecture created (distance techniques).\n');
    fprintf('Network expects input size: [1 x 8001 x 1] (single microphone signals)\n');
    return; % Exit if just defining architecture
end

%% Training Mode - Data Available
fprintf('Starting azimuth estimation training with distance techniques...\n');

% Validate input data format
fprintf('Input data validation:\n');
fprintf('  localization_data size: [%d x %d]\n', size(localization_data, 1), size(localization_data, 2));
fprintf('  localization_labels.azimuth length: %d\n', length(localization_labels.azimuth));

% Ensure data consistency
num_samples = size(localization_data, 1);
if length(localization_labels.azimuth) ~= num_samples
    error('Mismatch: %d signals but %d azimuth labels', num_samples, length(localization_labels.azimuth));
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
    learning_rate = 5e-4; % Default
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

% Get azimuth labels and ensure column vector
Y_train = localization_labels.azimuth(:);

fprintf('Training data prepared: %d samples\n', num_samples);
fprintf('Azimuth range: [%.1f, %.1f] degrees\n', min(Y_train), max(Y_train));

%% Data Augmentation (controlled by GUI)
if use_data_augmentation
    fprintf('Applying data augmentation (enabled in GUI)...\n');
    
    % Conservative noise factor like distance estimation
    noise_factor = 0.005; % Smaller noise factor for better feature preservation
    X_augmented = X_train;
    Y_augmented = Y_train;
    
    for i = 1:num_samples
        % Add Gaussian noise
        noisy_signal = X_train(1, :, 1, i) + noise_factor * randn(1, 8001);
        
        % Append to dataset
        X_augmented = cat(4, X_augmented, reshape(noisy_signal, [1, 8001, 1, 1]));
        Y_augmented = [Y_augmented; Y_train(i)];
    end
    
    % Use augmented data
    X_train = X_augmented;
    Y_train = Y_augmented;
    num_samples = length(Y_train);
    
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
    YTrain = Y_train;
    XVal = [];
    YVal = [];
    use_validation = false;
    fprintf('Training with %d samples (no validation)\n', num_samples);
else
    % Normal validation split for larger datasets
    cv = cvpartition(num_samples, 'HoldOut', validation_split);
    XTrain = X_train(:,:,:,cv.training);
    YTrain = Y_train(cv.training);
    XVal = X_train(:,:,:,cv.test);
    YVal = Y_train(cv.test);
    use_validation = true;
    fprintf('Data split: %d training, %d validation samples\n', ...
        size(XTrain, 4), size(XVal, 4));
end

%% Update Network Architecture for Training
layers = [
    % Input Layer
    imageInputLayer([1 8001 1], 'Name', 'input')
    
    %% 1. Frequency Analysis Layer
    convolution2dLayer([1 256], 64, 'Name', 'freq_conv', 'Padding', 'same')
    batchNormalizationLayer('Name', 'freq_bn')
    reluLayer('Name', 'freq_relu')
    maxPooling2dLayer([1 2], 'Name', 'freq_pool', 'Stride', [1 2])
    
    %% 2. Azimuth Feature Extraction Layers - Enhanced with distance techniques
    convolution2dLayer([1 64], 128, 'Name', 'azi_conv1', 'Padding', 'same')
    batchNormalizationLayer('Name', 'azi_bn1')
    reluLayer('Name', 'azi_relu1')
    maxPooling2dLayer([1 4], 'Name', 'azi_pool1', 'Stride', [1 4])
    
    convolution2dLayer([1 32], 128, 'Name', 'azi_conv2', 'Padding', 'same')
    batchNormalizationLayer('Name', 'azi_bn2')
    reluLayer('Name', 'azi_relu2')
    maxPooling2dLayer([1 4], 'Name', 'azi_pool2', 'Stride', [1 4])
    
    convolution2dLayer([1 16], 96, 'Name', 'azi_conv3', 'Padding', 'same')
    batchNormalizationLayer('Name', 'azi_bn3')
    reluLayer('Name', 'azi_relu3')
    maxPooling2dLayer([1 2], 'Name', 'azi_pool3', 'Stride', [1 2])
    
    %% 3. Fully Connected Layers
    flattenLayer('Name', 'flatten')
    
    fullyConnectedLayer(512, 'Name', 'fc1')
    batchNormalizationLayer('Name', 'fc1_bn')
    reluLayer('Name', 'fc1_relu')
    dropoutLayer(0.5, 'Name', 'fc1_dropout')
    
    fullyConnectedLayer(256, 'Name', 'fc2')
    batchNormalizationLayer('Name', 'fc2_bn')
    reluLayer('Name', 'fc2_relu')
    dropoutLayer(0.3, 'Name', 'fc2_dropout')
    
    fullyConnectedLayer(128, 'Name', 'fc3')
    reluLayer('Name', 'fc3_relu')
    dropoutLayer(0.2, 'Name', 'fc3_dropout')
    
    fullyConnectedLayer(1, 'Name', 'fc_out')
    regressionLayer('Name', 'output')
];

trainedNetwork_azimuth = layerGraph(layers);

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
        'LearnRateDropPeriod', 30, ...                   % 🔒 HARDCODED
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
fprintf('Starting azimuth estimation network training (distance techniques)...\n');
tic;
[trainedNetwork_azimuth, info] = trainNetwork(XTrain, YTrain, trainedNetwork_azimuth, options);
training_time = toc;
fprintf('Training completed in %.2f seconds\n', training_time);

%% Evaluate on validation set (if available)
if use_validation
    fprintf('Evaluating azimuth estimation on validation set...\n');
    predictions = predict(trainedNetwork_azimuth, XVal);
    eval_data = YVal;
    eval_type = 'Validation';
else
    fprintf('No validation set - evaluating on training data...\n');
    predictions = predict(trainedNetwork_azimuth, XTrain);
    eval_data = YTrain;
    eval_type = 'Training';
end

% Calculate metrics
rmse = sqrt(mean((predictions - eval_data).^2));
mae = mean(abs(predictions - eval_data));
max_error = max(abs(predictions - eval_data));
r_squared = 1 - sum((eval_data - predictions).^2) / sum((eval_data - mean(eval_data)).^2);

fprintf('\n%s Results:\n', eval_type);
fprintf('  RMSE: %.3f degrees\n', rmse);     % 3 decimal places like distance
fprintf('  MAE: %.3f degrees\n', mae);       % 3 decimal places like distance
fprintf('  Max Error: %.3f degrees\n', max_error); % 3 decimal places like distance
fprintf('  R-squared: %.4f\n', r_squared);

%% Visualize Results
figure('Position', [100, 100, 1200, 400]);

% Subplot 1: Predictions vs Truth
subplot(1, 3, 1);
scatter(eval_data, predictions, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;
plot([min(eval_data) max(eval_data)], [min(eval_data) max(eval_data)], 'r--', 'LineWidth', 2);
xlabel('True Azimuth (degrees)');
ylabel('Predicted Azimuth (degrees)');
title(sprintf('%s Azimuth Prediction\n(RMSE: %.3f°, R²: %.3f)', eval_type, rmse, r_squared));
grid on;
axis equal;
legend('Predictions', 'Perfect', 'Location', 'best');

% Subplot 2: Error Distribution
subplot(1, 3, 2);
errors = predictions - eval_data;
histogram(errors, min(20, length(errors)), 'FaceColor', 'blue', 'EdgeColor', 'black');
xlabel('Prediction Error (degrees)');
ylabel('Frequency');
title(sprintf('Error Distribution\n(MAE: %.3f°)', mae));
grid on;
xline(0, 'r--', 'LineWidth', 2);

% Subplot 3: Training Progress
subplot(1, 3, 3);
plot(info.TrainingLoss, 'b-', 'LineWidth', 2);
hold on;
if isfield(info, 'ValidationLoss')
    plot(info.ValidationLoss, 'r-', 'LineWidth', 2);
    legend('Training', 'Validation', 'Location', 'best');
else
    legend('Training', 'Location', 'best');
end
xlabel('Iteration');
ylabel('Loss');
title('Training Progress');
grid on;

%% Save Trained Network and Results
results = struct();
results.rmse = rmse;
results.mae = mae;
results.max_error = max_error;
results.r_squared = r_squared;
results.training_time = training_time;
results.training_info = info;
results.validation_predictions = predictions;
results.validation_truth = eval_data;
results.dataset_type = eval_type;

save('azimuth_model_distance_techniques.mat', 'trainedNetwork_azimuth', 'results');
fprintf('\nTrained azimuth estimation model saved to: azimuth_model_distance_techniques.mat\n');

%% Display Final Summary
fprintf('\n=================================================\n');
fprintf('AZIMUTH ESTIMATION TRAINING COMPLETE (DISTANCE TECHNIQUES)\n');
fprintf('=================================================\n');
fprintf('Training samples: %d\n', size(XTrain, 4));
if use_validation
    fprintf('Validation samples: %d\n', size(XVal, 4));
else
    fprintf('Validation samples: 0 (too few total samples)\n');
end
fprintf('Training time: %.1f seconds\n', training_time);
fprintf('Final RMSE: %.3f degrees\n', rmse);
fprintf('Final MAE: %.3f degrees\n', mae);
fprintf('R-squared: %.4f\n', r_squared);
fprintf('Model saved to: azimuth_model_distance_techniques.mat\n');
fprintf('=================================================\n');