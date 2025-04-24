% simulate_microphone_signals.m
% Step 3: Simulate received signals at each microphone (with and without noise)
clc;
clear;

%% Load precomputed data
load('room_impulse_response.mat'); % impulse_responses, fs, t, etc.
load('room_geometry_setup.mat'); % micPositions, numMics, etc.

% Display debug info about loaded data
disp('Loaded data information:');
whos impulse_responses
fprintf('Number of microphones: %d\n', numMics);
fprintf('Sampling frequency: %d Hz\n', fs);

%% Generate source signal (10 ms chirp from 500 Hz to 2000 Hz)
Tsig = 0.01; % Duration of signal in seconds
t_sig = 0:1/fs:Tsig;
source_signal = chirp(t_sig, 500, Tsig, 2000);

% Normalize source signal
source_signal = source_signal / max(abs(source_signal));

fprintf('Source signal length: %d samples\n', length(source_signal));

%% Convolve signal with impulse responses
signal_length = length(source_signal);
ir_length = length(impulse_responses(1,:)); % Should be same as n_samples from IR script
output_length = signal_length + ir_length - 1;

% Preallocate for full convolution
full_convolution = zeros(numMics, output_length);

% Perform convolution for each microphone
disp('Convolving source signal with impulse responses...');
for mic = 1:numMics
    % Get impulse response for this microphone
    ir = impulse_responses(mic, :);
    
    % Perform convolution
    full_convolution(mic, :) = conv(source_signal, ir);
    
    % Show progress for every 10th microphone
    if mod(mic, 10) == 0
        fprintf('Processed microphone %d of %d\n', mic, numMics);
    end
end

% Crop to a reasonable length for visualization and further processing
% We'll keep the same length as the original time vector t from the IR calculation
if length(t) < size(full_convolution, 2)
    clean_signals = full_convolution(:, 1:length(t));
else
    % Pad with zeros if needed
    clean_signals = [full_convolution zeros(numMics, length(t) - size(full_convolution, 2))];
end

fprintf('Clean signals size: [%d, %d]\n', size(clean_signals, 1), size(clean_signals, 2));

%% Add noise (simulate low SNR condition)
% The paper used SNR = 10 dB (relative amplitude of noise = 0.1)
SNR_dB = 10; % Set to 10 dB as specified in the paper
noisy_signals = zeros(size(clean_signals));

disp(['Adding Gaussian noise at ', num2str(SNR_dB), ' dB SNR...']);

for mic = 1:numMics
    % Calculate signal power
    signal_power = mean(clean_signals(mic,:).^2);
    
    % Calculate required noise power for the target SNR
    noise_power = signal_power / (10^(SNR_dB/10));
    
    % Generate noise with the correct power
    noise = sqrt(noise_power) * randn(1, size(clean_signals, 2));
    
    % Add noise to the signal
    noisy_signals(mic, :) = clean_signals(mic, :) + noise;
end

%% Plot results for center microphone
center_mic = ceil(numMics / 2);

figure;
subplot(2,1,1);
plot(t, clean_signals(center_mic,:));
title(sprintf('Clean Received Signal at Center Mic (Mic #%d)', center_mic));
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t, noisy_signals(center_mic,:));
title(sprintf('Noisy Signal at Center Mic (SNR = %d dB)', SNR_dB));
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Visualize sound pressure distribution (as in Fig. 3 of the paper)
% Extract a time slice at t = 0.1s (or closest point)
[~, time_idx] = min(abs(t - 0.1));
fprintf('Visualizing sound pressure at t = %.3f seconds\n', t(time_idx));

% Reshape the microphone signals to a 2D array for plotting
array_elements = sqrt(numMics); % Should be 11 for 11x11 array
sound_pressure_clean = reshape(clean_signals(:, time_idx), array_elements, array_elements);
sound_pressure_noisy = reshape(noisy_signals(:, time_idx), array_elements, array_elements);

% Plot clean sound pressure distribution
figure;
subplot(1,2,1);
imagesc(sound_pressure_clean);
title('Clean Sound Pressure Distribution (t = 0.1s)');
colorbar;
axis equal tight;
colormap jet;

% Plot noisy sound pressure distribution
subplot(1,2,2);
imagesc(sound_pressure_noisy);
title(sprintf('Noisy Sound Pressure (SNR = %d dB)', SNR_dB));
colorbar;
axis equal tight;
colormap jet;

%% Save simulated signals for use in beamforming
save('simulated_microphone_signals.mat', ...
    'clean_signals', 'noisy_signals', 'source_signal', 'SNR_dB', 'fs', 't');

disp('Microphone signal simulation completed and saved.');