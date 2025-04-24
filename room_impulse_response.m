% Room Impulse Response Implementation for SC-DAMAS Simulation
% This script implements equations (3), (4), (5), and (6) from the paper
% to calculate the room impulse response function using the mirror source method

% Load the room geometry setup
load('room_geometry_setup.mat');

%% Define simulation parameters
c = 343;             % Speed of sound (m/s)
fs = 16000;          % Sampling frequency (Hz)
max_order = 3;       % Maximum reflection order
freq_analysis = 2000; % Frequency for analysis (Hz) - can be changed to match paper's examples

% Time vector for impulse response
t_max = 0.5;         % Maximum time for impulse response (s)
t = 0:1/fs:t_max;    % Time vector
n_samples = length(t);

% Initialize impulse response matrix for all microphones
h = zeros(numMics, n_samples);

% Find the reflection coefficients at the analysis frequency
[~, freq_idx] = min(abs(freq - freq_analysis));
beta_x1 = reflection_coef(1, freq_idx); % Front wall
beta_x2 = reflection_coef(3, freq_idx); % Right wall
beta_y1 = reflection_coef(6, freq_idx); % Back wall
beta_y2 = reflection_coef(1, freq_idx); % Left wall (using front wall coefficient)
beta_z1 = reflection_coef(5, freq_idx); % Ground
beta_z2 = reflection_coef(4, freq_idx); % Ceiling

%% Calculate impulse responses using mirror source method
disp('Calculating room impulse responses...');

% For each microphone
for mic = 1:numMics
    mic_pos = micPositions(mic, :);
    
    % Direct sound contribution
    r_direct = norm(mic_pos - sourcePos);
    delay_samples = round(r_direct/c * fs);
    
    % Make sure delay is within range
    if delay_samples < n_samples
        amplitude = 1/(4*pi*r_direct);
        h(mic, delay_samples+1) = h(mic, delay_samples+1) + amplitude;
    end
    
    % Mirror sources contribution - implementing Eq. (3) from the paper
    for nx = -max_order:max_order
        for ny = -max_order:max_order
            for nz = -max_order:max_order
                % Skip the direct sound (already handled)
                if nx == 0 && ny == 0 && nz == 0
                    continue;
                end
                
                % Calculate the mirror source position
                mirror_x = sourcePos(1) + 2*nx*roomDim(1);
                mirror_y = sourcePos(2) + 2*ny*roomDim(2);
                mirror_z = sourcePos(3) + 2*nz*roomDim(3);
                
                % Calculate the distance from mirror source to microphone
                r = norm([mirror_x - mic_pos(1), mirror_y - mic_pos(2), mirror_z - mic_pos(3)]);
                
                if r > 0  % Avoid singularity
                    % Calculate reflection coefficients based on reflections
                    % For the x-axis: use beta_x1 for reflections from wall near origin,
                    % use beta_x2 for reflections from wall far from origin
                    % Similarly for y and z axes
                    if nx < 0
                        refl_x = beta_x1^abs(nx); % Wall near origin (x=0)
                    else
                        refl_x = beta_x2^abs(nx); % Wall far from origin (x=Lx)
                    end
                    
                    if ny < 0
                        refl_y = beta_y1^abs(ny); % Wall near origin (y=0)
                    else
                        refl_y = beta_y2^abs(ny); % Wall far from origin (y=Ly)
                    end
                    
                    if nz < 0
                        refl_z = beta_z1^abs(nz); % Wall near origin (z=0)
                    else
                        refl_z = beta_z2^abs(nz); % Wall far from origin (z=Lz)
                    end
                    
                    % Combined reflection coefficient
                    refl_coef = refl_x * refl_y * refl_z;
                    
                    delay_samples = round(r/c * fs);
                    
                    % Make sure delay is within range
                    if delay_samples < n_samples
                        amplitude = refl_coef/(4*pi*r);
                        h(mic, delay_samples+1) = h(mic, delay_samples+1) + amplitude;
                    end
                end
            end
        end
    end
    
    % Show progress
    if mod(mic, 10) == 0
        fprintf('Processed %d of %d microphones\n', mic, numMics);
    end
end

%% Calculate frequency response using FFT (Eq. 4)
% Calculate FFT for each impulse response
H = fft(h, n_samples, 2);
freq_vector = (0:n_samples-1)*fs/n_samples;

%% Visualize results
% Plot impulse response for center microphone
center_mic = ceil(numMics/2);

figure;
subplot(2,1,1);
plot(t*1000, h(center_mic,:));
title(sprintf('Impulse Response at Center Microphone (x=%.1f, y=%.1f, z=%.1f)', ...
    micPositions(center_mic,1), micPositions(center_mic,2), micPositions(center_mic,3)));
xlabel('Time (ms)');
ylabel('Amplitude');
grid on;

% Plot magnitude of frequency response
subplot(2,1,2);
% Only plot up to Nyquist frequency
plot_freq = freq_vector(1:n_samples/2);
plot_H = abs(H(center_mic, 1:n_samples/2));
semilogx(plot_freq, 20*log10(plot_H));
title('Frequency Response at Center Microphone');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
grid on;
xlim([20, fs/2]);

% Plot impulse responses for all microphones as an image
figure;
imagesc(t*1000, 1:numMics, h);
title('Impulse Responses for All Microphones');
xlabel('Time (ms)');
ylabel('Microphone #');
colorbar;
colormap jet;

% Plot 3D visualization of reflection paths for the center microphone
figure;
hold on;

% Plot room
x = [0 0 0 0 roomDim(1) roomDim(1) roomDim(1) roomDim(1)];
y = [0 0 roomDim(2) roomDim(2) 0 0 roomDim(2) roomDim(2)];
z = [0 roomDim(3) 0 roomDim(3) 0 roomDim(3) 0 roomDim(3)];
faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
patch('Vertices', [x(:) y(:) z(:)], 'Faces', faces, ...
    'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

% Plot sound source position
scatter3(sourcePos(1), sourcePos(2), sourcePos(3), 100, 'r', 'filled');
text(sourcePos(1)+0.1, sourcePos(2), sourcePos(3), 'Sound Source', 'Color', 'r');

% Plot center microphone position
scatter3(micPositions(center_mic,1), micPositions(center_mic,2), micPositions(center_mic,3), 80, 'b', 'filled');
text(micPositions(center_mic,1)+0.1, micPositions(center_mic,2), micPositions(center_mic,3), 'Center Mic', 'Color', 'b');

impulse_responses = h;  % Create a copy with a different name
% Plot first-order mirror sources and their paths (for visualization clarity)
% We'll only show a few primary reflections for better visualization
for nx = -1:1
    for ny = -1:1
        for nz = -1:1
            % Skip the direct sound
            if nx == 0 && ny == 0 && nz == 0
                continue;
            end
            
            % Calculate mirror source position
            mirror_x = sourcePos(1) + 2*nx*roomDim(1);
            mirror_y = sourcePos(2) + 2*ny*roomDim(2);
            mirror_z = sourcePos(3) + 2*nz*roomDim(3);
            
            % Plot mirror source
            scatter3(mirror_x, mirror_y, mirror_z, 50, 'm', 'filled', 'MarkerFaceAlpha', 0.5);
            
            % Plot path from mirror source to microphone
            plot3([mirror_x, micPositions(center_mic,1)], ...
                  [mirror_y, micPositions(center_mic,2)], ...
                  [mirror_z, micPositions(center_mic,3)], 'm--', 'LineWidth', 0.5);
            % Set transparency (compatible with all MATLAB versions)
            h = findobj(gca, 'Type', 'line', 'Color', 'm');
            if ~isempty(h)
                set(h(end), 'Color', [1 0 1 0.3]);  % magenta with 0.3 alpha
            end
        end
    end
end

% Plot direct path
plot3([sourcePos(1), micPositions(center_mic,1)], ...
      [sourcePos(2), micPositions(center_mic,2)], ...
      [sourcePos(3), micPositions(center_mic,3)], 'g-', 'LineWidth', 2);

% Set plot properties
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Room with Sound Source, Microphone, and First-order Mirror Sources');
grid on;
axis equal;
view(30, 20);

% Save the calculated impulse and frequency responses
% IMPORTANT: Make a copy of h to avoid saving the plot handle

save('room_impulse_response.mat', 'impulse_responses', 'H', 'freq_vector', 't', 'c', 'fs', 'freq_analysis');

disp('Room impulse responses calculated and saved successfully.');