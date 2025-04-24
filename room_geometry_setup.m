% Room Geometry Setup for SC-DAMAS Simulation
% This script implements the first step of the simulation workflow
% for acoustic imaging of low SNR sound sources in reverberant fields

clear all;
close all;
clc;

%% Define Room Dimensions (in meters)
roomDim = [3, 4, 3];  % [Length (x), Width (y), Height (z)]

%% Define Sound Source Position (in meters)
sourcePos = [1.5, 3, 1.5];  % Position of the monopole source

%% Define Sensor Array
array_dist = 1;     % Distance from sound source (m)
array_size = 1;     % Array size (m x m)
array_elements = 11; % Number of elements in each direction (11x11)
element_spacing = 0.1; % Spacing between array elements (m)

% Calculate array position based on source position (1m in front of source)
% Assuming array is placed at y=2 (1m from source at y=3)
arrayCenter = [sourcePos(1), sourcePos(2)-array_dist, sourcePos(3)];

% Calculate array element positions
[x_grid, z_grid] = meshgrid(-(array_size/2):element_spacing:(array_size/2), ...
                           -(array_size/2):element_spacing:(array_size/2));
                       
array_x = arrayCenter(1) + x_grid;
array_y = repmat(arrayCenter(2), size(x_grid));
array_z = arrayCenter(3) + z_grid;

% Store all microphone positions
numMics = array_elements^2;
micPositions = zeros(numMics, 3);
counter = 1;
for i = 1:size(array_x, 1)
    for j = 1:size(array_x, 2)
        micPositions(counter, :) = [array_x(i,j), array_y(i,j), array_z(i,j)];
        counter = counter + 1;
    end
end

%% Define Acoustic Properties of Walls
% Absorption coefficients from Table 2 in the paper
freq = [125, 250, 500, 1000, 2000, 4000]; % Frequencies in Hz

% Define absorption coefficients for each wall
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

% Calculate reflection coefficients using formula β = ±√(1-γ)
% We'll use positive values for reflection coefficients
reflection_coef = sqrt(1 - absorption_coef);

%% Visualization of the Geometry
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

% Plot microphone array
scatter3(micPositions(:,1), micPositions(:,2), micPositions(:,3), 25, 'b', 'filled');
text(arrayCenter(1), arrayCenter(2)-0.2, arrayCenter(3), 'Microphone Array', 'Color', 'b');

% Add room dimensions text
text(roomDim(1)/2, -0.3, 0, ['Length: ' num2str(roomDim(1)) ' m'], 'Color', 'k');
text(-0.3, roomDim(2)/2, 0, ['Width: ' num2str(roomDim(2)) ' m'], 'Color', 'k');
text(0, -0.3, roomDim(3)/2, ['Height: ' num2str(roomDim(3)) ' m'], 'Color', 'k');

% Set plot properties
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Room Geometry for SC-DAMAS Simulation');
grid on;
axis equal;
view(30, 20);

% Set axis limits with some padding
axis([-0.5 roomDim(1)+0.5 -0.5 roomDim(2)+0.5 -0.5 roomDim(3)+0.5]);

% Display a summary of the setup
fprintf('Room Dimensions: %.1f m × %.1f m × %.1f m\n', roomDim(1), roomDim(2), roomDim(3));
fprintf('Sound Source Position: (%.1f, %.1f, %.1f)\n', sourcePos(1), sourcePos(2), sourcePos(3));
fprintf('Microphone Array: %d × %d elements (%.1f m × %.1f m)\n', ...
    array_elements, array_elements, array_size, array_size);
fprintf('Number of microphones: %d\n', numMics);
fprintf('Distance from source to array: %.1f m\n', array_dist);

% Save the workspace variables for next steps
save('room_geometry_setup.mat');

disp('Room geometry setup completed and saved.');