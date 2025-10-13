% 3D Image Source Method Visualization
% Shows the cube room and virtual source positions as dots

clear; clc; close all;

% Fixed room dimensions (meters)
Lx = 4.0;  % Length
Ly = 3.0;  % Width  
Lz = 2.5;  % Height

% Source position (center of room)
source_pos = [Lx/2, Ly/2, Lz/2];  % [2.0, 1.5, 1.25]

% Create figure
figure('Position', [100, 100, 1200, 800]);

% Plot the room as a wireframe cube
hold on;

% Define cube vertices
x_room = [0 0 0 0 Lx Lx Lx Lx];
y_room = [0 0 Ly Ly 0 0 Ly Ly];
z_room = [0 Lz 0 Lz 0 Lz 0 Lz];

% Define cube faces
faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];

% Plot room wireframe
patch('Vertices', [x_room(:) y_room(:) z_room(:)], 'Faces', faces, ...
      'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2);

% Plot original source
plot3(source_pos(1), source_pos(2), source_pos(3), 'ro', ...
      'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
text(source_pos(1)+0.2, source_pos(2), source_pos(3), 'Original Source', ...
     'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');

% Calculate and plot virtual sources for different orders
colors = {'b', 'g', 'm', 'c', 'y'};  % Different colors for different orders
order_labels = {};
order_counts = [];

for nmax = 1:1
    virtual_sources = [];
    path_info = {};
    
    % Calculate all virtual sources for this order
    path_counter = 1;
    for nx = -nmax:nmax
        for ny = -nmax:nmax  
            for nz = -nmax:nmax
                
                % Skip direct path [0,0,0]
                if nx == 0 && ny == 0 && nz == 0
                    continue;
                end
                
                % Skip paths exceeding maximum order
                current_order = abs(nx) + abs(ny) + abs(nz);
                if current_order > nmax
                    continue;
                end
                
                % Calculate spatial offset (THE KEY FORMULA)
                offset = [2*nx*Lx, 2*ny*Ly, 2*nz*Lz];
                
                % Virtual source position  
                virtual_source_pos = source_pos + offset;
                
                % Store for plotting
                virtual_sources = [virtual_sources; virtual_source_pos, current_order];
                path_info{path_counter} = sprintf('[%d,%d,%d]', nx, ny, nz);
                path_counter = path_counter + 1;
            end
        end
    end
    
    % Plot virtual sources by order
    for order = 1:nmax
        order_mask = virtual_sources(:,4) == order;
        if any(order_mask)
            sources_this_order = virtual_sources(order_mask, 1:3);
            
            plot3(sources_this_order(:,1), sources_this_order(:,2), sources_this_order(:,3), ...
                  'o', 'Color', colors{order}, 'MarkerSize', 8, ...
                  'MarkerFaceColor', colors{order}, 'MarkerEdgeColor', 'k');
            
            order_labels{end+1} = sprintf('Order %d', order);
            order_counts(end+1) = size(sources_this_order, 1);
        end
    end
end

% Set axis properties
xlabel('X (meters)', 'FontSize', 12);
ylabel('Y (meters)', 'FontSize', 12);
zlabel('Z (meters)', 'FontSize', 12);
title('3D Image Source Method - Virtual Source Positions', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
axis equal;

% Set viewing limits to show virtual sources
x_lim = [-10, 14];
y_lim = [-7, 10];
z_lim = [-5, 8];
xlim(x_lim); ylim(y_lim); zlim(z_lim);

% Add legend
legend_entries = {'Original Source', order_labels{:}};
legend(legend_entries, 'Location', 'best', 'FontSize', 10);

% Set viewing angle
view(45, 30);

% Add room dimension annotations
text(Lx/2, -0.5, -0.5, sprintf('%.1fm', Lx), 'HorizontalAlignment', 'center', 'FontSize', 10);
text(-0.5, Ly/2, -0.5, sprintf('%.1fm', Ly), 'HorizontalAlignment', 'center', 'FontSize', 10);
text(-0.5, -0.5, Lz/2, sprintf('%.1fm', Lz), 'HorizontalAlignment', 'center', 'FontSize', 10);

% Create a second figure showing specific examples
figure('Position', [200, 200, 1400, 600]);

% Subplot 1: Order 1 paths only
subplot(1,3,1);
hold on;

% Plot room
patch('Vertices', [x_room(:) y_room(:) z_room(:)], 'Faces', faces, ...
      'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

% Plot original source
plot3(source_pos(1), source_pos(2), source_pos(3), 'ro', ...
      'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Calculate order 1 virtual sources with labels
examples_1st_coords = [-1, 0, 0; 1, 0, 0; 0, -1, 0; 0, 1, 0; 0, 0, -1; 0, 0, 1];
examples_1st_labels = {'Left Wall', 'Right Wall', 'Back Wall', 'Front Wall', 'Floor', 'Ceiling'};

for i = 1:size(examples_1st_coords, 1)
    nx = examples_1st_coords(i, 1);
    ny = examples_1st_coords(i, 2);
    nz = examples_1st_coords(i, 3);
    label = examples_1st_labels{i};
    
    offset = [2*nx*Lx, 2*ny*Ly, 2*nz*Lz];
    virtual_pos = source_pos + offset;
    
    plot3(virtual_pos(1), virtual_pos(2), virtual_pos(3), 'bo', ...
          'MarkerSize', 8, 'MarkerFaceColor', 'b');
    text(virtual_pos(1), virtual_pos(2), virtual_pos(3)+0.3, ...
         sprintf('%s\n[%d,%d,%d]', label, nx, ny, nz), ...
         'HorizontalAlignment', 'center', 'FontSize', 8);
end

title('Order 1 Reflections (6 paths)', 'FontSize', 12);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
xlim([-10, 14]); ylim([-7, 10]); zlim([-6, 8]);
grid on; axis equal; view(45, 30);

% Subplot 2: Specific examples with coordinates
subplot(1,3,2);
hold on;

% Plot room
patch('Vertices', [x_room(:) y_room(:) z_room(:)], 'Faces', faces, ...
      'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);

% Plot original source
plot3(source_pos(1), source_pos(2), source_pos(3), 'ro', ...
      'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Specific examples
examples_coords = [-1, 0, 0; -2, 0, 0; 1, 1, 0; -1, 1, -1];
examples_colors = {'b', 'g', 'm', 'c'};

for i = 1:size(examples_coords, 1)
    nx = examples_coords(i, 1);
    ny = examples_coords(i, 2);
    nz = examples_coords(i, 3);
    color = examples_colors{i};
    
    offset = [2*nx*Lx, 2*ny*Ly, 2*nz*Lz];
    virtual_pos = source_pos + offset;
    
    plot3(virtual_pos(1), virtual_pos(2), virtual_pos(3), 'o', ...
          'Color', color, 'MarkerSize', 10, 'MarkerFaceColor', color);
    text(virtual_pos(1), virtual_pos(2), virtual_pos(3)+0.3, ...
         sprintf('[%d,%d,%d]\n(%.1f,%.1f,%.1f)', nx, ny, nz, ...
         virtual_pos(1), virtual_pos(2), virtual_pos(3)), ...
         'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', color);
end

title('Key Examples with Coordinates', 'FontSize', 12);
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
xlim([-10, 14]); ylim([-7, 10]); zlim([-6, 8]);
grid on; axis equal; view(45, 30);

% Subplot 3: 2D Cross-section (X-Y plane at Z = source height)
subplot(1,3,3);
hold on;

% Plot room outline (top view)
rectangle('Position', [0, 0, Lx, Ly], 'EdgeColor', 'k', 'LineWidth', 2);

% Plot original source
plot(source_pos(1), source_pos(2), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

% Plot virtual sources in X-Y plane
for nx = -2:2
    for ny = -2:2
        if nx == 0 && ny == 0, continue; end
        if abs(nx) + abs(ny) > 2, continue; end
        
        offset = [2*nx*Lx, 2*ny*Ly, 0];
        virtual_pos = source_pos + offset;
        
        order = abs(nx) + abs(ny);
        plot(virtual_pos(1), virtual_pos(2), 'o', ...
             'Color', colors{order}, 'MarkerSize', 8, 'MarkerFaceColor', colors{order});
        text(virtual_pos(1)+0.2, virtual_pos(2), sprintf('[%d,%d,0]', nx, ny), ...
             'FontSize', 8);
    end
end

title('Top View (X-Y plane)', 'FontSize', 12);
xlabel('X (meters)'); ylabel('Y (meters)');
xlim([-10, 14]); ylim([-7, 10]);
grid on; axis equal;

% Add wall labels
text(Lx/2, -0.5, 'Front Wall', 'HorizontalAlignment', 'center');
text(Lx/2, Ly+0.3, 'Back Wall', 'HorizontalAlignment', 'center');
text(-0.8, Ly/2, 'Left Wall', 'HorizontalAlignment', 'center', 'Rotation', 90);
text(Lx+0.3, Ly/2, 'Right Wall', 'HorizontalAlignment', 'center', 'Rotation', 90);

% Print numerical results
fprintf('=== NUMERICAL RESULTS ===\n');
fprintf('Room: %.1f × %.1f × %.1f meters\n', Lx, Ly, Lz);
fprintf('Source: (%.1f, %.1f, %.1f)\n\n', source_pos(1), source_pos(2), source_pos(3));

key_examples_coords = [-1, 0, 0; 1, 0, 0; -2, 0, 0; 1, 1, 0; -1, 1, -1];
key_examples_desc = {
    'One reflection off LEFT wall';
    'One reflection off RIGHT wall';
    'Two reflections off LEFT wall';
    'Corner: RIGHT + FRONT walls';
    'Triple: LEFT + FRONT + FLOOR'
};

for i = 1:size(key_examples_coords, 1)
    nx = key_examples_coords(i, 1);
    ny = key_examples_coords(i, 2);
    nz = key_examples_coords(i, 3);
    description = key_examples_desc{i};
    
    offset = [2*nx*Lx, 2*ny*Ly, 2*nz*Lz];
    virtual_pos = source_pos + offset;
    distance = norm(virtual_pos - source_pos);
    
    fprintf('%s:\n', description);
    fprintf('  [nx,ny,nz] = [%d,%d,%d]\n', nx, ny, nz);
    fprintf('  Offset = [%.1f, %.1f, %.1f]\n', offset(1), offset(2), offset(3));
    fprintf('  Virtual source = (%.1f, %.1f, %.1f)\n', virtual_pos(1), virtual_pos(2), virtual_pos(3));
    fprintf('  Distance from original = %.1f meters\n\n', distance);
end