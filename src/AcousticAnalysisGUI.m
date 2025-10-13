
function AcousticAnalysisGUI()
% ACOUSTICANALYSISGUI - Comprehensive GUI for Acoustic Simulation & Analysis System
% 
% This GUI provides a complete interface for:
% - Room acoustic simulation with Doppler effects
% - Deep learning source localization 
% - Impulse response analysis and comparison
% - Advanced beamforming algorithms
% - Real-time measurement and recording
%
% The GUI orchestrates all existing modules without requiring code changes
% Author: Advanced Acoustic Analysis System
% Version: 2.0 Professional


    % Initialize main figure
    fig = uifigure('Name', 'Advanced Acoustic Analysis System v2.0', ...
                   'Position', [100, 100, 1400, 900], ...
                   'Color', [0.94, 0.94, 0.94], ...
                   'Resize', 'on');
    
    % Create main layout
    mainGrid = uigridlayout(fig, [1, 2]);
    mainGrid.ColumnWidth = {220, '1x'};
    
    % Left panel for navigation and controls
    leftPanel = uipanel(mainGrid, 'Title', 'Control Panel', ...
                        'BackgroundColor', [0.9, 0.9, 0.9], ...
                        'FontWeight', 'bold');
    
    % Right panel for main content
    rightPanel = uipanel(mainGrid, 'Title', 'Analysis Workspace', ...
                         'BackgroundColor', 'white', ...
                         'FontWeight', 'bold');
    
    % Create navigation buttons
    navGrid = uigridlayout(leftPanel, [8, 1]);
    navGrid.RowHeight = {50, 50, 50, 50, 50, 50, 50, '1x'};
    
    % Navigation buttons
    btnSimulation = uibutton(navGrid, 'Text', '🎯 Simulation Setup', ...
                             'FontSize', 12, 'FontWeight', 'bold', ...
                             'BackgroundColor', [0.2, 0.6, 0.9], ...
                             'FontColor', 'white');
    
    btnDNN = uibutton(navGrid, 'Text', '🧠 Deep Learning', ...
                      'FontSize', 12, 'FontWeight', 'bold', ...
                      'BackgroundColor', [0.9, 0.4, 0.2], ...
                      'FontColor', 'white');
    
    btnAnalysis = uibutton(navGrid, 'Text', '📊 IR Analysis', ...
                           'FontSize', 12, 'FontWeight', 'bold', ...
                           'BackgroundColor', [0.2, 0.8, 0.4], ...
                           'FontColor', 'white');
    
    btnBeamforming = uibutton(navGrid, 'Text', '📡 Beamforming', ...
                              'FontSize', 12, 'FontWeight', 'bold', ...
                              'BackgroundColor', [0.8, 0.2, 0.8], ...
                              'FontColor', 'white');
    
    
    
 
    
 
    
    % Status panel
    statusPanel = uipanel(navGrid, 'Title', 'System Status', ...
                          'BackgroundColor', [0.95, 0.95, 0.95]);
    statusGrid = uigridlayout(statusPanel, [4, 1]);
    statusGrid.RowHeight = {'fit', 'fit', 'fit', '1x'};
    
    statusText = uilabel(statusGrid, 'Text', '✅ System Ready', ...
                         'FontColor', [0, 0.6, 0], 'FontWeight', 'bold');
    % progressBar = uiprogressdlg(fig, 'Title', 'Processing...', 'Indeterminate', 'on', 'Visible', 'off');
    progressBar = [];  % Initialize as empty, create when needed
    
    % Create content panels (initially hidden)
    contentPanels = containers.Map();
    
    %% ============= SIMULATION SETUP PANEL =============
    simPanel = createSimulationPanel(rightPanel);
    contentPanels('simulation') = simPanel;
    
    %% ============= DEEP LEARNING PANEL =============
    dnnPanel = createDNNPanel(rightPanel);
    contentPanels('dnn') = dnnPanel;
    
    %% ============= IR ANALYSIS PANEL =============
    analysisPanel = createAnalysisPanel(rightPanel);
    contentPanels('analysis') = analysisPanel;
    
    %% ============= BEAMFORMING PANEL =============
    beamPanel = createBeamformingPanel(rightPanel);
    contentPanels('beamforming') = beamPanel;
    
    


    % Initially show simulation panel
    showPanel('simulation');
    % Initially show simulation panel
    showPanel('simulation');
    
    % Set button callbacks
    btnSimulation.ButtonPushedFcn = @(~,~) showPanel('simulation');
    btnDNN.ButtonPushedFcn = @(~,~) showPanel('dnn');
    btnAnalysis.ButtonPushedFcn = @(~,~) showPanel('analysis');
    btnBeamforming.ButtonPushedFcn = @(~,~) showPanel('beamforming');


    
    % Store GUI data
    guiData = struct();
    guiData.figure = fig;
    guiData.statusText = statusText;
    guiData.progressBar = progressBar;
    guiData.panels = contentPanels;
    fig.UserData = guiData;
    
    %% ============= HELPER FUNCTIONS =============
    
    function showPanel(panelName)
        fprintf('Showing panel: %s\n', panelName);
        fprintf('Panel exists: %d\n', isKey(contentPanels, panelName));
        % Hide all panels
        panels = values(contentPanels);
        for i = 1:length(panels)
            panels{i}.Visible = 'off';
        end
        
        % Show selected panel
        if isKey(contentPanels, panelName)
            panel = contentPanels(panelName);
            panel.Visible = 'on';
            % ADD THESE DEBUG LINES:
            fprintf('Panel visibility: %s\n', panel.Visible);
            fprintf('Panel position: [%g %g %g %g]\n', panel.Position);
            fprintf('Panel parent class: %s\n', class(panel.Parent));
        end
        
        % Update status
        switch panelName
            case 'simulation'
                updateStatus('🎯 Simulation Setup Active', [0, 0.6, 0]);
            case 'dnn'
                updateStatus('🧠 Deep Learning Module Active', [0.8, 0.4, 0]);
            case 'analysis'
                updateStatus('📊 IR Analysis Module Active', [0, 0.6, 0.2]);
            case 'beamforming'
                updateStatus('📡 Beamforming Module Active', [0.6, 0, 0.6]);
            
            
        end
    end
    
    function updateStatus(message, color)
        statusText.Text = message;
        statusText.FontColor = color;
    end
    
    function showProgress(title, message)
        if isempty(progressBar) || ~isvalid(progressBar)
            progressBar = uiprogressdlg(fig, 'Title', title, 'Message', message, 'Indeterminate', 'on');
        else
            progressBar.Title = title;
            progressBar.Message = message;
        end
        drawnow;
    end
    
    function hideProgress()
        if ~isempty(progressBar) && isvalid(progressBar)
            close(progressBar);
            progressBar = [];
        end
        drawnow;
    end
end




%% ============= SIMULATION SETUP PANEL =============
function panel = createSimulationPanel(parent)
    panel = uipanel(parent, 'BackgroundColor', 'white', 'Visible', 'off');
    
    grid = uigridlayout(panel, [4, 2]);
    grid.RowHeight = {120, 220, '1x', 60}; 
    grid.ColumnWidth = {'1x', '1x'};
    
    % Title and description
    titlePanel = uipanel(grid, 'BackgroundColor', [0.9, 0.95, 1]);
    titlePanel.Layout.Row = 1;
    titlePanel.Layout.Column = [1, 2];
    
    titleGrid = uigridlayout(titlePanel, [3, 1]);
    uilabel(titleGrid, 'Text', '🎯 ACOUSTIC SIMULATION SETUP', ...
            'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    uilabel(titleGrid, 'Text', 'Configure room geometry, microphone arrays, and simulation parameters', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    uilabel(titleGrid, 'Text', 'Integrates with FastHybridDopplerReverbSimulator for realistic acoustic modeling', ...
        'FontSize', 10, 'FontAngle', 'italic', 'HorizontalAlignment', 'center');
    
    % Room Configuration Panel
    roomPanel = uipanel(grid, 'Title', '🏠 Room Configuration', ...
                        'FontWeight', 'bold', 'BackgroundColor', [0.98, 0.98, 1]);
    roomPanel.Layout.Row = 2;
    roomPanel.Layout.Column = 1;
    
    roomGrid = uigridlayout(roomPanel, [6, 2]);
    roomGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    
    uilabel(roomGrid, 'Text', 'Room Type:');
    roomTypeDD = uidropdown(roomGrid, 'Items', {'Rectangular', 'Cylindrical'}, 'Value', 'Rectangular');
    
    uilabel(roomGrid, 'Text', 'Length (m):');
    roomLengthSpinner = uispinner(roomGrid, 'Value', 4, 'Limits', [1, 50], 'Step', 0.1);
    
    uilabel(roomGrid, 'Text', 'Width (m):');
    roomWidthSpinner = uispinner(roomGrid, 'Value', 4, 'Limits', [1, 50], 'Step', 0.1);
    
    uilabel(roomGrid, 'Text', 'Height (m):');
    roomHeightSpinner = uispinner(roomGrid, 'Value', 4, 'Limits', [1, 20], 'Step', 0.1);
    
    uilabel(roomGrid, 'Text', 'Temperature (°C):');
    tempSpinner = uispinner(roomGrid, 'Value', 20, 'Limits', [-20, 50], 'Step', 1);
    
    uilabel(roomGrid, 'Text', 'Humidity (%):');
    humiditySpinner = uispinner(roomGrid, 'Value', 50, 'Limits', [0, 100], 'Step', 1);
    
    % Microphone Array Panel
    micPanel = uipanel(grid, 'Title', '🎙️ Microphone Array', ...
                       'FontWeight', 'bold', 'BackgroundColor', [0.98, 1, 0.98]);
    micPanel.Layout.Row = 2;
    micPanel.Layout.Column = 2;
    
    micGrid = uigridlayout(micPanel, [6, 2]);
    micGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    
    uilabel(micGrid, 'Text', 'Array Type:');
    arrayTypeDD = uidropdown(micGrid, 'Items', {'Linear', 'Rectangular Grid', 'Circular', 'Spiral', 'Spherical'}, 'Value', 'Linear');
    
    uilabel(micGrid, 'Text', 'Number of Mics:');
    numMicsSpinner = uispinner(micGrid, 'Value', 2, 'Limits', [1, 64], 'Step', 1);
    
    % Single Array Center field
    uilabel(micGrid, 'Text', 'Array Center [X,Y,Z] (m):');
    arrayCenterEdit = uieditfield(micGrid, 'Value', '0.5, 0.5, 0.5', ...
                                  'BackgroundColor', [1, 1, 1], ...
                                  'Tooltip', 'Enter coordinates as: X, Y, Z (e.g., 0.5, 0.5, 0.5)');
    
    uilabel(micGrid, 'Text', 'Orientation:');
    orientationDD = uidropdown(micGrid, 'Items', {'X-axis', 'Y-axis', 'Z-axis'}, 'Value', 'Y-axis');
    
    % Source Configuration Panel
    sourcePanel = uipanel(grid, 'Title', '🔊 Source Configuration', ...
                          'FontWeight', 'bold', 'BackgroundColor', [1, 0.98, 0.98]);
    sourcePanel.Layout.Row = 3;
    sourcePanel.Layout.Column = 1;
    
    sourceGrid = uigridlayout(sourcePanel, [7, 2]);
    sourceGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    
    uilabel(sourceGrid, 'Text', 'Signal Type:');
    signalTypeDD = uidropdown(sourceGrid, 'Items', {'Measurement Sweep', 'Chirp', 'Siren'}, 'Value', 'Measurement Sweep');
    
    % uilabel(sourceGrid, 'Text', 'Start Position X (m):');
    % startXSpinner = uispinner(sourceGrid, 'Value', 2, 'Limits', [0, 10], 'Step', 0.1);
    % 
    % uilabel(sourceGrid, 'Text', 'Start Position Y (m):');
    % startYSpinner = uispinner(sourceGrid, 'Value', 2, 'Limits', [0, 10], 'Step', 0.1);
    % 
    % uilabel(sourceGrid, 'Text', 'Start Position Z (m):');
    % startZSpinner = uispinner(sourceGrid, 'Value', 2, 'Limits', [0, 10], 'Step', 0.1);

    % Single End Position field (initially disabled)
    uilabel(sourceGrid, 'Text', 'Start Position [X,Y,Z] (m):');
    startPositionEdit = uieditfield(sourceGrid, 'Value', '2, 2, 2', 'Enable', 'off', ...
                                  'BackgroundColor', [0.94, 0.94, 0.94], ...
                                  'Tooltip', 'Enter coordinates as: X, Y, Z (e.g., 3.5, 2.1, 1.8)');

    % Single End Position field (initially disabled)
    uilabel(sourceGrid, 'Text', 'End Position [X,Y,Z] (m):');
    endPositionEdit = uieditfield(sourceGrid, 'Value', '2, 2, 2', 'Enable', 'off', ...
                                  'BackgroundColor', [0.94, 0.94, 0.94], ...
                                  'Tooltip', 'Enter coordinates as: X, Y, Z (e.g., 3.5, 2.1, 1.8)');
    
    uilabel(sourceGrid, 'Text', 'Duration (s):');
    durationSpinner = uispinner(sourceGrid, 'Value', 5, 'Limits', [0.1, 10], 'Step', 0.1);
    
    uilabel(sourceGrid, 'Text', 'Source Power (W):');
    powerSpinner = uispinner(sourceGrid, 'Value', 0.5, 'Limits', [0.01, 10], 'Step', 0.01);
    
    % Moving source checkbox and end position
    movingSourceCB = uicheckbox(sourceGrid, 'Text', 'Moving Source');
    movingSourceCB.Layout.Column = [1, 2];
    
    % Simulation Options Panel
    simOptionsPanel = uipanel(grid, 'Title', '⚙️ Simulation Options', ...
                              'FontWeight', 'bold', 'BackgroundColor', [1, 1, 0.98]);
    simOptionsPanel.Layout.Row = 3;
    simOptionsPanel.Layout.Column = 2;
    
    simOptionsGrid = uigridlayout(simOptionsPanel, [7, 2]);
    simOptionsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    
    enableNoiseCB = uicheckbox(simOptionsGrid, 'Text', 'Enable Noise', 'Value', true);
    enableNoiseCB.Layout.Column = [1, 2];
    
    uilabel(simOptionsGrid, 'Text', 'Noise Type:');
    noiseTypeDD = uidropdown(simOptionsGrid, 'Items', {'White', 'Pink', 'Ambient'}, 'Value', 'White');
    
    uilabel(simOptionsGrid, 'Text', 'Max Reflections:');
    maxReflSpinner = uispinner(simOptionsGrid, 'Value', 1, 'Limits', [0, 5], 'Step', 1);
    
    uilabel(simOptionsGrid, 'Text', 'Sampling Rate (Hz):');
    fsSpinner = uispinner(simOptionsGrid, 'Value', 16000, 'Limits', [8000, 96000], 'Step', 1000);
    
    enableVisualizationCB = uicheckbox(simOptionsGrid, 'Text', 'Enable Visualization', 'Value', true);
    enableVisualizationCB.Layout.Column = [1, 2];
    
    freqDomainReflCB = uicheckbox(simOptionsGrid, 'Text', 'Frequency Domain Reflections');
    freqDomainReflCB.Layout.Column = [1, 2];
    
    vehicleBodyCB = uicheckbox(simOptionsGrid, 'Text', 'Vehicle Body Reflections');
    vehicleBodyCB.Layout.Column = [1, 2];
    
    % Control Buttons Panel
    controlPanel = uipanel(grid, 'BackgroundColor', [0.95, 0.95, 0.95]);
    controlPanel.Layout.Row = 4;
    controlPanel.Layout.Column = [1, 2];
    
    controlGrid = uigridlayout(controlPanel, [1, 4]);
    
    controlGrid.ColumnWidth = {'1x', 'fit', 'fit', '1x'};
    
    
    uilabel(controlGrid, 'Text', ''); % Spacer
    
    runSimBtn = uibutton(controlGrid, 'Text', '🚀 Run Simulation', ...
                         'FontSize', 14, 'FontWeight', 'bold', ...
                         'BackgroundColor', [0.2, 0.7, 0.3], ...
                         'FontColor', 'white');
    
    previewBtn = uibutton(controlGrid, 'Text', '👁️ Preview Setup', ...
                          'FontSize', 12, ...
                          'BackgroundColor', [0.3, 0.6, 0.9], ...
                          'FontColor', 'white');

    movingSourceCB.ValueChangedFcn = @(src, ~) toggleEndPosition(src.Value);
    startPositionEdit.ValueChangedFcn = @(src, ~) syncPositions();
    
    % Nested functions for callbacks
    function toggleEndPosition(isMoving)
        if isMoving
            % Enable end position field
            endPositionEdit.Enable = 'on';
            endPositionEdit.BackgroundColor = [1, 1, 1]; % White
            endPositionEdit.FontColor = [0, 0, 0]; % Black text
        else
            % Disable end position and auto-sync with start
            endPositionEdit.Enable = 'off';
            endPositionEdit.BackgroundColor = [0.94, 0.94, 0.94]; % Gray
            endPositionEdit.FontColor = [0.5, 0.5, 0.5]; % Gray text
            endPositionEdit.Value = startPositionEdit.Value; % Auto-sync
        end
    end
    
    function syncPositions()
        % If moving source is disabled, sync end position with start position
        if ~movingSourceCB.Value
            endPositionEdit.Value = startPositionEdit.Value;
        end
    end
    
    % Set button callbacks
    runSimBtn.ButtonPushedFcn = @(~,~) runSimulation();
    previewBtn.ButtonPushedFcn = @(~,~) previewSetup();
    
    % Store UI components for access in callbacks
    panel.UserData = struct('roomType', roomTypeDD, 'roomLength', roomLengthSpinner, ...
                       'roomWidth', roomWidthSpinner, 'roomHeight', roomHeightSpinner, ...
                       'temperature', tempSpinner, 'humidity', humiditySpinner, ...
                       'arrayType', arrayTypeDD, 'numMics', numMicsSpinner, ...
                       'arrayCenter', arrayCenterEdit, ... % REPLACED 3 array center spinners
                       'orientation', orientationDD, ...
                       'signalType', signalTypeDD, ...
                       'startPosition', startPositionEdit, ...
                       'endPosition', endPositionEdit, ...
                       'duration', durationSpinner, 'power', powerSpinner, ...
                       'movingSource', movingSourceCB, 'enableNoise', enableNoiseCB, ...
                       'noiseType', noiseTypeDD, 'maxRefl', maxReflSpinner, ...
                       'fs', fsSpinner, 'enableVis', enableVisualizationCB, ...
                       'freqDomainRefl', freqDomainReflCB, 'vehicleBody', vehicleBodyCB);
    
    function runSimulation()
        % Get GUI data
        mainGuiData = ancestor(panel, 'figure').UserData;
        mainGuiData.statusText.Text = '🚀 Running Acoustic Simulation...';
        mainGuiData.statusText.FontColor = [0.8, 0.4, 0];
        
        try
            % Show progress
            showProgress('Acoustic Simulation', 'Preparing simulation parameters...');
            
            % Collect all parameters from GUI
            params = collectSimulationParameters(panel.UserData);
            
            % Update progress
            showProgress('Acoustic Simulation', 'Running FastHybridDopplerReverbSimulator...');
            
            % Call the simulator with external parameters
            runExternalSimulatorComplete(params);
            
            % Hide progress
            hideProgress();
            
            % Update status
            mainGuiData.statusText.Text = '✅ Simulation Complete!';
            mainGuiData.statusText.FontColor = [0, 0.6, 0];
            
            % Show completion dialog
            uialert(ancestor(panel, 'figure'), ...
                   'Acoustic simulation completed successfully! Check workspace for results.', ...
                   'Simulation Complete', 'Icon', 'success');
                   
        catch ME
            hideProgress();
            mainGuiData.statusText.Text = '❌ Simulation Failed';
            mainGuiData.statusText.FontColor = [0.8, 0, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   sprintf('Simulation failed: %s', ME.message), ...
                   'Simulation Error', 'Icon', 'error');
        end
        
        function showProgress(title, message)
            mainGuiData.progressBar.Title = title;
            mainGuiData.progressBar.Message = message;
            mainGuiData.progressBar.Visible = 'on';
            drawnow;
        end
        
        function hideProgress()
            mainGuiData.progressBar.Visible = 'off';
            drawnow;
        end
    end
    
    function previewSetup()
        % Show a preview of the current setup
        params = collectSimulationParameters(panel.UserData);
        
        msg = sprintf(['Room: %s (%.1f × %.1f × %.1f m)\n' ...
                      'Array: %s with %d microphones\n' ...
                      'Source: %s at [%.1f, %.1f, %.1f]\n' ...
                      'Duration: %.1f seconds\n' ...
                      'Noise: %s (%s)'], ...
                      params.room_type, params.roomDim(1), params.roomDim(2), params.roomDim(3), ...
                      params.array_type, params.numMics, ...
                      params.signal_choice, params.start_pos(1), params.start_pos(2), params.start_pos(3), ...
                      params.duration, ...
                      char(string(params.enable_noise)), params.noise_type);
        
        uialert(ancestor(panel, 'figure'), msg, 'Setup Preview', 'Icon', 'info');
    end
end



%% ============= DEEP LEARNING PANEL =============
function panel = createDNNPanel(parent)
    panel = uipanel(parent, 'BackgroundColor', 'white', 'Visible', 'off');
    
    grid = uigridlayout(panel, [4, 2]);
    grid.RowHeight = {120, '1x', '1x', 60};  % Give dropdown 60px height
    grid.ColumnWidth = {'1x', '1x'};
    
    % Title and description
    titlePanel = uipanel(grid, 'BackgroundColor', [1, 0.95, 0.9]);
    titlePanel.Layout.Row = 1;
    titlePanel.Layout.Column = [1, 2];
    
    titleGrid = uigridlayout(titlePanel, [3, 1]);
    uilabel(titleGrid, 'Text', '🧠 DEEP LEARNING SOURCE LOCALIZATION', ...
            'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    uilabel(titleGrid, 'Text', 'Train and evaluate neural networks for 3D acoustic source localization', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    uilabel(titleGrid, 'Text', 'Supports individual (azimuth/distance/elevation) and combined 3D models', ...
            'FontSize', 10, 'FontAngle', 'italic', 'HorizontalAlignment', 'center');
    
    % Model Selection Panel
    modelPanel = uipanel(grid, 'Title', '🎯 Model Selection', ...
                         'FontWeight', 'bold', 'BackgroundColor', [1, 0.98, 0.95]);
    modelPanel.Layout.Row = 2;
    modelPanel.Layout.Column = 1;
    
    modelGrid = uigridlayout(modelPanel, [4, 1]);
    modelGrid.RowHeight = {'1x', 'fit'};  % Add this line
    modelGrid.RowSpacing = 5;             % Add this line
    
    % Individual models group - CHANGED TO CHECKBOXES for multiple selection
    individualBG = uipanel(modelGrid, 'Title', 'Individual Models');
    azimuthRB = uicheckbox(individualBG, 'Text', 'Azimuth Estimation', 'Position', [10, 60, 150, 22], 'Value', false);
    distanceRB = uicheckbox(individualBG, 'Text', 'Distance Estimation', 'Position', [10, 35, 150, 22], 'Value', false);
    elevationRB = uicheckbox(individualBG, 'Text', 'Elevation Estimation', 'Position', [10, 10, 150, 22], 'Value', false);
    
    % Combined model
    combinedCB = uicheckbox(modelGrid, 'Text', 'Combined 3D Model', 'FontWeight', 'bold');
    
    % Operation Selection Panel
    operationPanel = uipanel(grid, 'Title', '⚙️ Operation Mode', ...
                             'FontWeight', 'bold', 'BackgroundColor', [0.95, 0.98, 1]);
    operationPanel.Layout.Row = 2;
    operationPanel.Layout.Column = 2;
    
    operationGrid = uigridlayout(operationPanel, [2, 1]);
    operationGrid.RowHeight = {'fit', 'fit'};
    operationGrid.RowSpacing = 10;
    
    
    % operationBG = uibuttongroup(operationGrid, 'Title', 'Select Operation');
    % trainRB = uiradiobutton(operationBG, 'Text', 'Train New Models', 'Value', true, 'Position', [10, 35, 150, 22]);
    % predictRB = uiradiobutton(operationBG, 'Text', 'Predict Using Existing Models', 'Position', [10, 10, 150, 22]);
    
    % Data source selection
    uilabel(operationGrid, 'Text', 'Data Source:', 'FontWeight', 'bold');
    % In createDNNPanel function, around line 318:
    dataSourceDD = uidropdown(operationGrid, 'Items', {'Use Current Simulation', 'Load from WAV Files'}, ...
                          'Value', 'Use Current Simulation');
    
    % Training Parameters Panel
    trainPanel = uipanel(grid, 'Title', '📊 Training Parameters', ...
                         'FontWeight', 'bold', 'BackgroundColor', [0.98, 1, 0.98]);
    trainPanel.Layout.Row = 3;
    trainPanel.Layout.Column = 1;
    
    trainGrid = uigridlayout(trainPanel, [6, 2]);
    trainGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    
    uilabel(trainGrid, 'Text', 'Epochs:');
    epochsSpinner = uispinner(trainGrid, 'Value', 100, 'Limits', [10, 1000], 'Step', 10);
    
    uilabel(trainGrid, 'Text', 'Batch Size:');
    batchSizeSpinner = uispinner(trainGrid, 'Value', 32, 'Limits', [4, 256], 'Step', 4);
    
    uilabel(trainGrid, 'Text', 'Learning Rate:');
    learningRateSpinner = uispinner(trainGrid, 'Value', 0.001, 'Limits', [0.0001, 0.1], 'Step', 0.0001, 'ValueDisplayFormat', '%.4f');
    
    uilabel(trainGrid, 'Text', 'Validation Split:');
    valSplitSpinner = uispinner(trainGrid, 'Value', 0.2, 'Limits', [0.1, 0.5], 'Step', 0.05, 'ValueDisplayFormat', '%.2f');
    
    dataAugmentationCB = uicheckbox(trainGrid, 'Text', 'Data Augmentation', 'Value', true);
    dataAugmentationCB.Layout.Column = [1, 2];
    
    gpuTrainingCB = uicheckbox(trainGrid, 'Text', 'GPU Training (if available)', 'Value', true);
    gpuTrainingCB.Layout.Column = [1, 2];
    
    % Model Performance Panel
    perfPanel = uipanel(grid, 'Title', '📈 Model Performance', ...
                        'FontWeight', 'bold', 'BackgroundColor', [1, 0.98, 1]);
    perfPanel.Layout.Row = 3;
    perfPanel.Layout.Column = 2;
    
    perfGrid = uigridlayout(perfPanel, [6, 1]);
    
    uilabel(perfGrid, 'Text', 'Model Status:', 'FontWeight', 'bold');
    modelStatusLabel = uilabel(perfGrid, 'Text', 'No models loaded', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(perfGrid, 'Text', 'Last Training Results:', 'FontWeight', 'bold');
    resultsTextArea = uitextarea(perfGrid, 'Value', {'Click "Check Models" to see available models'}, ...
                                'Editable', 'off', 'FontName', 'Courier');
    resultsTextArea.Layout.Row = [4, 6];
    
    checkModelsBtn = uibutton(perfGrid, 'Text', '🔍 Check Models', ...
                             'BackgroundColor', [0.3, 0.6, 0.9], 'FontColor', 'white');
    
    % Control Buttons Panel
    dnnControlPanel = uipanel(grid, 'BackgroundColor', [0.95, 0.95, 0.95]);
    dnnControlPanel.Layout.Row = 4;
    dnnControlPanel.Layout.Column = [1, 2];
    
    dnnControlGrid = uigridlayout(dnnControlPanel, [2, 4]);
    dnnControlGrid.RowHeight = {'1x', 'fit'};
    dnnControlGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
    
    % Progress and status area
    dnnProgressText = uilabel(dnnControlGrid, 'Text', 'Ready to train/predict', ...
                             'FontColor', [0, 0.6, 0], 'FontWeight', 'bold');
    dnnProgressText.Layout.Row = 1;
    dnnProgressText.Layout.Column = [1, 4];
    
    uilabel(dnnControlGrid, 'Text', ''); % Spacer
    
    trainDNNBtn = uibutton(dnnControlGrid, 'Text', '🎓 Train Models', ...
                          'FontSize', 14, 'FontWeight', 'bold', ...
                          'BackgroundColor', [0.8, 0.3, 0.3], 'FontColor', 'white');
    
    predictDNNBtn = uibutton(dnnControlGrid, 'Text', '🔮 Make Predictions', ...
                            'FontSize', 14, 'FontWeight', 'bold', ...
                            'BackgroundColor', [0.3, 0.7, 0.8], 'FontColor', 'white');
    
    generateDataBtn = uibutton(dnnControlGrid, 'Text', '📊 Generate Dataset', ...
                              'FontSize', 12, ...
                              'BackgroundColor', [0.6, 0.4, 0.8], 'FontColor', 'white');
    
    % Set button callbacks
    trainDNNBtn.ButtonPushedFcn = @(~,~) runDNNTraining();
    predictDNNBtn.ButtonPushedFcn = @(~,~) runDNNPredictionNew();
    generateDataBtn.ButtonPushedFcn = @(~,~) runGenerateDataset();
    checkModelsBtn.ButtonPushedFcn = @(~,~) checkExistingModels();
    
    % Store UI components for access in callbacks
    panel.UserData = struct('azimuthRB', azimuthRB, 'distanceRB', distanceRB, ...
                           'elevationRB', elevationRB, 'combinedCB', combinedCB, ...
                           'dataSource', dataSourceDD, 'epochs', epochsSpinner, ...
                           'batchSize', batchSizeSpinner, 'learningRate', learningRateSpinner, ...
                           'valSplit', valSplitSpinner, 'dataAugmentation', dataAugmentationCB, ...
                           'gpuTraining', gpuTrainingCB, 'modelStatus', modelStatusLabel, ...
                           'resultsText', resultsTextArea, 'progressText', dnnProgressText);
    
    function runDNNTraining()
        % Training logic here
        mainGuiData = ancestor(panel, 'figure').UserData;
        try
            showProgress('DNN Training', 'Preparing training data...');
            
            % Determine which models to train
            params = collectDNNParameters(panel.UserData);
            
            % Update progress
            showProgress('DNN Training', 'Training neural networks...');
            
            % Call the DNN training scripts
            runDNNWorkflow(params);
            
            hideProgress();
            mainGuiData.statusText.Text = '✅ DNN Training Complete!';
            mainGuiData.statusText.FontColor = [0, 0.6, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   'DNN training completed successfully!', 'Training Complete', 'Icon', 'success');
                   
        catch ME
            hideProgress();
            mainGuiData.statusText.Text = '❌ DNN Training Failed';
            mainGuiData.statusText.FontColor = [0.8, 0, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   sprintf('Training failed: %s', ME.message), ...
                   'Training Error', 'Icon', 'error');
        end
        
        function showProgress(title, message)
            mainGuiData.progressBar.Title = title;
            mainGuiData.progressBar.Message = message;
            mainGuiData.progressBar.Visible = 'on';
            drawnow;
        end
        
        function hideProgress()
            mainGuiData.progressBar.Visible = 'off';
            drawnow;
        end
    end

    function runDNNPredictionNew()
        % Get selected models
        params = collectDNNParameters(panel.UserData);
        
        % Check which models are selected
        selected_models = {};
        missing_models = {};
        
       if params.train_azimuth
            if exist('azimuth_model_distance_techniques.mat', 'file') || exist('azimuth_model.mat', 'file')
                selected_models{end+1} = 'Azimuth';
            else
                missing_models{end+1} = 'azimuth_model.mat';
            end
        end
        
        if params.train_distance
            if exist('distance_model.mat', 'file')
                selected_models{end+1} = 'Distance';
            else
                missing_models{end+1} = 'distance_model.mat';
            end
        end
        
        if params.train_elevation
            if exist('elevation_model.mat', 'file')
                selected_models{end+1} = 'Elevation';
            else
                missing_models{end+1} = 'elevation_model.mat';
            end
        end
        
        if params.train_combined
            if exist('combined_3D_model.mat', 'file')
                selected_models{end+1} = 'Combined 3D';
            else
                missing_models{end+1} = 'combined_3D_model.mat';
            end
        end
        
        % Count WAV files
        wav_count = 0;
        if exist('acoustic_wavs', 'dir')
            wav_files = dir('acoustic_wavs/*.wav');
            wav_count = length(wav_files);
        end
        
        % Show configuration dialog
        prediction_config = showPredictionConfigDialog(selected_models, missing_models, wav_count);
        
        if ~isempty(prediction_config)
            % Pass config to params
            params.use_wav_files = prediction_config.use_wav_files;
            params.use_current_sim = prediction_config.use_current_sim;
            
            % Start prediction process
            runPredictionProcess(prediction_config, params);
        end
    end

    function runGenerateDataset()
        % Generate Dataset button callback
        mainGuiData = ancestor(panel, 'figure').UserData;
        try
            showProgress('Dataset Generation', 'Preparing dataset generation...');
            
            % Show dataset configuration dialog
            dataset_config = showDatasetConfigDialog();
            
            if isempty(dataset_config)
                return; % User cancelled
            end
            
            % Collect DNN parameters from GUI
            params = collectDNNParameters(panel.UserData);
            
            % Add dataset configuration
            params.generate_dataset_requested = true;
            params.dataset_config = dataset_config;
            
            % Call the unified DNN workflow
            runDNNWorkflow(params);
            
            hideProgress();
            mainGuiData.statusText.Text = '✅ Dataset Generation Complete!';
            mainGuiData.statusText.FontColor = [0, 0.6, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   'Dataset generation and DNN workflow completed successfully!', ...
                   'Generation Complete', 'Icon', 'success');
                   
        catch ME
            hideProgress();
            mainGuiData.statusText.Text = '❌ Dataset Generation Failed';
            mainGuiData.statusText.FontColor = [0.8, 0, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   sprintf('Dataset generation failed: %s', ME.message), ...
                   'Generation Error', 'Icon', 'error');
        end
        
        function showProgress(title, message)
            mainGuiData.progressBar.Title = title;
            mainGuiData.progressBar.Message = message;
            mainGuiData.progressBar.Visible = 'on';
            drawnow;
        end
        
        function hideProgress()
            mainGuiData.progressBar.Visible = 'off';
            drawnow;
        end
    end
    
    
    
    function checkExistingModels()
        % Check for existing model files
        modelFiles = {'azimuth_model.mat', 'distance_model.mat', 'elevation_model.mat', 'combined_3D_model.mat'};
        status = {};
        
        for i = 1:length(modelFiles)
            if exist(modelFiles{i}, 'file')
                status{end+1} = sprintf('✅ %s found', modelFiles{i});
            else
                status{end+1} = sprintf('❌ %s missing', modelFiles{i});
            end
        end
        
        panel.UserData.resultsText.Value = status;
        
        % Update status label
        foundModels = sum(cellfun(@(x) exist(x, 'file'), modelFiles));
        if foundModels == 0
            panel.UserData.modelStatus.Text = 'No models found';
            panel.UserData.modelStatus.FontColor = [0.8, 0, 0];
        elseif foundModels == length(modelFiles)
            panel.UserData.modelStatus.Text = 'All models available';
            panel.UserData.modelStatus.FontColor = [0, 0.6, 0];
        else
            panel.UserData.modelStatus.Text = sprintf('%d/%d models found', foundModels, length(modelFiles));
            panel.UserData.modelStatus.FontColor = [0.8, 0.6, 0];
        end
    end
end

function config = showDatasetConfigDialog()
    % Create modal dialog
    fig = uifigure('Name', 'Dataset Generation Configuration', ...
                   'Position', [300, 300, 500, 400], ...
                   'WindowStyle', 'modal');
    
    grid = uigridlayout(fig, [4, 1]);
    grid.RowHeight = {150, 120, 80, 50};
    
    % Position Configuration Panel
    posPanel = uipanel(grid, 'Title', 'Position Grid Configuration');
    posGrid = uigridlayout(posPanel, [4, 4]);
    
    uilabel(posGrid, 'Text', 'Min Position [X,Y,Z]:');
    minPosEdit = uieditfield(posGrid, 'Value', '1.5, 1.5, 1.5');
    uilabel(posGrid, 'Text', 'Max Position [X,Y,Z]:');
    maxPosEdit = uieditfield(posGrid, 'Value', '3.5, 2.5, 2.5');
    
    uilabel(posGrid, 'Text', 'Grid Points [X,Y,Z]:');
    gridPointsEdit = uieditfield(posGrid, 'Value', '3, 2, 2');
    uilabel(posGrid, 'Text', 'Total Positions:');
    totalPosLabel = uilabel(posGrid, 'Text', '12', 'FontWeight', 'bold');
    
    % Microphone Selection Panel
    micPanel = uipanel(grid, 'Title', 'Microphone Selection');
    micGrid = uigridlayout(micPanel, [3, 4]);
    
    uilabel(micGrid, 'Text', 'Selection Mode:');
    selectionModeDD = uidropdown(micGrid, 'Items', {'all', 'range', 'custom', 'random'}, 'Value', 'range');
    
    uilabel(micGrid, 'Text', 'Range [start:end]:');
    rangeEdit = uieditfield(micGrid, 'Value', '1:4');
    
    uilabel(micGrid, 'Text', 'Custom Indices:');
    customEdit = uieditfield(micGrid, 'Value', '[1,2,3,4]', 'Enable', 'off');
    
    uilabel(micGrid, 'Text', 'Random Count:');
    randomSpinner = uispinner(micGrid, 'Value', 4, 'Limits', [1, 16], 'Enable', 'off');
    
    % Array Override Panel  
    arrayPanel = uipanel(grid, 'Title', 'Array Override (Optional)');
    arrayGrid = uigridlayout(arrayPanel, [2, 4]);
    
    overrideCB = uicheckbox(arrayGrid, 'Text', 'Override Array Center');
    arrayCenterEdit = uieditfield(arrayGrid, 'Value', '0.5, 0.5, 0.5', 'Enable', 'off');
    
    % Buttons
    buttonPanel = uipanel(grid);
    buttonGrid = uigridlayout(buttonPanel, [1, 3]);
    buttonGrid.ColumnWidth = {'1x', 'fit', 'fit'};
    
    uilabel(buttonGrid, 'Text', '');
    cancelBtn = uibutton(buttonGrid, 'Text', 'Cancel');
    okBtn = uibutton(buttonGrid, 'Text', 'Generate Dataset', 'FontWeight', 'bold');
    
    % ==================== ALL CALLBACKS ==================== 
    
    % Position callbacks - update total when any position field changes
    minPosEdit.ValueChangedFcn = @(~,~) updateTotalPositions();
    maxPosEdit.ValueChangedFcn = @(~,~) updateTotalPositions();
    gridPointsEdit.ValueChangedFcn = @(~,~) updateTotalPositions();
    
    % Microphone selection callback - enable/disable fields based on mode
    selectionModeDD.ValueChangedFcn = @(~,~) updateMicrophoneFields();
    
    % Array override callback - enable/disable array center field
    overrideCB.ValueChangedFcn = @(~,~) toggleArrayCenter();
    
    % Button callbacks
    config = [];
    okBtn.ButtonPushedFcn = @(~,~) acceptConfig();
    cancelBtn.ButtonPushedFcn = @(~,~) close(fig);
    
    % Initialize field states
    updateMicrophoneFields(); % Set initial microphone field states
    
    % Wait for user input
    uiwait(fig);
    
    % ==================== CALLBACK FUNCTIONS ==================== 
    
    function updateTotalPositions()
        try
            gridPoints = str2num(gridPointsEdit.Value);
            if length(gridPoints) == 3 && all(gridPoints > 0)
                total = gridPoints(1) * gridPoints(2) * gridPoints(3);
                totalPosLabel.Text = sprintf('%d', total);
                totalPosLabel.FontColor = [0, 0, 0]; % Black for valid
            else
                totalPosLabel.Text = 'Invalid';
                totalPosLabel.FontColor = [1, 0, 0]; % Red for invalid
            end
        catch
            totalPosLabel.Text = 'Error';
            totalPosLabel.FontColor = [1, 0, 0]; % Red for error
        end
    end
    
    function updateMicrophoneFields()
        mode = selectionModeDD.Value;
        
        % Disable all first
        rangeEdit.Enable = 'off';
        customEdit.Enable = 'off';
        randomSpinner.Enable = 'off';
        
        % Enable based on selection
        switch mode
            case 'range'
                rangeEdit.Enable = 'on';
            case 'custom'
                customEdit.Enable = 'on';
            case 'random'
                randomSpinner.Enable = 'on';
            case 'all'
                % All disabled - use all available mics
        end
    end
    
    function toggleArrayCenter()
        if overrideCB.Value
            arrayCenterEdit.Enable = 'on';
            arrayCenterEdit.BackgroundColor = [1, 1, 1]; % White background
        else
            arrayCenterEdit.Enable = 'off';
            arrayCenterEdit.BackgroundColor = [0.94, 0.94, 0.94]; % Gray background
        end
    end
    
    function acceptConfig()
        config = struct();
        
        % Parse positions
        try
            minPos = str2num(minPosEdit.Value);
            maxPos = str2num(maxPosEdit.Value);
            gridPoints = str2num(gridPointsEdit.Value);
            
            if length(minPos) ~= 3 || length(maxPos) ~= 3 || length(gridPoints) ~= 3
                uialert(fig, 'Position values must be 3-element vectors [X,Y,Z]', 'Invalid Input');
                return;
            end
            
            config.min_position = minPos;
            config.max_position = maxPos;
            config.grid_points = gridPoints;
        catch
            uialert(fig, 'Invalid position format. Use: X, Y, Z', 'Input Error');
            return;
        end
        
        % Mic selection
        config.mic_selection_mode = selectionModeDD.Value;
        config.mic_range = rangeEdit.Value;
        config.mic_custom = customEdit.Value;
        config.mic_random_count = randomSpinner.Value;
        
        % Array override
        config.override_array = overrideCB.Value;
        if overrideCB.Value
            try
                config.array_center = str2num(arrayCenterEdit.Value);
                if length(config.array_center) ~= 3
                    uialert(fig, 'Array center must be 3-element vector [X,Y,Z]', 'Invalid Input');
                    return;
                end
            catch
                uialert(fig, 'Invalid array center format. Use: X, Y, Z', 'Input Error');
                return;
            end
        end
        
        close(fig);
    end
end

%% ============= IR ANALYSIS PANEL =============
function panel = createAnalysisPanel(parent)
    panel = uipanel(parent, 'BackgroundColor', 'white', 'Visible', 'off');
    
    grid = uigridlayout(panel, [4, 2]);
    grid.RowHeight = {120, 150, '1x', 60};
    grid.ColumnWidth = {'1x', '1x'};
    
    % Title and description
    titlePanel = uipanel(grid, 'BackgroundColor', [0.9, 1, 0.9]);
    titlePanel.Layout.Row = 1;
    titlePanel.Layout.Column = [1, 2];
    
    titleGrid = uigridlayout(titlePanel, [3, 1]);
    uilabel(titleGrid, 'Text', '📊 IMPULSE RESPONSE ANALYSIS & COMPARISON', ...
            'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    uilabel(titleGrid, 'Text', 'Compare measured vs simulated impulse responses with advanced acoustic metrics', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    uilabel(titleGrid, 'Text', 'RT60, clarity metrics, cross-correlation, and frequency domain analysis', ...
            'FontSize', 10, 'FontAngle', 'italic', 'HorizontalAlignment', 'center');
    
    % Data Source Panel
    dataPanel = uipanel(grid, 'Title', '📁 Data Sources', ...
                        'FontWeight', 'bold', 'BackgroundColor', [0.98, 1, 0.98]);
    dataPanel.Layout.Row = 2;
    dataPanel.Layout.Column = 1;
    
    dataGrid = uigridlayout(dataPanel, [4, 3]);
    dataGrid.ColumnWidth = {'fit', '1x', 'fit'};
    
    uilabel(dataGrid, 'Text', 'Measured Data:', 'FontWeight', 'bold');
    measuredStatusLabel = uilabel(dataGrid, 'Text', 'No file selected', 'FontColor', [0.8, 0, 0]);
    browseMeasuredBtn = uibutton(dataGrid, 'Text', '📁 Browse', 'BackgroundColor', [0.3, 0.6, 0.9], 'FontColor', 'white');
    
    uilabel(dataGrid, 'Text', 'Simulated Data:', 'FontWeight', 'bold');
    simulatedStatusLabel = uilabel(dataGrid, 'Text', 'No file selected', 'FontColor', [0.8, 0, 0]);
    browseSimulatedBtn = uibutton(dataGrid, 'Text', '📁 Browse', 'BackgroundColor', [0.3, 0.6, 0.9], 'FontColor', 'white');
    
    % Analysis Options Panel
    optionsPanel = uipanel(grid, 'Title', '⚙️ Analysis Options', ...
                           'FontWeight', 'bold', 'BackgroundColor', [0.98, 0.98, 1]);
    optionsPanel.Layout.Row = 2;
    optionsPanel.Layout.Column = 2;
    
    optionsGrid = uigridlayout(optionsPanel, [4, 1]);
    
    enableRT60CB = uicheckbox(optionsGrid, 'Text', 'RT60 Analysis', 'Value', true, 'FontWeight', 'bold');
    enableClarityMS = uicheckbox(optionsGrid, 'Text', 'Clarity Metrics (C50, C80)', 'Value', true, 'FontWeight', 'bold');
    enableFreqAnalysisCB = uicheckbox(optionsGrid, 'Text', 'Frequency Domain Analysis', 'Value', true, 'FontWeight', 'bold');
    enableStatsCB = uicheckbox(optionsGrid, 'Text', 'Statistical Analysis', 'Value', true, 'FontWeight', 'bold');
    
    % Results Display Panel
    resultsDisplayPanel = uipanel(grid, 'Title', '📈 Analysis Results', ...
                                  'FontWeight', 'bold', 'BackgroundColor', [1, 1, 0.98]);
    resultsDisplayPanel.Layout.Row = 3;
    resultsDisplayPanel.Layout.Column = [1, 2];
    
    resultsGrid = uigridlayout(resultsDisplayPanel, [1, 2]);
    resultsGrid.ColumnWidth = {'1x', 300};
    
    % Results text area
    analysisResultsText = uitextarea(resultsGrid, 'Value', {'No analysis performed yet.'}, ...
                                    'Editable', 'off', 'FontName', 'Courier');
    
    % Quick stats panel
    statsPanel = uipanel(resultsGrid, 'Title', 'Quick Statistics', 'BackgroundColor', [0.95, 0.95, 1]);
    statsGrid = uigridlayout(statsPanel, [6, 2]);
    statsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    
    uilabel(statsGrid, 'Text', 'Avg. Correlation:');
    corrLabel = uilabel(statsGrid, 'Text', '-', 'FontWeight', 'bold');
    
    uilabel(statsGrid, 'Text', 'Avg. NRMSE:');
    nrmseLabel = uilabel(statsGrid, 'Text', '-', 'FontWeight', 'bold');
    
    uilabel(statsGrid, 'Text', 'Avg. SNR (dB):');
    snrLabel = uilabel(statsGrid, 'Text', '-', 'FontWeight', 'bold');
    
    uilabel(statsGrid, 'Text', 'Alignment Quality:');
    alignLabel = uilabel(statsGrid, 'Text', '-', 'FontWeight', 'bold');
    
    uilabel(statsGrid, 'Text', 'Overall Quality:');
    qualityLabel = uilabel(statsGrid, 'Text', 'Not analyzed', 'FontWeight', 'bold');
    
    uilabel(statsGrid, 'Text', 'Microphones:');
    micsLabel = uilabel(statsGrid, 'Text', '-', 'FontWeight', 'bold');
    
    % Control Buttons Panel
    analysisControlPanel = uipanel(grid, 'BackgroundColor', [0.95, 0.95, 0.95]);
    analysisControlPanel.Layout.Row = 4;
    analysisControlPanel.Layout.Column = [1, 2];
    
    analysisControlGrid = uigridlayout(analysisControlPanel, [1, 4]);
    analysisControlGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
    
    uilabel(analysisControlGrid, 'Text', ''); % Spacer
    
    checkDataBtn = uibutton(analysisControlGrid, 'Text', '🔍 Check Data', ...
                           'BackgroundColor', [0.4, 0.6, 0.8], 'FontColor', 'white');
    
    runAnalysisBtn = uibutton(analysisControlGrid, 'Text', '📊 Run Analysis', ...
                             'FontSize', 14, 'FontWeight', 'bold', ...
                             'BackgroundColor', [0.2, 0.7, 0.3], 'FontColor', 'white');
    
    viewPlotsBtn = uibutton(analysisControlGrid, 'Text', '📈 View Plots', ...
                           'BackgroundColor', [0.6, 0.3, 0.8], 'FontColor', 'white');
    
    % Set button callbacks
    checkDataBtn.ButtonPushedFcn = @(~,~) checkAvailableData();
    runAnalysisBtn.ButtonPushedFcn = @(~,~) runIRAnalysis();
    viewPlotsBtn.ButtonPushedFcn = @(~,~) viewAnalysisPlots();
    browseMeasuredBtn.ButtonPushedFcn = @(~,~) browseMeasuredFile();
    browseSimulatedBtn.ButtonPushedFcn = @(~,~) browseSimulatedFile();
    
    % Store UI components
    panel.UserData = struct('measuredStatus', measuredStatusLabel, 'simulatedStatus', simulatedStatusLabel, ...
                       'browseMeasuredBtn', browseMeasuredBtn, 'browseSimulatedBtn', browseSimulatedBtn, ...
                       'measuredFile', '', 'simulatedFile', '', ...
                       'enableRT60', enableRT60CB, 'enableClarity', enableClarityMS, ...
                       'enableFreq', enableFreqAnalysisCB, 'enableStats', enableStatsCB, ...
                       'resultsText', analysisResultsText, 'corrLabel', corrLabel, ...
                       'nrmseLabel', nrmseLabel, 'snrLabel', snrLabel, ...
                       'alignLabel', alignLabel, 'qualityLabel', qualityLabel, ...
                       'micsLabel', micsLabel);
    
    function checkAvailableData()
        % Check for analysis results in latest IR comparison folder
        
        % Find latest IR comparison results folder
        result_folders = dir('IR_comparison_results_*');
        if isempty(result_folders)
            uialert(ancestor(panel, 'figure'), ...
                   'No IR comparison results folder found. Run analysis first.', ...
                   'No Results Folder', 'Icon', 'error');
            return;
        end
        
        % Sort by date (newest first)
        [~, idx] = sort([result_folders.datenum], 'descend');
        latest_folder = result_folders(idx(1)).name;
        
        % Check for expected files from IRs_comparison_analysis_tool3
        expected_files = {
            'IR_comparison_results_*.mat',           % Main results file
            'Time_Domain_Comparison_*.png',          % Time domain plots
            'Frequency_Domain_Comparison_*.png',     % Frequency domain plots  
            'Acoustic_Metrics_Comparison_*.png',     % Acoustic metrics plots
            'template_trimmed_recorded.wav',         % Aligned recorded signal
            'template_trimmed_simulated.wav',        % Aligned simulated signal
            'reference_template.wav'                 % Reference template
        };
        
        missing_files = {};
        found_files = {};
        
        for i = 1:length(expected_files)
            file_pattern = expected_files{i};
            files_found = dir(fullfile(latest_folder, file_pattern));
            
            if isempty(files_found)
                missing_files{end+1} = file_pattern;
            else
                found_files{end+1} = file_pattern;
            end
        end
        
        % Update status display
        status_message = sprintf('Latest folder: %s\n✅ Found: %d files\n❌ Missing: %d files', ...
                               latest_folder, length(found_files), length(missing_files));
        
        if isempty(missing_files)
            panel.UserData.resultsText.Value = {
                sprintf('✅ Complete analysis results found in: %s', latest_folder);
                '';
                'All expected files are present:';
                '• Main results (.mat file)';
                '• All comparison plots (.png files)';
                '• Aligned audio signals (.wav files)';
                '';
                'Ready for viewing plots and results!'
            };
            
            uialert(ancestor(panel, 'figure'), ...
                   sprintf('✅ Complete analysis results found!\n\nFolder: %s\nAll %d expected files are present.', ...
                          latest_folder, length(found_files)), ...
                   'Results Check Complete', 'Icon', 'success');
        else
            panel.UserData.resultsText.Value = {
                sprintf('⚠ Incomplete results in: %s', latest_folder);
                '';
                sprintf('Found %d files, missing %d files:', length(found_files), length(missing_files));
                sprintf('Missing: %s', strjoin(missing_files, ', '));
                '';
                'Consider re-running the analysis.'
            };
            
            uialert(ancestor(panel, 'figure'), ...
                   sprintf('⚠ Incomplete analysis results\n\nFolder: %s\nMissing files:\n%s', ...
                          latest_folder, strjoin(missing_files, '\n• ')), ...
                   'Incomplete Results', 'Icon', 'warning');
        end
        
        % Store latest folder for other functions to use
        panel.UserData.latestResultsFolder = latest_folder;
    end
    
    function runIRAnalysis()
        % Run the impulse response comparison analysis
        mainGuiData = ancestor(panel, 'figure').UserData;
        try
            % Check if files are selected
            if isempty(panel.UserData.measuredFile) || isempty(panel.UserData.simulatedFile)
                hideProgress();
                uialert(ancestor(panel, 'figure'), ...
                       'Please select both measured and simulated WAV files using the Browse buttons.', ...
                       'Files Required', 'Icon', 'warning');
                return;
            end
            
            % Set the file paths for IRs_comparison_analysis_tool3 to use
            assignin('base', 'selected_measured_file', panel.UserData.measuredFile);
            assignin('base', 'selected_simulated_file', panel.UserData.simulatedFile);
            
            showProgress('IR Analysis', 'Loading selected WAV files...');
            
            % Update progress
            showProgress('IR Analysis', 'Running comprehensive analysis...');
            
            % Call the IR comparison analysis tool
            runExternalIRAnalysisComplete();
            
            % Update progress
            showProgress('IR Analysis', 'Processing results...');
            
            % Load and display results
            loadAnalysisResultsComplete(panel);
            
            hideProgress();
            mainGuiData.statusText.Text = '✅ IR Analysis Complete!';
            mainGuiData.statusText.FontColor = [0, 0.6, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   'Impulse response analysis completed successfully!', ...
                   'Analysis Complete', 'Icon', 'success');
                   
        catch ME
            hideProgress();
            mainGuiData.statusText.Text = '❌ IR Analysis Failed';
            mainGuiData.statusText.FontColor = [0.8, 0, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   sprintf('Analysis failed: %s', ME.message), ...
                   'Analysis Error', 'Icon', 'error');
        end
        
        function showProgress(title, message)
            mainGuiData.progressBar.Title = title;
            mainGuiData.progressBar.Message = message;
            mainGuiData.progressBar.Visible = 'on';
            drawnow;
        end
        
        function hideProgress()
            mainGuiData.progressBar.Visible = 'off';
            drawnow;
        end
    end
    
    function viewAnalysisPlots()
        % View plots based on selected analysis options
        
        % Check if we have a results folder
        if ~isfield(panel.UserData, 'latestResultsFolder') || isempty(panel.UserData.latestResultsFolder)
            % Try to find latest folder
            result_folders = dir('IR_comparison_results_*');
            if isempty(result_folders)
                uialert(ancestor(panel, 'figure'), ...
                       'No analysis results found. Run analysis first.', ...
                       'No Results', 'Icon', 'error');
                return;
            end
            [~, idx] = sort([result_folders.datenum], 'descend');
            panel.UserData.latestResultsFolder = result_folders(idx(1)).name;
        end
        
        results_folder = panel.UserData.latestResultsFolder;
        
        % Get user's analysis option preferences
        show_time_domain = panel.UserData.enableRT60.Value;  % Using RT60 checkbox as proxy for time domain
        show_frequency = panel.UserData.enableFreq.Value;
        show_metrics = panel.UserData.enableClarity.Value;  % Using clarity as proxy for metrics
        show_stats = panel.UserData.enableStats.Value;
        
        plots_opened = 0;
        
        % Open Time Domain plots
        if show_time_domain
            time_plots = dir(fullfile(results_folder, 'Time_Domain_Comparison_*.png'));
            for i = 1:length(time_plots)
                plot_path = fullfile(results_folder, time_plots(i).name);
                if exist(plot_path, 'file')
                    img = imread(plot_path);
                    figure('Name', 'Time Domain Analysis', 'Position', [100, 100, 1200, 600]);
                    imshow(img);
                    plots_opened = plots_opened + 1;
                end
            end
        end
        
        % Open Frequency Domain plots  
        if show_frequency
            freq_plots = dir(fullfile(results_folder, 'Frequency_Domain_Comparison_*.png'));
            for i = 1:length(freq_plots)
                plot_path = fullfile(results_folder, freq_plots(i).name);
                if exist(plot_path, 'file')
                    img = imread(plot_path);
                    figure('Name', 'Frequency Domain Analysis', 'Position', [150, 150, 1200, 600]);
                    imshow(img);
                    plots_opened = plots_opened + 1;
                end
            end
        end
        
        % Open Acoustic Metrics plots
        if show_metrics
            metrics_plots = dir(fullfile(results_folder, 'Acoustic_Metrics_Comparison_*.png'));
            for i = 1:length(metrics_plots)
                plot_path = fullfile(results_folder, metrics_plots(i).name);
                if exist(plot_path, 'file')
                    img = imread(plot_path);
                    figure('Name', 'Acoustic Metrics Analysis', 'Position', [200, 200, 1000, 600]);
                    imshow(img);
                    plots_opened = plots_opened + 1;
                end
            end
        end
        
        % Load and display statistical results if requested
        if show_stats
            results_file = dir(fullfile(results_folder, 'IR_comparison_results_*.mat'));
            if ~isempty(results_file)
                try
                    data = load(fullfile(results_folder, results_file(1).name));
                    if isfield(data, 'comparison_results')
                        results = data.comparison_results;
                        
                        % Display statistical summary
                        stats_text = {
                            '=== STATISTICAL ANALYSIS SUMMARY ===';
                            '';
                            sprintf('Correlation: %.3f', results.statistical_analysis.correlation_coefficient);
                            sprintf('NRMSE: %.4f', results.statistical_analysis.nrmse);
                            sprintf('SNR: %.1f dB', results.statistical_analysis.snr_db);
                            sprintf('RMSE: %.4f', results.statistical_analysis.rmse);
                            '';
                            sprintf('Analysis date: %s', char(results.metadata.analysis_timestamp));
                            sprintf('Results folder: %s', results_folder);
                        };
                        
                        panel.UserData.resultsText.Value = stats_text;
                        plots_opened = plots_opened + 1;
                    end
                catch ME
                    fprintf('Warning: Could not load statistical results: %s\n', ME.message);
                end
            end
        end
        
        % Show summary message
        if plots_opened > 0
            uialert(ancestor(panel, 'figure'), ...
                   sprintf('Opened %d plot(s) based on your analysis options.\n\nResults from: %s', ...
                          plots_opened, results_folder), ...
                   'Plots Opened', 'Icon', 'success');
        else
            uialert(ancestor(panel, 'figure'), ...
                   'No plots opened. Please select analysis options (RT60, Frequency, Clarity, Stats) and ensure analysis has been run.', ...
                   'No Plots', 'Icon', 'warning');
        end
    end

    function browseMeasuredFile()
            [file, path] = uigetfile('*.wav', 'Select Measured WAV File');
            if file ~= 0
                panel.UserData.measuredFile = fullfile(path, file);
                panel.UserData.measuredStatus.Text = ['✅ ' file];
                panel.UserData.measuredStatus.FontColor = [0, 0.6, 0];
            end
        end
        
        function browseSimulatedFile()
            [file, path] = uigetfile('*.wav', 'Select Simulated WAV File');
            if file ~= 0
                panel.UserData.simulatedFile = fullfile(path, file);
                panel.UserData.simulatedStatus.Text = ['✅ ' file];
                panel.UserData.simulatedStatus.FontColor = [0, 0.6, 0];
            end
        end


    
end

%% ============= BEAMFORMING PANEL =============
function panel = createBeamformingPanel(parent)
    panel = uipanel(parent, 'BackgroundColor', 'white', 'Visible', 'off');
    
    grid = uigridlayout(panel, [4, 2]);
    grid.RowHeight = {120, 200, '1x', 60};
    grid.ColumnWidth = {'1x', '1x'};
    
    % Title and description
    titlePanel = uipanel(grid, 'BackgroundColor', [0.95, 0.9, 1]);
    titlePanel.Layout.Row = 1;
    titlePanel.Layout.Column = [1, 2];
    
    titleGrid = uigridlayout(titlePanel, [3, 1]);
    uilabel(titleGrid, 'Text', '📡 ADVANCED BEAMFORMING & SOURCE IMAGING', ...
            'FontSize', 18, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    uilabel(titleGrid, 'Text', 'SC-DAMAS, OMP-DAMAS, and conventional beamforming for reverberant environments', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    uilabel(titleGrid, 'Text', 'Real-time source localization and acoustic imaging', ...
            'FontSize', 10, 'FontAngle', 'italic', 'HorizontalAlignment', 'center');
    
    % Algorithm Selection Panel
    algoPanel = uipanel(grid, 'Title', '🎯 Algorithm Selection', ...
                        'FontWeight', 'bold', 'BackgroundColor', [1, 0.98, 1]);
    algoPanel.Layout.Row = 2;
    algoPanel.Layout.Column = 1;
    
    algoGrid = uigridlayout(algoPanel, [5, 1]);
    
    cbfCB = uicheckbox(algoGrid, 'Text', 'Conventional Beamforming (CBF)', 'Value', true, 'FontWeight', 'bold');
    ompDamasCB = uicheckbox(algoGrid, 'Text', 'OMP-DAMAS', 'Value', true, 'FontWeight', 'bold');
    scDamasCB = uicheckbox(algoGrid, 'Text', 'SC-DAMAS (Advanced)', 'Value', true, 'FontWeight', 'bold');
    
    uilabel(algoGrid, 'Text', 'Comparison Mode:', 'FontWeight', 'bold');
    comparisonCB = uicheckbox(algoGrid, 'Text', 'Generate Algorithm Comparison', 'Value', true);
    
    % Processing Parameters Panel
    paramPanel = uipanel(grid, 'Title', '⚙️ Processing Parameters', ...
                         'FontWeight', 'bold', 'BackgroundColor', [0.98, 1, 1]);
    paramPanel.Layout.Row = 2;
    paramPanel.Layout.Column = 2;
    
    paramGrid = uigridlayout(paramPanel, [6, 2]);
    paramGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
    paramGrid.RowSpacing = 3;  % Add this line to reduce spacing

    uilabel(paramGrid, 'Text', 'Frequency (Hz):');
    freqSpinner = uispinner(paramGrid, 'Value', 1000, 'Limits', [100, 8000], 'Step', 100);
    
    uilabel(paramGrid, 'Text', 'Grid Size:');
    gridSizeSpinner = uispinner(paramGrid, 'Value', 15, 'Limits', [5, 51], 'Step', 2);
    
    uilabel(paramGrid, 'Text', 'Scan Area (m):');
    scanAreaSpinner = uispinner(paramGrid, 'Value', 1, 'Limits', [0.5, 5], 'Step', 0.1);
    
    uilabel(paramGrid, 'Text', 'Max Iterations:');
    maxIterSpinner = uispinner(paramGrid, 'Value', 50, 'Limits', [10, 200], 'Step', 10);
    
    uilabel(paramGrid, 'Text', 'Convergence:');
    convergenceSpinner = uispinner(paramGrid, 'Value', 1e-6, 'Limits', [1e-8, 1e-3], 'Step', 1e-7, 'ValueDisplayFormat', '%.0e');
    
    uilabel(paramGrid, 'Text', 'Sparsity Step:');
    sparsitySpinner = uispinner(paramGrid, 'Value', 3, 'Limits', [1, 10], 'Step', 1);
    
    % Visualization Panel
    vizPanel = uipanel(grid, 'Title', '📊 Visualization Options', ...
                       'FontWeight', 'bold', 'BackgroundColor', [1, 1, 0.98]);
    vizPanel.Layout.Row = 3;
    vizPanel.Layout.Column = 1;
    
    vizGrid = uigridlayout(vizPanel, [6, 1]);
    
    plot2DCB = uicheckbox(vizGrid, 'Text', '2D Heat Maps', 'Value', true, 'FontWeight', 'bold');
    plot3DCB = uicheckbox(vizGrid, 'Text', '3D Surface Plots', 'Value', true, 'FontWeight', 'bold');
    plotComparisonCB = uicheckbox(vizGrid, 'Text', 'Side-by-side Comparison', 'Value', true, 'FontWeight', 'bold');
    saveResultsCB = uicheckbox(vizGrid, 'Text', 'Save Results to File', 'Value', true);
    
    uilabel(vizGrid, 'Text', 'Export Format:', 'FontWeight', 'bold');
    exportFormatDD = uidropdown(vizGrid, 'Items', {'MATLAB (.mat)', 'PNG Images', 'Excel (.xlsx)', 'All Formats'}, 'Value', 'All Formats');
    
    % Results Display Panel
    beamResultsPanel = uipanel(grid, 'Title', '📈 Localization Results', ...
                              'FontWeight', 'bold', 'BackgroundColor', [0.98, 0.98, 1]);
    beamResultsPanel.Layout.Row = 3;
    beamResultsPanel.Layout.Column = 2;
    
    beamResultsGrid = uigridlayout(beamResultsPanel, [1, 2]);
    beamResultsGrid.ColumnWidth = {'1x', 150};
    
    % Results table
    beamResultsTable = uitable(beamResultsGrid, 'Data', {}, ...
                              'ColumnName', {'Algorithm', 'Est. X', 'Est. Z', 'Error (m)'}, ...
                              'ColumnEditable', false);
    
    % Quick stats
    beamStatsPanel = uipanel(beamResultsGrid, 'Title', 'Statistics', 'BackgroundColor', [0.95, 0.95, 1]);
    beamStatsGrid = uigridlayout(beamStatsPanel, [5, 1]);
    beamStatsGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', '1x'};
    
    uilabel(beamStatsGrid, 'Text', 'Source Position:', 'FontWeight', 'bold');
    sourcePositionLabel = uilabel(beamStatsGrid, 'Text', 'Not set', 'FontSize', 10);
    
    uilabel(beamStatsGrid, 'Text', 'Best Algorithm:', 'FontWeight', 'bold');
    bestAlgoLabel = uilabel(beamStatsGrid, 'Text', '-', 'FontSize', 10, 'FontColor', [0, 0.6, 0]);
    
    uilabel(beamStatsGrid, 'Text', 'Min Error:', 'FontWeight', 'bold');
    minErrorLabel = uilabel(beamStatsGrid, 'Text', '-', 'FontSize', 10, 'FontColor', [0, 0.6, 0]);
    
    % Control Buttons Panel
    beamControlPanel = uipanel(grid, 'BackgroundColor', [0.95, 0.95, 0.95]);
    beamControlPanel.Layout.Row = 4;
    beamControlPanel.Layout.Column = [1, 2];
    
    beamControlGrid = uigridlayout(beamControlPanel, [2, 4]);
    beamControlGrid.RowHeight = {'1x', 'fit'};
    beamControlGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
    
    % Progress and status area
    beamProgressText = uilabel(beamControlGrid, 'Text', 'Ready for beamforming analysis', ...
                              'FontColor', [0, 0.6, 0], 'FontWeight', 'bold');
    beamProgressText.Layout.Row = 1;
    beamProgressText.Layout.Column = [1, 4];
    
    uilabel(beamControlGrid, 'Text', ''); % Spacer
    
    checkSignalsBtn = uibutton(beamControlGrid, 'Text', '🔍 Check Signals', ...
                              'BackgroundColor', [0.4, 0.6, 0.8], 'FontColor', 'white');
    
    runBeamformingBtn = uibutton(beamControlGrid, 'Text', '📡 Run Beamforming', ...
                                'FontSize', 14, 'FontWeight', 'bold', ...
                                'BackgroundColor', [0.6, 0.2, 0.8], 'FontColor', 'white');
    
    viewBeamPlotsBtn = uibutton(beamControlGrid, 'Text', '📊 View Results', ...
                               'BackgroundColor', [0.3, 0.7, 0.5], 'FontColor', 'white');
    
    % Set button callbacks
    checkSignalsBtn.ButtonPushedFcn = @(~,~) checkBeamformingSignals();
    runBeamformingBtn.ButtonPushedFcn = @(~,~) runBeamformingAnalysis();
    viewBeamPlotsBtn.ButtonPushedFcn = @(~,~) viewBeamformingPlots();
    
    % Store UI components
    panel.UserData = struct('cbfCB', cbfCB, 'ompDamasCB', ompDamasCB, 'scDamasCB', scDamasCB, ...
                           'comparisonCB', comparisonCB, 'freqSpinner', freqSpinner, ...
                           'gridSizeSpinner', gridSizeSpinner, 'scanAreaSpinner', scanAreaSpinner, ...
                           'maxIterSpinner', maxIterSpinner, 'convergenceSpinner', convergenceSpinner, ...
                           'sparsitySpinner', sparsitySpinner, 'plot2DCB', plot2DCB, 'plot3DCB', plot3DCB, ...
                           'plotComparisonCB', plotComparisonCB, 'saveResultsCB', saveResultsCB, ...
                           'exportFormatDD', exportFormatDD, 'resultsTable', beamResultsTable, ...
                           'sourcePositionLabel', sourcePositionLabel, 'bestAlgoLabel', bestAlgoLabel, ...
                           'minErrorLabel', minErrorLabel, 'progressText', beamProgressText);
    
    function checkBeamformingSignals()
        % Check if simulation signals are available for beamforming
        if evalin('base', 'exist(''final_signals'', ''var'')')
            panel.UserData.progressText.Text = '✅ Signals available for beamforming';
            panel.UserData.progressText.FontColor = [0, 0.6, 0];
            
            % Try to get source position
            try
                sourcePos = evalin('base', 'sourcePos');
                panel.UserData.sourcePositionLabel.Text = sprintf('[%.2f, %.2f, %.2f]', sourcePos(1), sourcePos(2), sourcePos(3));
            catch
                panel.UserData.sourcePositionLabel.Text = 'Unknown';
            end
            
            uialert(ancestor(panel, 'figure'), ...
                   'Simulation signals found and ready for beamforming analysis!', ...
                   'Signals Check', 'Icon', 'success');
        else
            panel.UserData.progressText.Text = '❌ No simulation signals found';
            panel.UserData.progressText.FontColor = [0.8, 0, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   'No simulation signals found. Please run the acoustic simulation first.', ...
                   'No Signals', 'Icon', 'error');
        end
    end
    
    function runBeamformingAnalysis()
        % Run the beamforming analysis
        mainGuiData = ancestor(panel, 'figure').UserData;
        try
            showProgress('Beamforming Analysis', 'Checking simulation data...');
            
            % Check if signals are available
            if ~evalin('base', 'exist(''final_signals'', ''var'')')
                hideProgress();
                uialert(ancestor(panel, 'figure'), ...
                       'No simulation signals found. Run acoustic simulation first.', ...
                       'No Data', 'Icon', 'error');
                return;
            end
            
            % Collect parameters
            params = collectBeamformingParameters(panel.UserData);
            
            % Update progress
            showProgress('Beamforming Analysis', 'Running beamforming algorithms...');
            
            % Call the beamforming analysis
            runExternalBeamformingComplete(params, panel);
            
            % Update progress
            showProgress('Beamforming Analysis', 'Processing results...');
            
            % Load and display results
            loadBeamformingResultsComplete(panel);
            
            hideProgress();
            mainGuiData.statusText.Text = '✅ Beamforming Complete!';
            mainGuiData.statusText.FontColor = [0, 0.6, 0];
            
            uialert(ancestor(panel, 'figure'), ...
                   'Beamforming analysis completed successfully!', ...
                   'Analysis Complete', 'Icon', 'success');
                   
        catch ME
            hideProgress();
            mainGuiData.statusText.Text = '❌ Beamforming Failed';
            mainGuiData.statusText.FontColor = [0.8, 0, 0];
            
            % DETAILED ERROR DEBUGGING
            fprintf('🐛 BEAMFORMING ERROR DEBUG:\n');
            fprintf('   Error: %s\n', ME.message);
            fprintf('   Identifier: %s\n', ME.identifier);
            fprintf('   Stack trace:\n');
            for i = 1:length(ME.stack)
                fprintf('     %s (line %d) in %s\n', ME.stack(i).name, ME.stack(i).line, ME.stack(i).file);
            end
            
            uialert(ancestor(panel, 'figure'), ...
                   sprintf('Beamforming failed: %s\n\nCheck Command Window for detailed error trace.', ME.message), ...
                   'Beamforming Error', 'Icon', 'error');
        end
        
        function showProgress(title, message)
            mainGuiData.progressBar.Title = title;
            mainGuiData.progressBar.Message = message;
            mainGuiData.progressBar.Visible = 'on';
            drawnow;
        end
        
        function hideProgress()
            mainGuiData.progressBar.Visible = 'off';
            drawnow;
        end
    end
    
    %% NEW - Controlled visualization function
    function viewBeamformingPlots()
        % Get GUI parameters for what to show
        params = collectBeamformingParameters(panel.UserData);
        
        % Check for results folder
        if evalin('base', 'exist(''latest_beamforming_folder'', ''var'')')
            results_folder = evalin('base', 'latest_beamforming_folder');
            plots_opened = 0;
            
            % Look for PNG files in figures subdirectory
            figures_dir = fullfile(results_folder, 'figures');
            
            if exist(figures_dir, 'dir')
                % Show 2D heat maps if requested
                if params.plot_2d
                    heat_map_files = dir(fullfile(figures_dir, 'beamforming_2D_*.png'));
                    for i = 1:length(heat_map_files)
                        img_path = fullfile(figures_dir, heat_map_files(i).name);
                        img = imread(img_path);
                        figure('Name', sprintf('2D Beamforming - %s', heat_map_files(i).name));
                        imshow(img);
                        plots_opened = plots_opened + 1;
                    end
                end
                
                % Show 3D surface plots if requested  
                if params.plot_3d
                    surface_files = dir(fullfile(figures_dir, 'beamforming_3D_*.png'));
                    for i = 1:length(surface_files)
                        img_path = fullfile(figures_dir, surface_files(i).name);
                        img = imread(img_path);
                        figure('Name', sprintf('3D Beamforming - %s', surface_files(i).name));
                        imshow(img);
                        plots_opened = plots_opened + 1;
                    end
                end
                
                if plots_opened > 0
                    fprintf('✅ Opened %d beamforming plot(s)\n', plots_opened);
                else
                    fprintf('⚠ No plots to show (check visualization options)\n');
                    uialert(ancestor(panel, 'figure'), ...
                           'No plots selected. Enable 2D or 3D visualization options.', ...
                           'No Plots', 'Icon', 'warning');
                end
            else
                fprintf('⚠ Figures directory not found: %s\n', figures_dir);
                uialert(ancestor(panel, 'figure'), ...
                       sprintf('Figures directory not found:\n%s', figures_dir), ...
                       'No Figures', 'Icon', 'error');
            end
        else
            uialert(ancestor(panel, 'figure'), ...
                   'No beamforming results found. Run analysis first.', ...
                   'No Results', 'Icon', 'error');
        end
    end

end

%% ============= HELPER FUNCTIONS FOR EXTERNAL MODULE CALLS =============

function params = collectSimulationParameters(uiData)
    % Room configuration
    params.room_type = lower(uiData.roomType.Value);
    params.room_length = uiData.roomLength.Value;
    params.room_width = uiData.roomWidth.Value; 
    params.room_height = uiData.roomHeight.Value;
    params.temperature = uiData.temperature.Value;
    params.humidity = uiData.humidity.Value;
    
    % Array configuration - PARSE ARRAY CENTER
    params.array_type = lower(strrep(uiData.arrayType.Value, ' ', '_'));
    params.numMics = uiData.numMics.Value;
    params.array_orientation = lower(uiData.orientation.Value);
    
    % Parse array center from text field
    arrayCenterStr = uiData.arrayCenter.Value;
    fprintf('Debug: arrayCenterStr = "%s"\n', arrayCenterStr);
    try
        centerCoords = str2num(arrayCenterStr); %#ok<ST2NM>
        if length(centerCoords) == 3
            params.array_center_x = centerCoords(1);
            params.array_center_y = centerCoords(2);
            params.array_center_z = centerCoords(3);

            fprintf('Debug: Set params.array_center_x = %g\n', params.array_center_x);
            fprintf('Debug: Set params.array_center_y = %g\n', params.array_center_y);
            fprintf('Debug: Set params.array_center_z = %g\n', params.array_center_z);
        else
            error('Invalid array center format');
        end
    catch
        warning('Could not parse array center "%s". Using default [0.5, 0.5, 0.5].', arrayCenterStr);
        params.array_center_x = 0.5;
        params.array_center_y = 0.5;
        params.array_center_z = 0.5;
    end

    % Source configuration - PARSE POSITIONS
    startPosStr = uiData.startPosition.Value;
    try
        startCoords = str2num(startPosStr); %#ok<ST2NM>
        if length(startCoords) == 3
            params.start_pos = startCoords;
        else
            error('Invalid start position format');
        end
    catch
        warning('Could not parse start position "%s". Using default [2, 2, 2].', startPosStr);
        params.start_pos = [2, 2, 2];
    end
    
    endPosStr = uiData.endPosition.Value;
    try
        endCoords = str2num(endPosStr); %#ok<ST2NM>
        if length(endCoords) == 3
            params.end_pos = endCoords;
        else
            error('Invalid end position format');
        end
    catch
        warning('Could not parse end position "%s". Using start position.', endPosStr);
        params.end_pos = params.start_pos;
    end
    
    % Rest of parameters...
    params.signal_choice = lower(strrep(uiData.signalType.Value, ' ', '_'));
    params.duration = uiData.duration.Value;
    params.source_power = uiData.power.Value;
    params.enable_noise = uiData.enableNoise.Value;
    params.noise_type = lower(uiData.noiseType.Value);
    params.max_reflection_order = uiData.maxRefl.Value;
    params.fs = uiData.fs.Value;
    params.enable_visualization = uiData.enableVis.Value;
    params.use_freq_domain_reflections = uiData.freqDomainRefl.Value;
    params.use_vehicle_body = uiData.vehicleBody.Value;
end



function params = collectDNNParameters(uiData)
    % Collect DNN training/prediction parameters
    params = struct();
    
    % Model selection
    params.train_azimuth = uiData.azimuthRB.Value;
    params.train_distance = uiData.distanceRB.Value;
    params.train_elevation = uiData.elevationRB.Value;
    params.train_combined = uiData.combinedCB.Value;
    
    % Operation mode
 
    params.training_mode = false; % Always prediction mode now
    params.data_source = uiData.dataSource.Value;
    
    % Training parameters
    params.epochs = uiData.epochs.Value;
    params.batch_size = uiData.batchSize.Value;
    params.learning_rate = uiData.learningRate.Value;
    params.validation_split = uiData.valSplit.Value;
    params.data_augmentation = uiData.dataAugmentation.Value;
    params.gpu_training = uiData.gpuTraining.Value;
end


function params = collectBeamformingParameters(uiData)
    % Collect beamforming analysis parameters
    params = struct();
    
    % Algorithm selection
    params.run_cbf = uiData.cbfCB.Value;
    params.run_omp_damas = uiData.ompDamasCB.Value;
    params.run_sc_damas = uiData.scDamasCB.Value;
    params.comparison_mode = uiData.comparisonCB.Value;
    
    % Processing parameters
    params.frequency = uiData.freqSpinner.Value;
    params.grid_size = uiData.gridSizeSpinner.Value;
    params.scan_area_size = uiData.scanAreaSpinner.Value;
    params.max_iterations = uiData.maxIterSpinner.Value;
    params.delta_criteria = uiData.convergenceSpinner.Value;
    params.sparsity_step = uiData.sparsitySpinner.Value;
    
    % Visualization options
    params.plot_2d = uiData.plot2DCB.Value;
    params.plot_3d = uiData.plot3DCB.Value;
    params.plot_comparison = uiData.plotComparisonCB.Value;
    params.save_results = uiData.saveResultsCB.Value;
    params.export_format = uiData.exportFormatDD.Value;
end



function runExternalMeasurement()
    % Run the measurement and recording system
    evalin('base', 'run(''recording_and_analysis_session.m'')');
end

%% ============= TAB CREATION FUNCTIONS FOR RESULTS PANEL =============

function createSimulationResultsTab(parent)
    % Create simulation results visualization
    grid = uigridlayout(parent, [2, 2]);
    grid.RowHeight = {'1x', 'fit'};
    grid.ColumnWidth = {'1x', 300};
    
    % Main plot area
    plotPanel = uipanel(grid, 'Title', 'Simulation Results', 'BackgroundColor', 'white');
    plotPanel.Layout.Row = 1;
    plotPanel.Layout.Column = 1;
    
    % Results summary
    summaryPanel = uipanel(grid, 'Title', 'Summary', 'BackgroundColor', [0.98, 0.98, 1]);
    summaryPanel.Layout.Row = 1;
    summaryPanel.Layout.Column = 2;
    
    summaryGrid = uigridlayout(summaryPanel, [8, 1]);
    uilabel(summaryGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
    
    % Control buttons
    controlPanel = uipanel(grid, 'BackgroundColor', [0.95, 0.95, 0.95]);
    controlPanel.Layout.Row = 2;
    controlPanel.Layout.Column = [1, 2];
    
    controlGrid = uigridlayout(controlPanel, [1, 4]);
    controlGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
    
    uilabel(controlGrid, 'Text', ''); % Spacer
    
    loadResultsBtn = uibutton(controlGrid, 'Text', '📊 Load Results', ...
                         'BackgroundColor', [0.3, 0.6, 0.9], 'FontColor', 'white', ...
                         'ButtonPushedFcn', @(~,~) loadSimulationResults());
    uibutton(controlGrid, 'Text', '🔊 Play Audio', 'BackgroundColor', [0.2, 0.7, 0.3], 'FontColor', 'white');
    uibutton(controlGrid, 'Text', '💾 Export Data', 'BackgroundColor', [0.6, 0.3, 0.8], 'FontColor', 'white');
end

function createDNNResultsTab(parent)
    % Create DNN results visualization
    grid = uigridlayout(parent, [2, 2]);
    grid.RowHeight = {'1x', 'fit'};
    grid.ColumnWidth = {'1x', 300};
    
    % Training/validation plots
    plotPanel = uipanel(grid, 'Title', 'Training Progress & Validation', 'BackgroundColor', 'white');
    plotPanel.Layout.Row = 1;
    plotPanel.Layout.Column = 1;
    
    % Model performance summary
    perfPanel = uipanel(grid, 'Title', 'Model Performance', 'BackgroundColor', [0.98, 1, 0.98]);
    perfPanel.Layout.Row = 1;
    perfPanel.Layout.Column = 2;
    
    perfGrid = uigridlayout(perfPanel, [10, 2]);
    perfGrid.RowHeight = repmat({'fit'}, 1, 10);
    
    uilabel(perfGrid, 'Text', 'Azimuth Model:', 'FontWeight', 'bold');
    uilabel(perfGrid, 'Text', 'Not trained', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(perfGrid, 'Text', 'RMSE:');
    uilabel(perfGrid, 'Text', '- deg', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(perfGrid, 'Text', 'Distance Model:', 'FontWeight', 'bold');
    uilabel(perfGrid, 'Text', 'Not trained', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(perfGrid, 'Text', 'RMSE:');
    uilabel(perfGrid, 'Text', '- m', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(perfGrid, 'Text', 'Elevation Model:', 'FontWeight', 'bold');
    uilabel(perfGrid, 'Text', 'Not trained', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(perfGrid, 'Text', 'RMSE:');
    uilabel(perfGrid, 'Text', '- deg', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(perfGrid, 'Text', 'Combined 3D:', 'FontWeight', 'bold');
    uilabel(perfGrid, 'Text', 'Not trained', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(perfGrid, 'Text', 'Avg Error:');
    uilabel(perfGrid, 'Text', '- m', 'FontColor', [0.6, 0.6, 0.6]);
    
    % Control buttons
    controlPanel = uipanel(grid, 'BackgroundColor', [0.95, 0.95, 0.95]);
    controlPanel.Layout.Row = 2;
    controlPanel.Layout.Column = [1, 2];
    
    controlGrid = uigridlayout(controlPanel, [1, 4]);
    controlGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
    
    uilabel(controlGrid, 'Text', ''); % Spacer
    
    uibutton(controlGrid, 'Text', '📈 Load Models', 'BackgroundColor', [0.3, 0.6, 0.9], 'FontColor', 'white');
    uibutton(controlGrid, 'Text', '🎯 Test Prediction', 'BackgroundColor', [0.8, 0.3, 0.3], 'FontColor', 'white');
    uibutton(controlGrid, 'Text', '📊 Compare Models', 'BackgroundColor', [0.6, 0.3, 0.8], 'FontColor', 'white');
end

function loadSimulationResults()
    % Load and display simulation results
    try
        fprintf('Loading simulation results...\n');
        
        % Check if simulation data exists in workspace
        if evalin('base', 'exist(''simulation_results'', ''var'')')
            results = evalin('base', 'simulation_results');
            
            % Display results in command window (always works)
            fprintf('\n=== SIMULATION RESULTS LOADED ===\n');
            fprintf('Execution Time: %.2f seconds\n', results.execution_time);
            fprintf('Microphones: %d\n', results.data.array.num_mics);
            fprintf('Sampling Rate: %d Hz\n', results.data.simulation.sampling_rate);
            fprintf('Room Type: %s\n', results.data.room.type);
            fprintf('Array Type: %s\n', results.data.array.type);
            fprintf('Success: %s\n', string(results.success));
            fprintf('Output Directory: %s\n', results.directories.base_dir);
            fprintf('\nWorkspace Variables Created:\n');
            fprintf('• simulation_results (complete structure)\n');
            fprintf('• final_signals (microphone signals)\n');
            fprintf('• micPositions (array geometry)\n');
            fprintf('• fs (sampling rate)\n');
            fprintf('• sourcePos (source position)\n');
            fprintf('================================\n\n');
            
            % Try to show dialog, fallback to command window if it fails
            try
                % Find the main UI figure more reliably
                allFigs = findall(groot, 'Type', 'figure');
                uiFig = [];
                for i = 1:length(allFigs)
                    if contains(allFigs(i).Name, 'Advanced Acoustic Analysis')
                        uiFig = allFigs(i);
                        break;
                    end
                end
                
                if ~isempty(uiFig)
                    message = sprintf(['Simulation Results Loaded Successfully!\n\n' ...
                                      'Execution Time: %.2f seconds\n' ...
                                      'Microphones: %d\n' ...
                                      'Sampling Rate: %d Hz\n' ...
                                      'Room Type: %s\n' ...
                                      'Array Type: %s\n\n' ...
                                      'Check Command Window for details.\n' ...
                                      'Output Directory:\n%s'], ...
                                      results.execution_time, ...
                                      results.data.array.num_mics, ...
                                      results.data.simulation.sampling_rate, ...
                                      results.data.room.type, ...
                                      results.data.array.type, ...
                                      results.directories.base_dir);
                    
                    uialert(uiFig, message, 'Results Loaded', 'Icon', 'success');
                else
                    fprintf('✅ Results loaded successfully! (See details above)\n');
                end
            catch
                fprintf('✅ Results loaded successfully! (See details above)\n');
            end
            
        elseif evalin('base', 'exist(''final_signals'', ''var'')')
            % Legacy format
            final_signals = evalin('base', 'final_signals');
            [numMics, numSamples] = size(final_signals);
            
            if evalin('base', 'exist(''fs'', ''var'')')
                fs = evalin('base', 'fs');
                duration = numSamples / fs;
            else
                fs = 16000;
                duration = numSamples / fs;
            end
            
            fprintf('\n=== LEGACY SIMULATION DATA FOUND ===\n');
            fprintf('Microphones: %d\n', numMics);
            fprintf('Samples: %d\n', numSamples);
            fprintf('Duration: %.2f seconds\n', duration);
            fprintf('Est. Sampling Rate: %d Hz\n', fs);
            fprintf('=====================================\n\n');
            
        else
            fprintf('❌ No simulation results found in workspace.\n');
            fprintf('Please run a simulation first.\n');
        end
        
    catch ME
        fprintf('❌ Error loading results: %s\n', ME.message);
    end
end

function createAnalysisResultsTab(parent)
    % Create IR analysis results visualization
    grid = uigridlayout(parent, [2, 2]);
    grid.RowHeight = {'1x', 'fit'};
    grid.ColumnWidth = {'1x', 300};
    
    % Comparison plots
    plotPanel = uipanel(grid, 'Title', 'IR Comparison Plots', 'BackgroundColor', 'white');
    plotPanel.Layout.Row = 1;
    plotPanel.Layout.Column = 1;
    
    % Analysis metrics
    metricsPanel = uipanel(grid, 'Title', 'Analysis Metrics', 'BackgroundColor', [1, 0.98, 0.98]);
    metricsPanel.Layout.Row = 1;
    metricsPanel.Layout.Column = 2;
    
    metricsGrid = uigridlayout(metricsPanel, [12, 2]);
    metricsGrid.RowHeight = repmat({'fit'}, 1, 12);
    
    uilabel(metricsGrid, 'Text', 'Correlation:', 'FontWeight', 'bold');
    uilabel(metricsGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'NRMSE:');
    uilabel(metricsGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'SNR (dB):');
    uilabel(metricsGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'RT60 Error:');
    uilabel(metricsGrid, 'Text', '- %', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'C50 Diff:');
    uilabel(metricsGrid, 'Text', '- dB', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'C80 Diff:');
    uilabel(metricsGrid, 'Text', '- dB', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'Alignment:', 'FontWeight', 'bold');
    uilabel(metricsGrid, 'Text', 'Unknown', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'Time Shift:');
    uilabel(metricsGrid, 'Text', '- ms', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'Microphones:');
    uilabel(metricsGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'Quality:', 'FontWeight', 'bold');
    uilabel(metricsGrid, 'Text', 'Not analyzed', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'Analysis Date:');
    uilabel(metricsGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(metricsGrid, 'Text', 'Duration:');
    uilabel(metricsGrid, 'Text', '- sec', 'FontColor', [0.6, 0.6, 0.6]);
    
    % Control buttons
    controlPanel = uipanel(grid, 'BackgroundColor', [0.95, 0.95, 0.95]);
    controlPanel.Layout.Row = 2;
    controlPanel.Layout.Column = [1, 2];
    
    controlGrid = uigridlayout(controlPanel, [1, 4]);
    controlGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
    
    uilabel(controlGrid, 'Text', ''); % Spacer
    
    uibutton(controlGrid, 'Text', '📊 Load Analysis', 'BackgroundColor', [0.3, 0.6, 0.9], 'FontColor', 'white');
    uibutton(controlGrid, 'Text', '📈 View Plots', 'BackgroundColor', [0.2, 0.7, 0.3], 'FontColor', 'white');
    uibutton(controlGrid, 'Text', '📄 Generate Report', 'BackgroundColor', [0.6, 0.3, 0.8], 'FontColor', 'white');
end

function createBeamformingResultsTab(parent)
    % Create beamforming results visualization
    grid = uigridlayout(parent, [2, 2]);
    grid.RowHeight = {'1x', 'fit'};
    grid.ColumnWidth = {'1x', 300};
    
    % Beamforming maps
    plotPanel = uipanel(grid, 'Title', 'Acoustic Source Maps', 'BackgroundColor', 'white');
    plotPanel.Layout.Row = 1;
    plotPanel.Layout.Column = 1;
    
    % Localization results
    locPanel = uipanel(grid, 'Title', 'Localization Results', 'BackgroundColor', [0.98, 0.98, 1]);
    locPanel.Layout.Row = 1;
    locPanel.Layout.Column = 2;
    
    locGrid = uigridlayout(locPanel, [14, 2]);
    locGrid.RowHeight = repmat({'fit'}, 1, 14);
    
    uilabel(locGrid, 'Text', 'True Position:', 'FontWeight', 'bold');
    uilabel(locGrid, 'Text', '[-, -, -]', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'CBF Result:', 'FontWeight', 'bold');
    uilabel(locGrid, 'Text', '[-, -]', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'CBF Error:');
    uilabel(locGrid, 'Text', '- m', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'OMP-DAMAS:', 'FontWeight', 'bold');
    uilabel(locGrid, 'Text', '[-, -]', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'OMP Error:');
    uilabel(locGrid, 'Text', '- m', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'SC-DAMAS:', 'FontWeight', 'bold');
    uilabel(locGrid, 'Text', '[-, -]', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'SC Error:');
    uilabel(locGrid, 'Text', '- m', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'Best Algorithm:', 'FontWeight', 'bold');
    uilabel(locGrid, 'Text', 'Unknown', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'Min Error:');
    uilabel(locGrid, 'Text', '- m', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'Frequency:');
    uilabel(locGrid, 'Text', '- Hz', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'Grid Size:');
    uilabel(locGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'Processing:');
    uilabel(locGrid, 'Text', '- sec', 'FontColor', [0.6, 0.6, 0.6]);
    
    uilabel(locGrid, 'Text', 'SNR (dB):');
    uilabel(locGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
    
    % Control buttons
    controlPanel = uipanel(grid, 'BackgroundColor', [0.95, 0.95, 0.95]);
    controlPanel.Layout.Row = 2;
    controlPanel.Layout.Column = [1, 2];
    
    controlGrid = uigridlayout(controlPanel, [1, 4]);
    controlGrid.ColumnWidth = {'1x', 'fit', 'fit', 'fit'};
    
    uilabel(controlGrid, 'Text', ''); % Spacer
    
    uibutton(controlGrid, 'Text', '📡 Load Results', 'BackgroundColor', [0.3, 0.6, 0.9], 'FontColor', 'white');
    uibutton(controlGrid, 'Text', '🗺️ 3D Visualization', 'BackgroundColor', [0.6, 0.2, 0.8], 'FontColor', 'white');
    uibutton(controlGrid, 'Text', '📊 Compare Algorithms', 'BackgroundColor', [0.8, 0.6, 0.2], 'FontColor', 'white');
end



% % This completes the comprehensive GUI implementation
% % The GUI integrates all your existing modules without requiring any code changes
% % It provides a professional interface for the entire acoustic analysis workflowText', 'Latest Simulation:', 'FontWeight', 'bold');
%     uilabel(summaryGrid, 'Text', 'No simulation data', 'FontColor', [0.6, 0.6, 0.6]);
%     uilabel(summaryGrid, 'Text', 'Microphones:', 'FontWeight', 'bold');
%     uilabel(summaryGrid, 'Text', 'Latest Simulation:', 'FontWeight', 'bold');
%     uilabel(summaryGrid, 'Text', 'No simulation data', 'FontColor', [0.6, 0.6, 0.6]);
%     uilabel(summaryGrid, 'Text', 'Microphones:', 'FontWeight', 'bold');
%     uilabel(summaryGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
%     uilabel(summaryGrid, 'Text', 'Duration:', 'FontWeight', 'bold');
%     uilabel(summaryGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);
%     uilabel(summaryGrid, 'Text', 'Processing Time:', 'FontWeight', 'bold');
%     uilabel(summaryGrid, 'Text', '-', 'FontColor', [0.6, 0.6, 0.6]);

%% ============= MISSING CORE FUNCTIONS =============

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

function runDNNWorkflow(params)
% RUNDNNWORKFLOW - Unified DNN workflow for dataset generation, training, and prediction
% 
% This function handles the complete DNN workflow:
% 1. Dataset generation (if needed)
% 2. Data loading and preparation  
% 3. DNN training or prediction
% 4. Results processing and display
%
% INPUT:
%   params - Structure containing DNN parameters from GUI
%
% WORKFLOW MODES:
%   - Generate Dataset -> Train/Predict: Creates new dataset then proceeds
%   - Use Current Simulation: Uses existing simulation data
%   - Load from WAV Files: Loads existing dataset from acoustic_wavs/

fprintf('=== STARTING UNIFIED DNN WORKFLOW ===\n');

try
    %% STEP 1: DETERMINE DATA SOURCE AND WORKFLOW MODE
    
    % Check if this was triggered by "Generate Dataset" button
    if isfield(params, 'generate_dataset_requested') && params.generate_dataset_requested
        effective_data_source = 'Generate New Dataset';
        if params.training_mode
            modeStr = 'Training';
        else
            modeStr = 'Prediction';
        end
        fprintf('Mode: Generate Dataset -> %s\n', modeStr);

    else
        if isfield(params, 'use_wav_files') && params.use_wav_files
            effective_data_source = 'Load from WAV Files';
        else
            effective_data_source = 'Use Current Simulation';
        end
        
        if params.training_mode
            modeStr = 'Training';
        else
            modeStr = 'Prediction';
        end
        fprintf('Mode: %s -> %s\n', effective_data_source, modeStr);

    end
    
    %% STEP 2: DATA ACQUISITION PHASE
    
    if strcmp(effective_data_source, 'Generate New Dataset')
        %% DATASET GENERATION PATH
        fprintf('\n--- DATASET GENERATION PHASE ---\n');
        
        % INHERIT parameters from current GUI simulation setup
        if evalin('base', 'exist(''simulation_results'', ''var'')')
            sim_results = evalin('base', 'simulation_results');
            fprintf('✓ Inheriting parameters from current simulation\n');
            
            % Inherit room configuration
            room_params = sim_results.data.room;
            array_params = sim_results.data.array;
            source_params = sim_results.data.source;
            
        else
            fprintf('⚠ No current simulation found, using default parameters\n');
            % Use defaults if no simulation exists
            room_params = struct('type', 'rectangular', 'dimensions', [4, 4, 4], ...
                               'temperature', 20, 'humidity', 50);
            array_params = struct('type', 'linear', 'orientation', 'y-axis', ...
                                'num_mics', 4, 'center', [0.5, 0.5, 0.5]);
            source_params = struct('signal_type', 'measurement_sweep', 'power', 0.5, 'duration', 1);
        end
        
        % TODO: ADD GUI - Multiple source positions for dataset generation
        % Currently using default 3x2x2 grid - should be user configurable
        fprintf('Using default position grid (TODO: Add GUI control for custom positions)\n');
        [x_grid, y_grid, z_grid] = meshgrid([1.5, 2.5, 3.5], [1.5, 2.5], [1.5, 2.5]);
        positions = [x_grid(:), y_grid(:), z_grid(:)];
        
        % TODO: ADD GUI - Microphone selection mode configuration  
        % Currently using range mode with first 4 mics - should be user configurable
        mic_config = struct();
        mic_config.selection_mode = 'range'; % TODO: ADD GUI - 'all', 'custom', 'range', 'random'
        mic_config.range = [1, min(4, array_params.num_mics)]; % Use available mics
        mic_config.custom_indices = [1, 2, 3, 4]; % TODO: ADD GUI - Custom mic selection
        mic_config.num_random = 4; % TODO: ADD GUI - Number of random mics
        fprintf('Using mic selection mode: %s (TODO: Add GUI control)\n', mic_config.selection_mode);
        
        % INHERIT noise configuration from GUI (basic on/off)
        noise_config = struct();
        if isfield(params, 'enable_noise')
            noise_config.enable_noise = params.enable_noise; % From GUI  
        else
            noise_config.enable_noise = true; % Default
        end

        if isfield(params, 'noise_type')
            noise_config.noise_type = params.noise_type; % From GUI  
        else
            noise_config.noise_type = 'white'; % Default
        end
        fprintf('Noise config: %s (enabled: %s)\n', noise_config.noise_type, string(noise_config.enable_noise));
        
        % TODO: ADD GUI - Advanced noise configuration
        % - SNR levels, noise profiles, multiple noise types, etc.
        
        % Call dataset generation
        % Prepare dataset generation parameters
        % Call dataset generation
        fprintf('Calling generateAcousticDataset...\n');
        %% Debug: Print positions
        fprintf('\n=== DEBUG: POSITIONS ===\n');
        fprintf('Positions size: [%d x %d]\n', size(positions, 1), size(positions, 2));
        fprintf('Number of source positions: %d\n', size(positions, 1));
        fprintf('Position coordinates:\n');
        for i = 1:size(positions, 1)
            fprintf('  Position %d: [%.2f, %.2f, %.2f]\n', i, positions(i,1), positions(i,2), positions(i,3));
        end
        fprintf('=== END POSITIONS DEBUG ===\n\n');
        dataset_results = generateAcousticDataset(params,positions,mic_config);
        
        if ~dataset_results.success
            error('Dataset generation failed');
        end
        
        fprintf('✓ Dataset generation completed: %d WAV files\n', dataset_results.total_wav_files);
        
        % Set data source for next phase
        data_source_for_dnn = 'Load from WAV Files';
        dataset_directory = dataset_results.output_directory;
        
    elseif strcmp(effective_data_source, 'Use Current Simulation')
        %% CURRENT SIMULATION PATH  
        fprintf('\n--- USING CURRENT SIMULATION DATA ---\n');
        
        if ~evalin('base', 'exist(''final_signals'', ''var'')')
            error('No simulation data found. Run simulation first.');
        end
        
        data_source_for_dnn = 'Use Current Simulation';
        dataset_directory = '';
        fprintf('✓ Using current simulation data from workspace\n');
        
    elseif strcmp(effective_data_source, 'Load from WAV Files')
        %% LOAD EXISTING DATASET PATH
        fprintf('\n--- LOADING EXISTING DATASET ---\n');
        
        dataset_directory = 'acoustic_wavs';
        if ~exist(dataset_directory, 'dir')
            error('Dataset directory not found: %s. Generate a dataset first.', dataset_directory);
        end
        
        data_source_for_dnn = 'Load from WAV Files';
        fprintf('✓ Loading dataset from: %s\n', dataset_directory);
        
    else
        error('Unknown data source: %s', effective_data_source);
    end
    
    %% STEP 3: DATA PREPARATION PHASE
    fprintf('\n--- DATA PREPARATION PHASE ---\n');
    
    if strcmp(data_source_for_dnn, 'Use Current Simulation')
        %% PREPARE SIMULATION DATA
        fprintf('Preparing simulation data for DNN...\n');
        
        % Get simulation data from workspace
        final_signals = evalin('base', 'final_signals');
        micPositions = evalin('base', 'micPositions');
        
        % Get source position - try multiple variable names
        if evalin('base', 'exist(''external_start_pos'', ''var'')')
            sourcePos = evalin('base', 'external_start_pos');
        elseif evalin('base', 'exist(''start_pos'', ''var'')')
            sourcePos = evalin('base', 'start_pos');
        elseif evalin('base', 'exist(''sourcePos'', ''var'')')
            sourcePos = evalin('base', 'sourcePos');
        else
            sourcePos = [2, 2, 2]; % Default fallback
            fprintf('Warning: Using default source position [2, 2, 2]\n');
        end
        
        % Single position case
        source_positions = sourcePos;
        fprintf('Single source position: [%.2f, %.2f, %.2f]\n', sourcePos(1), sourcePos(2), sourcePos(3));
        
    else
        %% LOAD DATASET DATA
        fprintf('Loading dataset from WAV files...\n');
        
        [final_signals, micPositions, source_positions] = loadDatasetFromWAV(dataset_directory);
        
        if isempty(final_signals)
            error('No WAV files found in dataset directory: %s', dataset_directory);
        end
        
        fprintf('Loaded dataset: %d signals from %d positions\n', ...
            size(final_signals, 1), size(source_positions, 1));
    end
    
    %% STEP 4: DNN DATA STRUCTURE PREPARATION
    fprintf('\n--- DNN DATA STRUCTURE PREPARATION ---\n');
    
    [numMics, numSamples] = size(final_signals);
    fprintf('Processing DNN data: %d microphones, %d samples each\n', numMics, numSamples);
    
    % Ensure signal length matches expected DNN input (8001 samples)
    wav_file_length = 8001;
    if numSamples > wav_file_length
        final_signals = final_signals(:, 1:wav_file_length);
        numSamples = wav_file_length;
    elseif numSamples < wav_file_length
        % Pad with zeros
        final_signals(:, end+1:wav_file_length) = 0;
        numSamples = wav_file_length;
    end
    
    % Calculate ground truth labels for each microphone
    source_to_mic_distances = [];
    source_to_mic_azimuths = [];
    source_to_mic_elevations = [];
    position_mic_signals = [];
    
    % Handle single position vs multiple positions
    if size(source_positions, 1) == 1
        %% SINGLE POSITION CASE (Current Simulation)
        fprintf('Processing single source position: [%.2f, %.2f, %.2f]\n', ...
            source_positions(1), source_positions(2), source_positions(3));
        
        for mic = 1:numMics
            mic_pos = micPositions(mic, :);
            
            % Calculate ground truth labels
            distance = norm(source_positions - mic_pos);
            dx = source_positions(1) - mic_pos(1);
            dy = source_positions(2) - mic_pos(2);
            azimuth = atan2(dy, dx) * 180 / pi;
            dz = source_positions(3) - mic_pos(3);
            dxy = sqrt(dx^2 + dy^2);
            elevation = atan2(dz, dxy) * 180 / pi;
            
            % Append to accumulated arrays
            source_to_mic_distances = [source_to_mic_distances, distance];
            source_to_mic_azimuths = [source_to_mic_azimuths, azimuth];
            source_to_mic_elevations = [source_to_mic_elevations, elevation];
            
            % Get microphone signal and normalize
            mic_signal = final_signals(mic, :);
            if max(abs(mic_signal)) > 0
                mic_signal = 0.9 * mic_signal / max(abs(mic_signal));
            end
            
            % Collect signal data
            position_mic_signals = [position_mic_signals; mic_signal];
        end
        
    else
        %% MULTIPLE POSITIONS CASE (Generated/Loaded Dataset)
        fprintf('Processing multiple source positions: %d positions\n', size(source_positions, 1));
        
        num_positions = size(source_positions, 1);
        mics_per_position = numMics / num_positions;
        
        if mod(numMics, num_positions) ~= 0
            fprintf('Warning: Number of mics (%d) not evenly divisible by positions (%d)\n', ...
                numMics, num_positions);
            mics_per_position = floor(numMics / num_positions);
        end
        
        signal_idx = 1;
        for pos_idx = 1:num_positions
            current_source_pos = source_positions(pos_idx, :);
            
            for mic_in_pos = 1:mics_per_position
                if signal_idx > numMics
                    break;
                end
                
                % Get microphone position (assuming same array for all positions)
                mic_pos = micPositions(mic_in_pos, :);
                
                % Calculate ground truth labels
                distance = norm(current_source_pos - mic_pos);
                dx = current_source_pos(1) - mic_pos(1);
                dy = current_source_pos(2) - mic_pos(2);
                azimuth = atan2(dy, dx) * 180 / pi;
                dz = current_source_pos(3) - mic_pos(3);
                dxy = sqrt(dx^2 + dy^2);
                elevation = atan2(dz, dxy) * 180 / pi;
                
                % Append to accumulated arrays
                source_to_mic_distances = [source_to_mic_distances, distance];
                source_to_mic_azimuths = [source_to_mic_azimuths, azimuth];
                source_to_mic_elevations = [source_to_mic_elevations, elevation];
                
                % Get microphone signal and normalize
                mic_signal = final_signals(signal_idx, :);
                if max(abs(mic_signal)) > 0
                    mic_signal = 0.9 * mic_signal / max(abs(mic_signal));
                end
                
                % Collect signal data
                position_mic_signals = [position_mic_signals; mic_signal];
                
                signal_idx = signal_idx + 1;
            end
        end
    end
    
    % Create data structures that DNN scripts expect
    localization_labels = struct('azimuth', source_to_mic_azimuths, ...
                               'elevation', source_to_mic_elevations, ...
                               'distance', source_to_mic_distances);
    
    % Create localization frame [num_samples x wav_file_length x 1]
    % Create localization frame [1 x wav_file_length x num_samples] - correct format for DNN
    localization_data = zeros(1, size(position_mic_signals, 2), 1, size(position_mic_signals, 1));
    for i = 1:size(position_mic_signals, 1)
        localization_data(1, :, 1, i) = position_mic_signals(i, :);
    end
    
    fprintf('Created localization frame of size [%d x %d x 1]\n', ...
        size(localization_data, 1), size(localization_data, 2));
    
    %% STEP 5: SET WORKSPACE VARIABLES FOR DNN SCRIPTS
    fprintf('\n--- SETTING WORKSPACE VARIABLES ---\n');
    
    assignin('base', 'localization_data', localization_data);
    assignin('base', 'localization_labels', localization_labels);
    assignin('base', 'position_mic_signals', position_mic_signals);
    assignin('base', 'source_to_mic_distances', source_to_mic_distances);
    assignin('base', 'source_to_mic_azimuths', source_to_mic_azimuths);
    assignin('base', 'source_to_mic_elevations', source_to_mic_elevations);
    
    %% ADD THIS TO runDNNWorkflow() function in AcousticAnalysisGUI.m
    
    
    % Set operation mode for DNN scripts
    if params.training_mode
        assignin('base', 'operation', '1');
    else
        assignin('base', 'operation', '2');
    end
    
    % Pass training parameters to DNN scripts - COMPLETE IMPLEMENTATION
    assignin('base', 'gui_epochs', params.epochs);
    assignin('base', 'gui_batch_size', params.batch_size);
    assignin('base', 'gui_learning_rate', params.learning_rate);
    assignin('base', 'gui_validation_split', params.validation_split);
    assignin('base', 'gui_data_augmentation', params.data_augmentation);
    assignin('base', 'gui_gpu_training', params.gpu_training);
    
    fprintf('GUI parameters passed to DNN scripts:\n');
    fprintf('  Epochs: %d\n', params.epochs);
    fprintf('  Batch Size: %d\n', params.batch_size);
    fprintf('  Learning Rate: %.4f\n', params.learning_rate);
    fprintf('  Validation Split: %.2f\n', params.validation_split);
    fprintf('  Data Augmentation: %s\n', string(params.data_augmentation));
    fprintf('  GPU Training: %s\n', string(params.gpu_training));
    
    %% STEP 6: EXECUTE DNN TRAINING OR PREDICTION
    fprintf('\n--- DNN EXECUTION PHASE ---\n');
    
    % Initialize prediction variables
    predicted_azimuth = [];
    predicted_distance = [];
    predicted_elevation = [];
    
    if params.training_mode
        %% TRAINING MODE
        if params.train_combined
            %% Combined Training Mode
            assignin('base', 'operation', '3');
            fprintf('Starting combined 3D DNN training...\n');
            fprintf('Training combined azimuth+distance+elevation model...\n');
            evalin('base', 'run(''DNN_combined_3D_localization.m'')');
            fprintf('✓ Combined 3D model trained successfully!\n');
            
        else
            %% Individual Training Mode
            assignin('base', 'operation', '1');
            fprintf('Starting individual DNN training...\n');
            
            if params.train_azimuth
                fprintf('Training azimuth estimation model...\n');
                evalin('base', 'run(''DNN_azimuth_est_clean_noisy1.m'')');
                fprintf('✓ Azimuth model trained\n');
            end
            
            if params.train_distance
                fprintf('Training distance estimation model...\n');
                evalin('base', 'run(''DNN_distance_est_clean_noisy1.m'')');
                fprintf('✓ Distance model trained\n');
            end
            
            if params.train_elevation
                fprintf('Training elevation estimation model...\n');
                evalin('base', 'run(''DNN_elevation_est_clean_noisy1.m'')');
                fprintf('✓ Elevation model trained\n');
            end
            
            fprintf('✓ All individual DNN models trained successfully!\n');
        end
        
    else
        %% PREDICTION MODE
        if params.train_combined
            %% Combined Prediction Mode
            assignin('base', 'operation', '4');
            fprintf('Starting combined 3D DNN prediction...\n');
            
            % Check if combined model exists
            if ~exist('combined_3D_model.mat', 'file')
                error('Combined 3D model not found: combined_3D_model.mat. Please train the combined model first.');
            end
            
            % Load trained combined model
            load('combined_3D_model.mat', 'trainedNetwork_combined');
            assignin('base', 'trainedNetwork_combined', trainedNetwork_combined);
            
            % Make combined prediction
            % Predict each sample individually for combined model
            predicted_combined = zeros(size(localization_data, 4), 3);
            for i = 1:size(localization_data, 4)
                single_sample = localization_data(:, :, :, i);
                predicted_combined(i, :) = predict(trainedNetwork_combined, single_sample);
            end
            
            % Extract individual components from combined output
            predicted_azimuth = predicted_combined(:, 1);
            predicted_distance = predicted_combined(:, 2);
            predicted_elevation = predicted_combined(:, 3);
            
            fprintf('✓ Combined 3D predictions completed successfully!\n');
            
        else
            %% Individual Prediction Mode
            assignin('base', 'operation', '2');
            fprintf('Starting individual DNN prediction...\n');
            
            % Check if trained models exist
            model_files = {'azimuth_model.mat', 'distance_model.mat', 'elevation_model.mat'};
            models_exist = all(cellfun(@(x) exist(x, 'file'), model_files));
            
            if ~models_exist
                missing_models = model_files(~cellfun(@(x) exist(x, 'file'), model_files));
                error('Missing trained models: %s. Please train models first.', strjoin(missing_models, ', '));
            end
            
            % Load trained individual models
            load('azimuth_model.mat', 'trainedNetwork_azimuth');
            load('distance_model.mat', 'trainedNetwork_distance');
            load('elevation_model.mat', 'trainedNetwork_elevation');
            
            % Make individual predictions
            % Predict each sample individually
            predicted_azimuth = zeros(size(localization_data, 4), 1);
            for i = 1:size(localization_data, 4)
                single_sample = localization_data(:, :, :, i);
                predicted_azimuth(i) = predict(trainedNetwork_azimuth, single_sample);
            end
            
            % Predict each sample individually
            predicted_distance = zeros(size(localization_data, 4), 1);
            for i = 1:size(localization_data, 4)
                single_sample = localization_data(:, :, :, i);
                predicted_distance(i) = predict(trainedNetwork_distance, single_sample);
            end

            predicted_elevation = zeros(size(localization_data, 4), 1);
            for i = 1:size(localization_data, 4)
                single_sample = localization_data(:, :, :, i);
                predicted_elevation(i) = predict(trainedNetwork_elevation, single_sample);
            end
            
            fprintf('✓ Individual predictions completed successfully!\n');
        end
        
        %% STEP 7: PREDICTION RESULTS ANALYSIS (only for prediction modes)
        if ~isempty(predicted_azimuth)
            fprintf('\n--- PREDICTION RESULTS ANALYSIS ---\n');
            
            % Display prediction results
            fprintf('\nPrediction Results:\n');
            fprintf('%-8s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
                'Sample', 'True Az', 'Pred Az', 'True El', 'Pred El', 'True Dist', 'Pred Dist');
            fprintf('%-8s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
                '------', '--------', '--------', '--------', '--------', '--------', '--------');
            
            for i = 1:min(10, length(localization_labels.azimuth)) % Show first 10 results
                fprintf('%-8d %-12.1f %-12.1f %-12.1f %-12.1f %-12.2f %-12.2f\n', ...
                    i, localization_labels.azimuth(i), predicted_azimuth(i), ...
                    localization_labels.elevation(i), predicted_elevation(i), ...
                    localization_labels.distance(i), predicted_distance(i));
            end
            
            if length(localization_labels.azimuth) > 10
                fprintf('... (showing first 10 of %d results)\n', length(localization_labels.azimuth));
            end
            
            % Calculate individual parameter errors
            az_error = mean(abs(predicted_azimuth - localization_labels.azimuth'));
            el_error = mean(abs(predicted_elevation - localization_labels.elevation'));
            dist_error = mean(abs(predicted_distance - localization_labels.distance'));
            
            fprintf('\nMean Absolute Errors:\n');
            fprintf('  Azimuth: %.2f degrees\n', az_error);
            fprintf('  Elevation: %.2f degrees\n', el_error);
            fprintf('  Distance: %.3f meters\n', dist_error);
            
            % Save prediction results
            prediction_results = struct();
            prediction_results.true_azimuth = localization_labels.azimuth;
            prediction_results.predicted_azimuth = predicted_azimuth';
            prediction_results.true_elevation = localization_labels.elevation;
            prediction_results.predicted_elevation = predicted_elevation';
            prediction_results.true_distance = localization_labels.distance;
            prediction_results.predicted_distance = predicted_distance';
            prediction_results.errors = struct('azimuth', az_error, 'elevation', el_error, 'distance', dist_error);
            prediction_results.data_source = effective_data_source;
            prediction_results.num_samples = length(localization_labels.azimuth);
            
            % Save results to file
            if params.train_combined
                results_file = 'combined_dnn_prediction_results.mat';
            else
                results_file = 'individual_dnn_prediction_results.mat';
            end
            
            save(results_file, 'prediction_results');
            fprintf('✓ Prediction results saved to: %s\n', results_file);
        end
    end
    
    %% STEP 8: COMPLETION SUMMARY
    fprintf('\n=======================================================\n');
    fprintf('DNN WORKFLOW COMPLETED SUCCESSFULLY!\n');
    fprintf('=======================================================\n');
    fprintf('Data source: %s\n', effective_data_source);
    if params.training_mode
        operationStr = 'Training';
    else
        operationStr = 'Prediction';
    end
    fprintf('Operation: %s\n', operationStr);

    fprintf('Training samples: %d\n', length(source_to_mic_azimuths));
    
    if strcmp(effective_data_source, 'Generate New Dataset')
        fprintf('Dataset files: %d WAV files created\n', dataset_results.total_wav_files);
        fprintf('Dataset directory: %s\n', dataset_results.output_directory);
    end
    
    if params.train_combined
        fprintf('Model type: Combined 3D\n');
    else
        models_trained = {};
        if params.train_azimuth, models_trained{end+1} = 'Azimuth'; end
        if params.train_distance, models_trained{end+1} = 'Distance'; end
        if params.train_elevation, models_trained{end+1} = 'Elevation'; end
        fprintf('Model type: Individual (%s)\n', strjoin(models_trained, ', '));
    end
    
    fprintf('=======================================================\n');
    
catch ME
    fprintf('\n❌ DNN WORKFLOW FAILED ❌\n');
    fprintf('Error: %s\n', ME.message);
    fprintf('Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    rethrow(ME);
end

end

%% ============= HELPER FUNCTION =============

function [final_signals, micPositions, source_positions] = loadDatasetFromWAV(dataset_dir)
% Load dataset from WAV files and metadata
% (This function already exists in the original code, included here for completeness)

final_signals = [];
micPositions = [];
source_positions = [];

% Check if metadata file exists
metadata_file = fullfile(dataset_dir, 'sound_source_metadata.csv');
if ~exist(metadata_file, 'file')
    fprintf('Warning: Metadata file not found: %s\n', metadata_file);
    return;
end

% Read metadata
try
    metadata = readtable(metadata_file);
    fprintf('Loaded metadata with %d entries\n', height(metadata));
catch ME
    fprintf('Error reading metadata: %s\n', ME.message);
    return;
end

% Extract unique source positions
unique_sources = unique(metadata(:, {'source_x', 'source_y', 'source_z'}), 'rows');
source_positions = table2array(unique_sources);

% Extract unique microphone positions
unique_mics = unique(metadata(:, {'mic_x', 'mic_y', 'mic_z'}), 'rows');
micPositions = table2array(unique_mics);

fprintf('Found %d unique source positions and %d microphone positions\n', ...
    size(source_positions, 1), size(micPositions, 1));

% Load WAV files
wav_files = metadata.file_name;
num_files = length(wav_files);

if num_files == 0
    fprintf('No WAV files found in metadata\n');
    return;
end

% Read first file to get dimensions
first_file = fullfile(dataset_dir, wav_files{1});
if exist(first_file, 'file')
    [sample_signal, fs] = audioread(first_file);
    signal_length = length(sample_signal);
    fprintf('WAV file format: %d samples at %d Hz\n', signal_length, fs);
else
    fprintf('Error: First WAV file not found: %s\n', first_file);
    return;
end

% Initialize signal matrix
final_signals = zeros(num_files, signal_length);

% Load all WAV files
loaded_count = 0;
for i = 1:num_files
    wav_file = fullfile(dataset_dir, wav_files{i});
    
    if exist(wav_file, 'file')
        try
            [signal, ~] = audioread(wav_file);
            
            % Ensure consistent length
            if length(signal) == signal_length
                final_signals(i, :) = signal';
                loaded_count = loaded_count + 1;
            else
                fprintf('Warning: Signal length mismatch in %s\n', wav_files{i});
            end
        catch ME
            fprintf('Warning: Could not load %s: %s\n', wav_files{i}, ME.message);
        end
    else
        fprintf('Warning: WAV file not found: %s\n', wav_file);
    end
end

fprintf('Successfully loaded %d out of %d WAV files\n', loaded_count, num_files);

% Remove empty rows
if loaded_count < num_files
    non_zero_rows = any(final_signals, 2);
    final_signals = final_signals(non_zero_rows, :);
    fprintf('Removed %d empty signal rows\n', num_files - loaded_count);
end

end




function runExternalIRAnalysisComplete()
    % Complete implementation of IR analysis
    try
        fprintf('Running impulse response comparison analysis...\n');
        
        % Get the selected files from base workspace
        measured_file = evalin('base', 'selected_measured_file');
        simulated_file = evalin('base', 'selected_simulated_file');

        % Get the output folder name
        output_folder = IRs_comparison_analysis_tool3(measured_file, simulated_file);
        
        % Store it for later use
        assignin('base', 'latest_ir_analysis_folder', output_folder);
        
        fprintf('✅ IR analysis completed! Results in: %s\n', output_folder);
    catch ME
        error('IR analysis failed: %s', ME.message);
    end
end

function runExternalBeamformingComplete(params, panel)
    try
        fprintf('🔍 Starting beamforming...\n');
        
        % Get all data from base workspace
        final_signals = evalin('base', 'final_signals');
        micPositions = evalin('base', 'micPositions');
        fs = evalin('base', 'fs');
        sourcePos = evalin('base', 'sourcePos');
        
        % Set ALL base workspace variables
        assignin('base', 'noisy_signals', final_signals);
        assignin('base', 'numMics', size(final_signals, 1));
        assignin('base', 't', (0:size(final_signals,2)-1)/fs);
        assignin('base', 'freq_vector', linspace(0, fs/2, floor(size(final_signals,2)/2)+1));
        assignin('base', 'H', ones(size(final_signals, 1), floor(size(final_signals,2)/2)+1));
        assignin('base', 'roomDim', [4, 4, 4]);
        assignin('base', 'SNR_dB', 10);
        assignin('base', 'freq_to_analyze', params.frequency);
        assignin('base', 'grid_size', params.grid_size);
        assignin('base', 'scan_area_size', params.scan_area_size);
        assignin('base', 'max_iterations', params.max_iterations);
        assignin('base', 'delta_criteria', params.delta_criteria);
        assignin('base', 'sparsity_step', params.sparsity_step);
        assignin('base', 'plot_results', true);
        assignin('base', 'save_results', true);
        
        % External call handshake
        assignin('base', 'EXTERNAL_CALL_FLAG', true);
        try
            evalin('base', 'run(''room_beamforming_and_comparison.m'')');
        catch ME
            if ~strcmp(ME.message, 'Beamforming_check_completed')
                error('Beamforming check failed: %s', ME.message);
            end
        end
        
        % Real beamforming run
        evalin('base', 'clear EXTERNAL_CALL_FLAG');
        evalin('base', 'run(''room_beamforming_and_comparison.m'')');
        
        % Set results folder
        assignin('base', 'latest_beamforming_folder', 'beamforming_results');
        
        fprintf('✅ Beamforming completed!\n');
        
    catch ME
        error('Beamforming analysis failed: %s', ME.message);
    end
end


function loadAnalysisResultsComplete(panel)
    try
        % Look for latest comparison results in correct location
        resultFiles = dir('IR_comparison_results_*/IR_comparison_results_*.mat');
        if ~isempty(resultFiles)
            % Sort by date and load latest
            [~, idx] = sort([resultFiles.datenum], 'descend');
            latestFile = fullfile(resultFiles(idx(1)).folder, resultFiles(idx(1)).name);
            
            data = load(latestFile);
            if isfield(data, 'comparison_results')
                results = data.comparison_results;
                
                % Your data structure has single-level statistical analysis (not per-microphone)
                % Extract the actual values from the loaded data
                correlation = results.statistical_analysis.correlation_coefficient;
                nrmse = results.statistical_analysis.nrmse;
                snr = results.statistical_analysis.snr_db;
                
                % Update GUI labels with actual data
                panel.UserData.corrLabel.Text = sprintf('%.3f', correlation);
                panel.UserData.nrmseLabel.Text = sprintf('%.4f', nrmse);
                panel.UserData.snrLabel.Text = sprintf('%.1f', snr);
                
                % Set microphones to 1 since this is single-channel analysis
                panel.UserData.micsLabel.Text = '1';
                
                % Check for alignment quality (may not exist in current structure)
                if isfield(results.time_domain, 'cross_correlation')
                    alignment_quality = results.time_domain.cross_correlation.max_correlation;
                    panel.UserData.alignLabel.Text = sprintf('%.3f', alignment_quality);
                else
                    panel.UserData.alignLabel.Text = 'N/A';
                end
                
                % Determine overall quality based on actual thresholds
                if correlation > 0.8 && nrmse < 0.2
                    quality = 'EXCELLENT';
                    color = [0, 0.6, 0];
                elseif correlation > 0.6 && nrmse < 0.4
                    quality = 'GOOD';
                    color = [0, 0.6, 0];
                elseif correlation > 0.3 && nrmse < 0.8
                    quality = 'FAIR';
                    color = [0.8, 0.6, 0];
                else
                    quality = 'POOR';
                    color = [0.8, 0, 0];
                end
                
                panel.UserData.qualityLabel.Text = quality;
                panel.UserData.qualityLabel.FontColor = color;
                
                % Create comprehensive summary text
                summaryText = {
                    sprintf('✅ Analysis completed at: %s', char(results.metadata.analysis_timestamp));
                    sprintf('📁 Results folder: %s', resultFiles(idx(1)).folder);
                    '';
                    '📊 STATISTICAL METRICS:';
                    sprintf('   Correlation: %.3f', correlation);
                    sprintf('   NRMSE: %.4f', nrmse);
                    sprintf('   SNR: %.1f dB', snr);
                    sprintf('   RMSE: %.4f', results.statistical_analysis.rmse);
                    '';
                    '⏱️ TIME DOMAIN:';
                    sprintf('   Max Cross-Correlation: %.3f', results.time_domain.cross_correlation.max_correlation);
                    sprintf('   Time Delay: %.1f ms', results.time_domain.cross_correlation.time_delay * 1000);
                    sprintf('   RMS Ratio: %.3f (%.1f dB)', results.time_domain.rms.ratio, results.time_domain.rms.difference_db);
                    '';
                    '🔊 ACOUSTIC METRICS:';
                };
                
                % Add RT60 data if available
                if isfield(results.acoustic_metrics, 'rt60')
                    freq_bands = results.metadata.frequency_bands;
                    summaryText{end+1} = '   RT60 Analysis:';
                    for i = 1:length(freq_bands)
                        band_name = sprintf('f_%dHz', freq_bands(i));
                        if isfield(results.acoustic_metrics.rt60, band_name)
                            rt60_data = results.acoustic_metrics.rt60.(band_name);
                            summaryText{end+1} = sprintf('     %d Hz: %.2fs → %.2fs (%.1f%% error)', ...
                                freq_bands(i), rt60_data.measured, rt60_data.simulated, rt60_data.relative_error);
                        end
                    end
                    summaryText{end+1} = '';
                end
                
                % Add clarity metrics if available
                if isfield(results.acoustic_metrics, 'clarity')
                    summaryText{end+1} = '   Clarity Metrics:';
                    summaryText{end+1} = sprintf('     C50: %.1f dB → %.1f dB (%.1f dB diff)', ...
                        results.acoustic_metrics.clarity.c50.measured, ...
                        results.acoustic_metrics.clarity.c50.simulated, ...
                        results.acoustic_metrics.clarity.c50.difference);
                    summaryText{end+1} = sprintf('     C80: %.1f dB → %.1f dB (%.1f dB diff)', ...
                        results.acoustic_metrics.clarity.c80.measured, ...
                        results.acoustic_metrics.clarity.c80.simulated, ...
                        results.acoustic_metrics.clarity.c80.difference);
                    summaryText{end+1} = '';
                end
                
                % Add final instructions
                summaryText{end+1} = '📈 Ready for visualization!';
                summaryText{end+1} = 'Use "View Plots" to see detailed analysis graphs.';
                
                panel.UserData.resultsText.Value = summaryText;
                
                fprintf('✅ Analysis results loaded successfully!\n');
                fprintf('   Correlation: %.3f, NRMSE: %.4f, Quality: %s\n', correlation, nrmse, quality);
            else
                error('comparison_results field not found in loaded data');
            end
        else
            fprintf('⚠ No comparison results found.\n');
            panel.UserData.resultsText.Value = {
                '⚠ No analysis results found.';
                '';
                'Please run IR analysis first using the';
                '"Run Analysis" button.';
                '';
                'Make sure to select both measured and';
                'simulated WAV files before running.';
            };
        end
    catch ME
        fprintf('❌ Failed to load analysis results: %s\n', ME.message);
        fprintf('   Stack: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        
        % Update with helpful error message
        panel.UserData.resultsText.Value = {
            '❌ Failed to load analysis results';
            '';
            sprintf('Error: %s', ME.message);
            '';
            'Possible causes:';
            '• Analysis not yet completed';
            '• Results file corrupted';
            '• Path permissions issue';
            '';
            'Try running the analysis again.';
        };
        
        % Reset all labels to default
        panel.UserData.corrLabel.Text = '-';
        panel.UserData.nrmseLabel.Text = '-';
        panel.UserData.snrLabel.Text = '-';
        panel.UserData.alignLabel.Text = '-';
        panel.UserData.qualityLabel.Text = 'Error';
        panel.UserData.qualityLabel.FontColor = [0.8, 0, 0];
        panel.UserData.micsLabel.Text = '-';
    end
end

%% AFTER - Enhanced results loading with folder support
function loadBeamformingResultsComplete(panel)
% BULLETPROOF VERSION - NO MORE INDEX ERRORS!
fprintf('🔧 BULLETPROOF: Starting result loading...\n');

try
    % Check if results folder exists
    if ~evalin('base', 'exist(''latest_beamforming_folder'', ''var'')')
        fprintf('❌ No beamforming folder found\n');
        panel.UserData.progressText.Text = '❌ No results folder found';
        panel.UserData.progressText.FontColor = [0.8, 0, 0];
        return;
    end
    
    results_folder = evalin('base', 'latest_beamforming_folder');
    fprintf('📁 Checking folder: %s\n', results_folder);
    
    % Find .mat files
    results_files = dir(fullfile(results_folder, 'sc_damas_results_*.mat'));
    if isempty(results_files)
        fprintf('❌ No .mat files found\n');
        panel.UserData.progressText.Text = '❌ No result files found';
        panel.UserData.progressText.FontColor = [0.8, 0, 0];
        return;
    end
    
    % Load the file
    latest_file = fullfile(results_folder, results_files(1).name);
    fprintf('📂 Loading: %s\n', latest_file);
    data = load(latest_file);
    
    % BULLETPROOF table building - NO ASSUMPTIONS!
    tableData = {};
    fprintf('🔍 Building results table...\n');
    
    % CBF - SUPER SAFE
    try
        if isfield(data, 'cbf_peak_pos') && isfield(data, 'cbf_errors')
            pos = data.cbf_peak_pos;
            err = data.cbf_errors;
            
            % Check if we have data
            if ~isempty(pos) && ~isempty(err) && numel(pos) >= 2 && numel(err) >= 1
                tableData{end+1, 1} = 'CBF';
                tableData{end, 2} = sprintf('%.2f', pos(1));
                tableData{end, 3} = sprintf('%.2f', pos(2));
                tableData{end, 4} = sprintf('%.3f', err(1));
                fprintf('✅ CBF data added successfully\n');
            else
                fprintf('⚠️ CBF data exists but wrong size - skipping\n');
            end
        else
            fprintf('⚠️ CBF fields missing - skipping\n');
        end
    catch ME
        fprintf('❌ CBF processing failed: %s\n', ME.message);
    end
    
    % OMP-DAMAS - SUPER SAFE
    try
        if isfield(data, 'omp_peak_pos') && isfield(data, 'omp_errors')
            pos = data.omp_peak_pos;
            err = data.omp_errors;
            
            if ~isempty(pos) && ~isempty(err) && numel(pos) >= 2 && numel(err) >= 1
                tableData{end+1, 1} = 'OMP-DAMAS';
                tableData{end, 2} = sprintf('%.2f', pos(1));
                tableData{end, 3} = sprintf('%.2f', pos(2));
                tableData{end, 4} = sprintf('%.3f', err(1));
                fprintf('✅ OMP-DAMAS data added successfully\n');
            else
                fprintf('⚠️ OMP-DAMAS data exists but wrong size - skipping\n');
            end
        else
            fprintf('⚠️ OMP-DAMAS fields missing - skipping\n');
        end
    catch ME
        fprintf('❌ OMP-DAMAS processing failed: %s\n', ME.message);
    end
    
    % SC-DAMAS - SUPER SAFE
    try
        if isfield(data, 'sc_peak_pos') && isfield(data, 'sc_errors')
            pos = data.sc_peak_pos;
            err = data.sc_errors;
            
            if ~isempty(pos) && ~isempty(err) && numel(pos) >= 2 && numel(err) >= 1
                tableData{end+1, 1} = 'SC-DAMAS';
                tableData{end, 2} = sprintf('%.2f', pos(1));
                tableData{end, 3} = sprintf('%.2f', pos(2));
                tableData{end, 4} = sprintf('%.3f', err(1));
                fprintf('✅ SC-DAMAS data added successfully\n');
            else
                fprintf('⚠️ SC-DAMAS data exists but wrong size - skipping\n');
            end
        else
            fprintf('⚠️ SC-DAMAS fields missing - skipping\n');
        end
    catch ME
        fprintf('❌ SC-DAMAS processing failed: %s\n', ME.message);
    end
    
    % Update GUI table SAFELY
    try
        panel.UserData.resultsTable.Data = tableData;
        fprintf('📊 Table updated with %d rows\n', size(tableData, 1));
    catch ME
        fprintf('❌ Table update failed: %s\n', ME.message);
    end
    
    % Find best algorithm SAFELY
    try
        if ~isempty(tableData) && size(tableData, 2) >= 4
            errorValues = [];
            for i = 1:size(tableData, 1)
                try
                    errorValues(i) = str2double(tableData{i, 4});
                catch
                    errorValues(i) = inf; % Bad data = worst score
                end
            end
            
            if ~isempty(errorValues)
                [minError, minIdx] = min(errorValues);
                panel.UserData.bestAlgoLabel.Text = tableData{minIdx, 1};
                panel.UserData.minErrorLabel.Text = sprintf('%.3f m', minError);
                fprintf('🏆 Best algorithm: %s (%.3f m)\n', tableData{minIdx, 1}, minError);
            end
        end
    catch ME
        fprintf('❌ Best algorithm calculation failed: %s\n', ME.message);
    end
    
    % Update source position SAFELY
    try
        if isfield(data, 'beamforming_metadata') && isfield(data.beamforming_metadata, 'sourcePos')
            sourcePos = data.beamforming_metadata.sourcePos;
            if numel(sourcePos) >= 3
                panel.UserData.sourcePositionLabel.Text = sprintf('[%.2f, %.2f, %.2f]', sourcePos(1), sourcePos(2), sourcePos(3));
                fprintf('📍 Source position updated\n');
            end
        end
    catch ME
        fprintf('❌ Source position update failed: %s\n', ME.message);
    end
    
    % Final status
    panel.UserData.progressText.Text = '✅ Results loaded successfully';
    panel.UserData.progressText.FontColor = [0, 0.6, 0];
    fprintf('🎉 BULLETPROOF LOADING COMPLETE!\n');
    
catch ME
    fprintf('💥 TOTAL FAILURE: %s\n', ME.message);
    panel.UserData.progressText.Text = '❌ Failed to load results';
    panel.UserData.progressText.FontColor = [0.8, 0, 0];
end
end



function safeExecuteModule(moduleFunction, errorTitle)
    % Safe execution wrapper for all module calls
    try
        moduleFunction();
    catch ME
        fprintf('❌ %s failed: %s\n', errorTitle, ME.message);
        
        % Show user-friendly error dialog
        errorMsg = sprintf('%s failed:\n\n%s\n\nPlease check:\n• Required files are present\n• Previous steps completed\n• MATLAB path is correct', ...
                          errorTitle, ME.message);
        
        uialert(gcf, errorMsg, [errorTitle ' Error'], 'Icon', 'error');
        rethrow(ME);
    end
end

function config = showPredictionConfigDialog(selected_models, missing_models, wav_count)
   fig = uifigure('Name', 'Prediction Configuration', ...
                  'Position', [300, 300, 450, 350], ...
                  'WindowStyle', 'modal');
   
   grid = uigridlayout(fig, [4, 1]);
   grid.RowHeight = {100, 120, 80, 50};
   
   % Status Panel
   statusPanel = uipanel(grid, 'Title', 'Model & Data Status');
   statusGrid = uigridlayout(statusPanel, [3, 2]);
   
   uilabel(statusGrid, 'Text', 'Available Models:', 'FontWeight', 'bold');
   if isempty(selected_models)
       uilabel(statusGrid, 'Text', 'None selected/found', 'FontColor', [0.8, 0, 0]);
   else
       uilabel(statusGrid, 'Text', strjoin(selected_models, ', '), 'FontColor', [0, 0.6, 0]);
   end
   
   uilabel(statusGrid, 'Text', 'Missing Models:', 'FontWeight', 'bold');
   if isempty(missing_models)
       uilabel(statusGrid, 'Text', 'All required found', 'FontColor', [0, 0.6, 0]);
   else
       uilabel(statusGrid, 'Text', strjoin(missing_models, ', '), 'FontColor', [0.8, 0, 0]);
   end
   
   uilabel(statusGrid, 'Text', 'WAV Files Available:', 'FontWeight', 'bold');
   if wav_count > 0
       wav_color = [0, 0.6, 0];
   else
       wav_color = [0.8, 0, 0];
   end
   uilabel(statusGrid, 'Text', sprintf('%d files', wav_count), 'FontColor', wav_color);
   

   % Data Source Panel
    dataPanel = uipanel(grid, 'Title', 'Data Source Selection');
    dataGrid = uigridlayout(dataPanel, [1, 1]);
    
    dataBG = uibuttongroup(dataGrid, 'Title', 'Select Data Source');
    currentRB = uiradiobutton(dataBG, 'Text', 'Use Current Simulation', ...
        'Position', [10, 60, 200, 22], 'Value', true);
    
    if wav_count > 0
        wav_enable = 'on';
    else
        wav_enable = 'off';
    end
    wavRB = uiradiobutton(dataBG, 'Text', 'Load from WAV Files', ...
        'Position', [10, 30, 200, 22], 'Enable', wav_enable);
   
   % Output Options Panel  
   outputPanel = uipanel(grid, 'Title', 'Output Options');
   outputGrid = uigridlayout(outputPanel, [2, 2]);
   
   saveResultsCB = uicheckbox(outputGrid, 'Text', 'Save Results to File', 'Value', true);
   showPlotsCB = uicheckbox(outputGrid, 'Text', 'Show Prediction Plots', 'Value', true);
   
   % Buttons
   buttonPanel = uipanel(grid);
   buttonGrid = uigridlayout(buttonPanel, [1, 3]);
   buttonGrid.ColumnWidth = {'1x', 'fit', 'fit'};
   
   uilabel(buttonGrid, 'Text', '');
   cancelBtn = uibutton(buttonGrid, 'Text', 'Cancel');
   predictBtn = uibutton(buttonGrid, 'Text', 'Start Prediction', 'FontWeight', 'bold');
   
   % Enable predict button only if models and data available
   can_predict = ~isempty(selected_models) && (wav_count > 0 || true); % current sim always available
   if can_predict
       predictBtn.Enable = 'on';
   else
       predictBtn.Enable = 'off';
   end
   
   config = [];
   predictBtn.ButtonPushedFcn = @(~,~) acceptPrediction();
   cancelBtn.ButtonPushedFcn = @(~,~) close(fig);
   
   uiwait(fig);
   
   function acceptPrediction()
       config = struct();
       config.selected_models = selected_models;
       config.missing_models = missing_models;
       config.wav_count = wav_count;
       config.use_current_sim = currentRB.Value;
       config.use_wav_files = wavRB.Value;
       config.save_results = saveResultsCB.Value;
       config.show_plots = showPlotsCB.Value;
       close(fig);
   end
end

function runPredictionProcess(prediction_config, params)
    % This function will handle the actual prediction workflow
    fprintf('Starting prediction with config:\n');
    fprintf('  Selected models: %s\n', strjoin(prediction_config.selected_models, ', '));
    fprintf('  Using current sim: %s\n', string(prediction_config.use_current_sim));
    fprintf('  WAV files: %d\n', prediction_config.wav_count);
    
    % Set prediction mode parameters
    params.training_mode = false; % This is prediction, not training
    params.show_plots = prediction_config.show_plots; % Pass plot flag
    
    % Call the unified DNN workflow
    try
        runDNNWorkflow(params);
        
        if prediction_config.show_plots
            fprintf('Creating prediction plots...\n');
            
            % Create plots if prediction results exist
            if exist('individual_dnn_prediction_results.mat', 'file') || exist('combined_dnn_prediction_results.mat', 'file')
                createPredictionPlots();
            else
                fprintf('No prediction results found to plot.\n');
            end
        end
        
        if prediction_config.save_results
            fprintf('Results saved to file.\n');
        end
        
    catch ME
        uialert(ancestor(gcf, 'figure'), sprintf('Prediction failed: %s', ME.message), ...
               'Prediction Error', 'Icon', 'error');
    end
end

function createPredictionPlots()
    % Load prediction results and create plots
    if exist('individual_dnn_prediction_results.mat', 'file')
        load('individual_dnn_prediction_results.mat', 'prediction_results');
        
        figure('Position', [100, 100, 1200, 400]);
        
        % Azimuth plot
        subplot(1, 3, 1);
        scatter(prediction_results.true_azimuth, prediction_results.predicted_azimuth, 50, 'b', 'filled');
        hold on;
        plot([min(prediction_results.true_azimuth) max(prediction_results.true_azimuth)], ...
             [min(prediction_results.true_azimuth) max(prediction_results.true_azimuth)], 'r--', 'LineWidth', 2);
        xlabel('True Azimuth (degrees)');
        ylabel('Predicted Azimuth (degrees)');
        title(sprintf('Azimuth (Error: %.1f°)', prediction_results.errors.azimuth));
        grid on;
        
        % Distance plot
        subplot(1, 3, 2);
        scatter(prediction_results.true_distance, prediction_results.predicted_distance, 50, 'g', 'filled');
        hold on;
        plot([min(prediction_results.true_distance) max(prediction_results.true_distance)], ...
             [min(prediction_results.true_distance) max(prediction_results.true_distance)], 'r--', 'LineWidth', 2);
        xlabel('True Distance (m)');
        ylabel('Predicted Distance (m)');
        title(sprintf('Distance (Error: %.2f m)', prediction_results.errors.distance));
        grid on;
        
        % Elevation plot
        subplot(1, 3, 3);
        scatter(prediction_results.true_elevation, prediction_results.predicted_elevation, 50, 'm', 'filled');
        hold on;
        plot([min(prediction_results.true_elevation) max(prediction_results.true_elevation)], ...
             [min(prediction_results.true_elevation) max(prediction_results.true_elevation)], 'r--', 'LineWidth', 2);
        xlabel('True Elevation (degrees)');
        ylabel('Predicted Elevation (degrees)');
        title(sprintf('Elevation (Error: %.1f°)', prediction_results.errors.elevation));
        grid on;
        
        fprintf('Prediction plots displayed successfully!\n');
    end
end

function debugBeamformingVariables()
% Debug function to check all beamforming-related variables in base workspace
fprintf('\n=== BEAMFORMING VARIABLES DEBUG ===\n');

% List of expected variables from beamforming
expected_vars = {
    'latest_beamforming_folder', 'beamforming_results', 'cbf_results', 
    'omp_damas_results', 'sc_damas_results', 'cbf_peak_pos', 'omp_peak_pos', 
    'sc_peak_pos', 'cbf_errors', 'omp_errors', 'sc_errors', 'beamforming_metadata'
};

fprintf('Checking %d expected variables...\n\n', length(expected_vars));

for i = 1:length(expected_vars)
    var_name = expected_vars{i};
    
    if evalin('base', sprintf('exist(''%s'', ''var'')', var_name))
        % Variable exists - get details
        try
            var_value = evalin('base', var_name);
            var_class = class(var_value);
            var_size = size(var_value);
            
            fprintf('✅ %s:\n', var_name);
            fprintf('   Class: %s\n', var_class);
            fprintf('   Size: [%s]\n', num2str(var_size));
            
            % Special handling for different variable types
            if isstruct(var_value)
                field_names = fieldnames(var_value);
                fprintf('   Fields: %s\n', strjoin(field_names, ', '));
                
                % Show critical field values
                if strcmp(var_name, 'beamforming_metadata') && isfield(var_value, 'sourcePos')
                    fprintf('   sourcePos: [%.2f, %.2f, %.2f]\n', var_value.sourcePos(1), var_value.sourcePos(2), var_value.sourcePos(3));
                end
                
            elseif isnumeric(var_value)
                if numel(var_value) <= 10
                    fprintf('   Value: %s\n', mat2str(var_value));
                else
                    fprintf('   Range: [%.3f to %.3f]\n', min(var_value(:)), max(var_value(:)));
                end
                
            elseif ischar(var_value) || isstring(var_value)
                fprintf('   Value: "%s"\n', char(var_value));
                
            end
            
        catch ME
            fprintf('❌ %s: ERROR reading - %s\n', var_name, ME.message);
        end
        
    else
        fprintf('❌ %s: NOT FOUND\n', var_name);
    end
    
    fprintf('\n');
end

% Check for .mat files in results folder
fprintf('=== CHECKING RESULTS FILES ===\n');
if evalin('base', 'exist(''latest_beamforming_folder'', ''var'')')
    results_folder = evalin('base', 'latest_beamforming_folder');
    fprintf('Results folder: %s\n', results_folder);
    
    if exist(results_folder, 'dir')
        % Check for .mat files
        mat_files = dir(fullfile(results_folder, '*.mat'));
        fprintf('MAT files found: %d\n', length(mat_files));
        for i = 1:length(mat_files)
            fprintf('  - %s (%.1f KB)\n', mat_files(i).name, mat_files(i).bytes/1024);
        end
        
        % Check for figures subdirectory
        figures_dir = fullfile(results_folder, 'figures');
        if exist(figures_dir, 'dir')
            png_files = dir(fullfile(figures_dir, '*.png'));
            fprintf('PNG files found: %d\n', length(png_files));
            for i = 1:length(png_files)
                fprintf('  - %s (%.1f KB)\n', png_files(i).name, png_files(i).bytes/1024);
            end
        else
            fprintf('❌ Figures directory not found: %s\n', figures_dir);
        end
    else
        fprintf('❌ Results folder does not exist: %s\n', results_folder);
    end
else
    fprintf('❌ latest_beamforming_folder not set\n');
end

% Check the specific error mentioned
fprintf('\n=== CHECKING FOR INDEX ERROR SOURCES ===\n');
problematic_vars = {'cbf_peak_pos', 'omp_peak_pos', 'sc_peak_pos', 'cbf_errors', 'omp_errors', 'sc_errors'};

for i = 1:length(problematic_vars)
    var_name = problematic_vars{i};
    if evalin('base', sprintf('exist(''%s'', ''var'')', var_name))
        var_value = evalin('base', var_name);
        fprintf('%s: size [%s], elements: %d\n', var_name, num2str(size(var_value)), numel(var_value));
        
        if size(var_value, 1) < 1
            fprintf('  ⚠️  WARNING: Empty first dimension!\n');
        end
        if length(var_value) < 1
            fprintf('  ⚠️  WARNING: Empty array!\n');
        end
    end
end

fprintf('\n=== DEBUG COMPLETE ===\n');
end




