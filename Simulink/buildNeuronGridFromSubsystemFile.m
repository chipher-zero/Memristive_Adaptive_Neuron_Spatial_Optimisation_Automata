%% Video-Ready Neuron Grid Generator (Floating Edges)
clear; clc;
modelName = 'Neuron_Grid_Model';

% --- 1. SETUP ---
% Source File Path
sourcePath = '/home/chipher/Desktop/GitHub_Projects/CBRAM_Modelling/Neuron_Subsystem.slx';
[sourceDir, sourceFile, ~] = fileparts(sourcePath);
addpath(sourceDir);

% Load the source (invisible)
if ~bdIsLoaded(sourceFile), load_system(sourcePath); end

% TARGET BLOCK:
% Ensure you have converted your file content to a Subsystem named "Neuron"
sourceBlock = [sourceFile, '/Neuron']; 

% Grid Configuration
N = 20; % Rows
M = 20; % Columns
spacing = [250, 200];

% Initialize Map
map = zeros(N, M); 
map1 = map;
% Reset/Create Model
if bdIsLoaded(modelName); close_system(modelName, 0); end
new_system(modelName);
open_system(modelName);

fprintf('Building %dx%d Grid (Floating Edges)...\n', N, M);

% --- 2. BUILD NODES ---
for r = 1:N
    for c = 1:M
        % Calculate Linear Index (k)
        k = (r-1)*M + c;
        map(r, c) = k;
        
        % Define Name & Path
        neuronName = sprintf('Neuron_%d', k);
        fullPath = [modelName, '/', neuronName];
        
        % Position
        x = c * spacing(1);
        y = r * spacing(2);
        
        % Add Block
        try
            add_block(sourceBlock, fullPath, 'Position', [x, y, x+80, y+80]);
        catch
            error('Could not find source "%s". Ensure you converted the file content to a Subsystem named "Neuron".', sourceBlock);
        end
        
        % --- CONSTANT INPUTS (In5 & In6) ---
        % In5: Synapse Resistance (Always connected)
        rsynName = sprintf('RSyn_%d', k);
        add_block('simulink/Sources/Constant', [modelName, '/', rsynName], ...
            'Value', '1e06', 'Position', [x-60, y+40, x-30, y+55]);
        add_line(modelName, [rsynName, '/1'], [neuronName, '/5'], 'autorouting', 'on');
        
        % In6: External Input (Always connected, default 0)
        extName = sprintf('Ext_%d', k);
        add_block('simulink/Sources/Constant', [modelName, '/', extName], ...
            'Value', '0', 'Position', [x-60, y+65, x-30, y+80]);
        add_line(modelName, [extName, '/1'], [neuronName, '/6'], 'autorouting', 'on');
        
        % Enable Logging (Out2)
        ph = get_param(fullPath, 'PortHandles');
        set_param(ph.Outport(2), 'DataLogging', 'on');
    end
end

% --- 3. CONNECTIVITY (Floating Edges) ---
for r = 1:N
    for c = 1:M
        currentK = (r-1)*M + c;
        currName = sprintf('Neuron_%d', currentK);
        
        % Define Neighbors [RowOffset, ColOffset, InputPort]
        neighbors = [-1, 0, 1;   % North -> In1
                      1, 0, 2;   % South -> In2
                      0, 1, 3;   % East  -> In3
                      0, -1, 4]; % West  -> In4
                      
        for i = 1:size(neighbors, 1)
            nr = r + neighbors(i, 1);
            nc = c + neighbors(i, 2);
            portIdx = neighbors(i, 3);
            
            % Check if neighbor exists
            if nr >= 1 && nr <= N && nc >= 1 && nc <= M
                % Calculate neighbor's linear index
                neighborK = (nr-1)*M + nc;
                neighborName = sprintf('Neuron_%d', neighborK);
                
                % Connect Neighbor Out2 -> Current In(PortIdx)
                add_line(modelName, [neighborName, '/2'], [currName, '/', num2str(portIdx)], 'autorouting', 'on');
            else
                % EDGE CASE: Neighbor does not exist.
                % PREVIOUSLY: We added a Ground block here.
                % NOW: We do nothing. The port remains unconnected (Floating).
            end
        end
    end
end

% --- 4. SAVE ---
save_system(modelName);
save('Neuron_Map.mat', 'map', 'N', 'M');

fprintf('Grid Complete.\n1. Model Saved: %s\n2. Map Saved: Neuron_Map.mat\n', modelName);