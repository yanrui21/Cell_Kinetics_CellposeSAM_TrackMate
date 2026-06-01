%% TrackMate Trajectory Analysis: ROI, Splits, and MSD
% -------------------------------------------------------------------------
% Description:
%   1. Loads Spots, Edges, and ROI files.
%   2. Filters for tracks that START within the ROI.
%   3. Analyzes Splitting (Frequency per track).
%   4. Decomposes trees into linear "Root-to-Leaf" paths.
%   5. Calculates Kinematics (Velocity, Persistence) and MSD.
% -------------------------------------------------------------------------

clear all; 

%% 1. USER SETTINGS
% =========================================================================
% File Paths (Change these to your actual files)
spotsFile = "1_Image_spots.csv";
edgesFile = "1_Image_edges.csv";
roiFile   = "1_left_epi_roi.csv";

% Physical Units
pixelSize = 0.36;     % microns (or pixels if set to 1)
frameInterval = 5; % minutes (or frames if set to 1)
xDirection = 1;      % 1 if positive X is "forward", -1 if negative X is forward
maxTime = 128; % Define the cutoff time in frame

%% 2. DATA IMPORT
% =========================================================================
disp('Importing data...');

% A. Load ROI
optsRoi = detectImportOptions(roiFile);
roiTable = readtable(roiFile, optsRoi);
roiPolyX = roiTable.X;
roiPolyY = roiTable.Y;

% B. Load Spots
% TrackMate CSVs have 3 lines of metadata; header is on line 4.
optsS = detectImportOptions(spotsFile);
optsS.DataLines = [4 Inf]; 
spots = readtable(spotsFile, optsS);

% Rename columns for easier access and handle different versions
nameMap = {'ID','SpotID'; 'TRACK_ID','TrackID'; 'POSITION_X','X'; ...
           'POSITION_Y','Y'; 'POSITION_T','T'; 'FRAME','Frame'};
for i = 1:size(nameMap,1)
    if ismember(nameMap{i,1}, spots.Properties.VariableNames)
        spots = renamevars(spots, nameMap{i,1}, nameMap{i,2});
    end
end

% Remove invalid rows (sometimes imports empty lines)
spots = spots(~isnan(spots.SpotID), :);

% C. Load Edges
optsE = detectImportOptions(edgesFile);
optsE.DataLines = [4 Inf];
edges = readtable(edgesFile, optsE);
edges = renamevars(edges, {'SPOT_SOURCE_ID', 'SPOT_TARGET_ID', 'TRACK_ID'}, ...
                          {'SourceID', 'TargetID', 'TrackID'});
edges = edges(~isnan(edges.SourceID), :);

%% 2.5 TIME FILTERING (Limit analysis to first x mins)
% =========================================================================

disp(['Filtering data to first ' num2str(maxTime) ' time units...']);

% Filter Spots
spots = spots(spots.T <= maxTime, :);

% Filter Edges (Keep only edges where both Source and Target exist in filtered spots)
validSpotIDs = spots.SpotID;
validSource = ismember(edges.SourceID, validSpotIDs);
validTarget = ismember(edges.TargetID, validSpotIDs);
edges = edges(validSource & validTarget, :);

disp(['Remaining Spots: ' num2str(height(spots))]);
disp(['Remaining Edges: ' num2str(height(edges))]);

%% 3. BUILD GRAPH
% =========================================================================
disp('Building trajectory graph...');

% Convert IDs to strings to ensure digraph treats them as Names, not Indices
G = digraph(string(edges.SourceID), string(edges.TargetID));

% We need to attach X, Y, T data to the graph nodes.
% 1. Create a map from SpotID -> Table Index
spotMap = containers.Map(spots.SpotID, 1:height(spots));

% 2. Get Node Names (Spot IDs) from graph
nodeNames = G.Nodes.Name;
nodeIDs = cellfun(@str2double, nodeNames);

% 3. Pre-allocate Node Data
numNodes = numnodes(G);
nodeX = nan(numNodes, 1);
nodeY = nan(numNodes, 1);
nodeT = nan(numNodes, 1);

% 4. Fill Data
for i = 1:numNodes
    sid = nodeIDs(i);
    if isKey(spotMap, sid)
        idx = spotMap(sid);
        nodeX(i) = spots.X(idx) * pixelSize;
        nodeY(i) = spots.Y(idx) * pixelSize;
        nodeT(i) = spots.T(idx) * frameInterval;
    end
end

% Add to Graph Table
G.Nodes.X = nodeX;
G.Nodes.Y = nodeY;
G.Nodes.T = nodeT;

%% 4. ROI FILTERING
% =========================================================================
disp('Filtering tracks starting in ROI...');

uniqueTracks = unique(spots.TrackID);
selectedTracks = []; % List of TrackIDs

% Plot Setup
figure('Color','w', 'Position', [50, 50, 600, 600]); hold on; axis equal;
set(gca, 'YDir', 'reverse'); % Images usually have Y down
plot(roiPolyX, roiPolyY, 'k-', 'LineWidth', 2);
title('ROI and Selected Start Points');
xlabel('X'); ylabel('Y');

for i = 1:length(uniqueTracks)
    tid = uniqueTracks(i);
    
    % Get all spots for this track
    trackMask = (spots.TrackID == tid);
    trackSpots = spots(trackMask, :);
    
    % Find Start (Minimum Time)
    [~, minIdx] = min(trackSpots.T);
    startX = trackSpots.X(minIdx);
    startY = trackSpots.Y(minIdx);
    
    % Check ROI
    if inpolygon(startX, startY, roiPolyX, roiPolyY)
        selectedTracks(end+1) = tid; 
        plot(startX, startY, 'b.', 'MarkerSize', 10);
    end
end

fprintf('Selected %d out of %d total tracks.\n', length(selectedTracks), length(uniqueTracks));

%% 5. ANALYSIS: SPLITS & DECOMPOSITION
% =========================================================================
disp('Analyzing tracks...');

splitStats = struct('TrackID', {}, 'Splits', {}, 'Duration', {}, 'Rate', {});
allLinearPaths = {}; % Store [T, X, Y] for every branch

for i = 1:length(selectedTracks)
    tid = selectedTracks(i);
    
    % 1. Extract Spot IDs for this track
    trackSpotIDs = spots.SpotID(spots.TrackID == tid);
    
    % 2. Find which of these spots are actually in the Graph
    % (Convert to string to match G.Nodes.Name)
    targetNames = string(trackSpotIDs);
    validNodeNames = intersect(targetNames, G.Nodes.Name);
    
    if isempty(validNodeNames), continue; end
    
    % 3. Create Subgraph
    % This inherits X, Y, T from the main graph G
    subG = subgraph(G, validNodeNames);
    
    % 4. Analyze Splits
    % A split is any node with Out-Degree > 1
    nSplits = sum(outdegree(subG) > 1);
    
    % Duration (Max T - Min T in this track)
    % We access the 'T' variable from the Nodes table of the subgraph
    if ismember('T', subG.Nodes.Properties.VariableNames)
        ts = subG.Nodes.T;
        dur = max(ts) - min(ts);
    else
        % Fallback if T didn't transfer for some reason
        dur = NaN; 
    end
    
    % Avoid division by zero
    if dur > 0
        rate = nSplits / dur;
    else
        rate = 0;
    end
    
    splitStats(end+1) = struct('TrackID', tid, 'Splits', nSplits, ...
                               'Duration', dur, 'Rate', rate); 
    
    % 5. Decompose into Linear Paths (Root -> Leaf)
    roots = find(indegree(subG) == 0);
    leaves = find(outdegree(subG) == 0);
    
    % Iterate all Root-Leaf combinations
    for r = 1:length(roots)
        for l = 1:length(leaves)
            % Find shortest path from specific root to specific leaf
            try
                pathIdx = shortestpath(subG, roots(r), leaves(l));
            catch
                continue;
            end
            
            if ~isempty(pathIdx) && length(pathIdx) > 2
                % Extract kinematics from subgraph nodes
                pT = subG.Nodes.T(pathIdx);
                pX = subG.Nodes.X(pathIdx);
                pY = subG.Nodes.Y(pathIdx);
                
                % Sort by time (crucial for valid MSD)
                [pT, sortI] = sort(pT);
                pX = pX(sortI);
                pY = pY(sortI);
                
                allLinearPaths{end+1} = [pT, pX, pY]; 
            end
        end
    end
end
% Assumes you have run the previous analysis script and have:
% - G (digraph)
% - selectedTracks (list of TrackIDs starting in ROI)
% - roiPolyX, roiPolyY

figure('Color','w'); hold on; axis equal;
set(gca, 'YDir', 'reverse'); % Image coordinates
colormap(jet); % or parula, viridis

% 1. Plot ROI
plot(roiPolyX.*pixelSize, roiPolyY.*pixelSize, 'k-', 'LineWidth', 2);

% 2. Prepare for Color Mapping
allTimes = G.Nodes.T;
minT = min(allTimes);
maxT = max(allTimes);

% 3. Iterate Tracks and Plot
for i = 1:length(selectedTracks)
    tid = selectedTracks(i);
    
    % --- CORRECTED SELECTION LOGIC ---
    % 1. Get the list of Spot IDs belonging to this track from the Table
    trackSpotIDs = spots.SpotID(spots.TrackID == tid);
    
    % 2. Convert to string (to match Graph format)
    targetNames = string(trackSpotIDs);
    
    % 3. Find which of these nodes exist in the Graph
    % (This prevents the array bounds error)
    validNodeNames = intersect(targetNames, string(G.Nodes.Name));
    
    if isempty(validNodeNames), continue; end
    
    % 4. Create the subgraph safely
    subG = subgraph(G, validNodeNames);
    % ---------------------------------
    
    % Get Edges
    edgesSub = subG.Edges;
    if isempty(edgesSub), continue; end
    
    % Plot each edge as a segment colored by Source Time
    for e = 1:height(edgesSub)
        sID = edgesSub.EndNodes{e, 1};
        tID = edgesSub.EndNodes{e, 2};
        
        % Lookup Coordinates in Main Graph (G)
        sIdx = findnode(G, sID);
        tIdx = findnode(G, tID);
        
        % Ensure nodes were found
        if sIdx == 0 || tIdx == 0, continue; end
        
        X = [G.Nodes.X(sIdx), G.Nodes.X(tIdx)];
        Y = [G.Nodes.Y(sIdx), G.Nodes.Y(tIdx)];
        T = G.Nodes.T(sIdx);
        
        % Plot Line Segment
        surface([X;X], [Y;Y], [0 0; 0 0], [T T; T T], ...
            'FaceColor','no', 'EdgeColor','interp', 'LineWidth', 1.5);
    end
    
    % 4. Highlight Splits
    splitNodes = find(outdegree(subG) > 1);
    if ~isempty(splitNodes)
        splitNodeNames = subG.Nodes.Name(splitNodes);
        
        % Map back to main graph for coordinates
        mainIndices = findnode(G, splitNodeNames);
        
        sx = G.Nodes.X(mainIndices);
        sy = G.Nodes.Y(mainIndices);
        
        scatter(sx, sy, 20, 'm', 'filled', 'MarkerEdgeColor', 'k');
    end
end

% Formatting
colorbar;
caxis([minT maxT]);
ylabel(colorbar, 'Time (min)');
title('Tracks in ROI: Time-Coded with Splits');
xlabel('X Position'); ylabel('Y Position');
%% 6. KINEMATICS & MSD (Principal Direction)
% =========================================================================
disp('Calculating Kinematics (Principal Direction) & MSD...');

% --- PASS 1: Calculate Global Principal Direction ---
totalDispX = 0;
totalDispY = 0;
totalDur   = 0;
trackCache = cell(length(selectedTracks), 1); % Cache to avoid recomputing graphs

for i = 1:length(selectedTracks)
    tid = selectedTracks(i);
    
    trackSpotIDs = spots.SpotID(spots.TrackID == tid);
    validNodeNames = intersect(string(trackSpotIDs), G.Nodes.Name);
    if isempty(validNodeNames), continue; end
    subG = subgraph(G, validNodeNames);
    
    roots = find(indegree(subG) == 0);
    leaves = find(outdegree(subG) == 0);
    if isempty(roots), continue; end
    
    rootNode = roots(1);
    
    sumDispX = 0;
    sumDispY = 0;
    sumDur = 0;
    thesePaths = {};
    
    for l = 1:length(leaves)
        leafNode = leaves(l);
        pathIdx = shortestpath(subG, rootNode, leafNode);
        if isempty(pathIdx), continue; end
        
        pT = subG.Nodes.T(pathIdx);
        pX = subG.Nodes.X(pathIdx);
        pY = subG.Nodes.Y(pathIdx);
        
        [pT, sortI] = sort(pT);
        pX = pX(sortI);
        pY = pY(sortI);
        
        % Calculate X and Y displacement for this branch
        dX = pX(end) - pX(1);
        dY = pY(end) - pY(1);
        dT = pT(end) - pT(1);
        
        if dT > 0
            sumDispX = sumDispX + dX;
            sumDispY = sumDispY + dY;
            sumDur   = sumDur + dT;
        end
        thesePaths{end+1} = [pT, pX, pY]; 
    end
    
    % Accumulate global totals
    totalDispX = totalDispX + sumDispX;
    totalDispY = totalDispY + sumDispY;
    totalDur   = totalDur + sumDur;
    
    % Save data to cache so we don't have to rebuild the graph in Pass 2
    trackCache{i} = struct('DispX', sumDispX, 'DispY', sumDispY, 'Dur', sumDur, 'Paths', {thesePaths});
end

% Compute the Unit Vector (uX, uY) for the Principal Direction
globalVx = totalDispX / totalDur;
globalVy = totalDispY / totalDur;
magV = sqrt(globalVx^2 + globalVy^2);
uX = globalVx / magV;
uY = globalVy / magV;

disp(['Principal Direction Vector (uX, uY): (', num2str(uX, '%.3f'), ', ', num2str(uY, '%.3f'), ')']);

% --- PASS 2: Project Velocities and Calculate MSD ---
kinematics = []; % [ProjectedVel, Persistence, Duration]
maxLag = 50; 
msdAccumulator = cell(maxLag, 1); 

for i = 1:length(selectedTracks)
    if isempty(trackCache{i}), continue; end
    tc = trackCache{i};
    
    % 1. Projected Velocity (Dot Product of Track Velocity and Unit Vector)
    if tc.Dur > 0
        trackVx = tc.DispX / tc.Dur;
        trackVy = tc.DispY / tc.Dur;
        projV = (trackVx * uX) + (trackVy * uY); 
    else
        projV = NaN;
    end
    
    % 2. Persistence Calculation
    totalP = 0; countP = 0;
    for p = 1:length(tc.Paths)
        P = tc.Paths{p};
        if size(P, 1) < 2, continue; end
        steps = diff(P(:,2:3));
        pathLen = sum(sqrt(sum(steps.^2, 2)));
        netDisp = sqrt((P(end,2)-P(1,2))^2 + (P(end,3)-P(1,3))^2);
        if pathLen > 0
            totalP = totalP + (netDisp / pathLen);
            countP = countP + 1;
        end
    end
    avgPersistence = totalP / max(countP, 1);
    
    kinematics(end+1, :) = [projV, avgPersistence, tc.Dur]; 
    
    % 3. MSD Calculation
    for p = 1:length(tc.Paths)
        P = tc.Paths{p};
        t_idx = round((P(:,1) - P(1,1)) / frameInterval) + 1;
        pos = P(:,2:3);
        
        for t1 = 1:length(t_idx)
            for t2 = (t1+1):length(t_idx)
                lag = t_idx(t2) - t_idx(t1);
                if lag <= maxLag && lag > 0
                    sqDisp = sum((pos(t2,:) - pos(t1,:)).^2);
                    msdAccumulator{lag} = [msdAccumulator{lag}; sqDisp];
                end
            end
        end
    end
end

% =========================================================================
% POST-PROCESSING MSD
% =========================================================================
msdMean = zeros(maxLag, 1);
msdSEM  = zeros(maxLag, 1);
lags    = (1:maxLag)' * frameInterval; 

for i = 1:maxLag
    data = msdAccumulator{i};
    if ~isempty(data)
        msdMean(i) = mean(data);
        msdSEM(i)  = std(data) / sqrt(length(data));
    else
        msdMean(i) = NaN;
    end
end

validFit = find(~isnan(msdMean) & (1:maxLag)' > 2); 
if length(validFit) > 4
    logT = log(lags(validFit));
    logM = log(msdMean(validFit));
    poly = polyfit(logT, logM, 1);
    alpha = poly(1);
    D_coeff = exp(poly(2))/4;
else
    alpha = NaN; D_coeff = NaN;
end

%% 7. VISUALIZATION & OUTPUT
% =========================================================================
figure('Name','Analysis Results','Color','w','Position',[100,100,1200,400]);


% A. Velocity (Principal Direction)
subplot(1, 3, 1);
histogram(kinematics(:,1), 10, 'FaceColor', [0.8500 0.3250 0.0980]); % Orange
xlabel('Velocity (Principal Dir)');
ylabel('Count');

meanProjV = mean(kinematics(:,1), 'omitnan');
title(['Mean V_{p}: ' num2str(meanProjV, '%.3f')]);

xline(0, 'r--', 'LineWidth', 1); % Zero reference
xline(meanProjV, 'b-', 'LineWidth', 2); % Mean line

% B. Persistence
subplot(1, 3, 2);
histogram(kinematics(:,2));
xlabel('Persistence (D/L)'); title(['Mean P: ' num2str(mean(kinematics(:,2)),'%.2f')]);
xlim([0 1]);

% C. MSD (Linear Scale)
subplot(1, 3, 3);
errorbar(lags, msdMean, msdSEM, 'ko', 'MarkerFaceColor','k');
% set(gca, 'XScale', 'log', 'YScale', 'log'); % <-- Commented out for linear scale
hold on;

if ~isnan(alpha)
    % The fitted curve will now visually appear as a curve rather than a straight line
    fitY = 4 * D_coeff * (lags .^ alpha);
    plot(lags, fitY, 'r-', 'LineWidth', 2);
    
    % Adjusting text position so it looks good on a linear plot
    textLocX = lags(round(end/2));
    textLocY = fitY(round(end/2)) * 1.2; 
    text(textLocX, textLocY, sprintf('\\alpha = %.2f', alpha), ...
        'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
end

xlabel('Time Lag (min)'); 
ylabel('MSD (\mum^2)');
title('Ensemble MSD (Linear)');
grid on;