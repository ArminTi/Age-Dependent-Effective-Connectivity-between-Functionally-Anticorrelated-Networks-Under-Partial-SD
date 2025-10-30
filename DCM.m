%% General Configuration
TR = 2.5; % Repetition time
TE = 0.03; % Echo time
rootDataPath = "C:\Users\Vida\Desktop\aida"; 
TIME = tic;
% Define ROIs (only one group)
desiredROIs = {'networks.DefaultMode.MPFC (1,55,-3)', 'networks.DefaultMode.LP (L) (-39,-77,33)', ...
    'networks.DefaultMode.LP (R) (47,-67,29)', 'networks.DefaultMode.PCC (1,-61,38)', ...
     'networks.DorsalAttention.FEF (L)  (-27,-9,64)', 'networks.DorsalAttention.FEF (R)  (30,-6,64)', ...
   'networks.DorsalAttention.IPS (L)  (-39,-43,52)', 'networks.DorsalAttention.IPS (R)  (39,-42,54)',...
  'networks.Salience.ACC (0,22,35)', 
     'networks.Salience.AInsula (L) (-44,13,1)', 
     'networks.Salience.AInsula (R) (47,14,0)', 
     'networks.Salience.RPFC (L) (-32,45,27)', 
     'networks.Salience.RPFC (R) (32,46,27)', 
     'networks.Salience.SMG (L) (-60,-39,31)', 
     'networks.Salience.SMG (R) (62,-35,32)'
    };
roiGroupName = 'DMN_DAN_SN'; % Combined name since it includes both DMN and DAN ROIs

% Define possible groups and sessions
possibleGroups = {'young', 'old'};
sessions = {'sleep_deprived', 'full_sleep'};

% Determine existing groups dynamically
existingGroups = {};
for g = 1:length(possibleGroups)
    groupFolder = possibleGroups{g};
    groupPath = fullfile(rootDataPath, groupFolder);
    if isfolder(groupPath) && ~isempty(dir(fullfile(groupPath, 'ROI_Subject*_Session*.mat')))
        existingGroups{end+1} = groupFolder;
    end
end

% Initialize GCM container as a structure only for existing groups
GCMs = struct();
for groupIdx = 1:length(existingGroups)
    for sessionIdx = 1:length(sessions)
        groupName = existingGroups{groupIdx};
        sessionName = sessions{sessionIdx};
        fieldName = sprintf('GCM_%s_%s_%s', groupName, sessionName, roiGroupName);
        GCMs.(fieldName) = {}; % Initialize as empty cell
    end
end

%% Process Data in One Pass
for groupIdx = 1:length(existingGroups)
    groupFolder = existingGroups{groupIdx};
    groupPath = fullfile(rootDataPath, groupFolder);
    
    subjectFiles = dir(fullfile(groupPath, 'ROI_Subject*_Session*.mat'));
    if isempty(subjectFiles)
        continue; % This should not happen due to earlier check, but kept for safety
    end
    
    fprintf("Processing ROIs Group (%s) for %s...\n", roiGroupName, groupFolder);
    
    for i = 1:length(subjectFiles)
        fileName = fullfile(groupPath, subjectFiles(i).name);
        subjInfo = split(subjectFiles(i).name, {'Subject', '_Session', '.mat'});
        subj = str2double(subjInfo{2});
        session = str2double(subjInfo{3});
        
        if isnan(subj) || isnan(session)
            fprintf('Invalid file name format: %s\n', subjectFiles(i).name);
            continue;
        end
        
        fprintf("Processing Subject %d, Session %d...\n", subj, session);
        
        % Determine session name
        sessionName = sessions{session};
        
        % Check for existing DCM
        saveFileName = fullfile(groupPath, sprintf("DCM_Subject%03d_%s_%s_%s.mat", ...
            subj, roiGroupName, groupFolder, sessionName));
        if isfile(saveFileName)
            fprintf("Loading existing DCM for Subject %d, Session %d...\n", subj, session);
            load(saveFileName, 'DCM');
            fieldName = sprintf('GCM_%s_%s_%s', groupFolder, sessionName, roiGroupName);
            GCMs.(fieldName){end+1, 1} = DCM;
            continue;
        end
        
        % Load and process data
        if ~isfile(fileName)
            fprintf("File %s not found.\n", fileName);
            continue;
        end
        dataStruct = load(fileName);
        if ~isfield(dataStruct, 'names') || ~isfield(dataStruct, 'data')
            fprintf('Required variables not found in file: %s\n', fileName);
            continue;
        end
        
        % Extract ROI time series
        [TimeSeries, ROI_name] = extractROIs(dataStruct.names, dataStruct.data, desiredROIs);
        if isempty(TimeSeries)
            fprintf("No matching ROIs found for Subject %d, Session %d.\n", subj, session);
            continue;
        end
        
        % Configure and estimate DCM
        DCM = configureDCM(TimeSeries, ROI_name, TR, TE, length(ROI_name));
        DCM = spm_dcm_fmri_csd(DCM);
        fprintf("DCM estimation completed for Subject %d, Session %d.\n", subj, session);
        
        % Save DCM
        save(saveFileName, 'DCM');
        
        % Store in appropriate GCM dynamically
        fieldName = sprintf('GCM_%s_%s_%s', groupFolder, sessionName, roiGroupName);
        GCMs.(fieldName){end+1, 1} = DCM;
    end
end

%% Save GCMs
if ~isempty(fieldnames(GCMs))
    save(fullfile(rootDataPath, "GCMs.mat"), 'GCMs');
    fprintf("DCM estimation and saving completed as 'GCMs.mat'.\n");
else
    fprintf("No data processed. GCMs is empty.\n");
end

%% Helper Function: Extract ROIs
function [TimeSeries, ROI_name] = extractROIs(names, data, desiredROIs)
    TimeSeries = [];
    ROI_name = {};
    for roiIdx = 1:length(desiredROIs)
        roiIndex = find(strcmp(names, desiredROIs{roiIdx}));
        if ~isempty(roiIndex)
            TimeSeries = [TimeSeries, data{1, roiIndex}];
            ROI_name{end + 1} = names{roiIndex};
        end
    end
end

%% Helper Function: Configure DCM
function DCM = configureDCM(TimeSeries, ROI_name, TR, TE, n)
    v = size(TimeSeries, 1);
    DCM.v = v;
    DCM.n = n;
    DCM.Y.dt = TR;
    DCM.Y.Q = spm_Ce(ones(1, n) * v);
    DCM.Y.y = TimeSeries;
    DCM.Y.name = ROI_name;
    DCM.U.u = zeros(v, 1);
    DCM.U.name = {'null'};
    DCM.a = ones(n, n); % Full connectivity
    DCM.b = zeros(n, n, 0);
    DCM.c = zeros(n, 0);
    DCM.d = zeros(n, n, 0);
    DCM.TE = TE;
    DCM.delays = repmat(TR, n, 1);
    DCM.model = 'spm_csd';
    DCM.options.nonlinear  = 0;
    DCM.options.two_state  = 0;
    DCM.options.stochastic = 0;
    DCM.options.analysis   = 'CSD';
    DCM.options.induced   = 1;
end

% -------------------------------
%  Author: Vida Feizi
%  Affiliation: M.Sc. Student, Computer Engineering, University of Tabriz
%  Email: vFeizii75@gmail.com

% -------------------------------
