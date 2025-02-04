clear;
clc;

%% General Configuration
TR = 2.5; % Repetition time
TE = 0.03; % Echo time
rootDataPath = "C:\Users\DELL\Desktop\Data_Analysis"; % Root directory

% Define ROI groups
ROIs_group_1 = { 'networks.DefaultMode.MPFC (1,55,-3)' , 'networks.DefaultMode.LP (L) (-39,-77,33)', 'networks.DefaultMode.LP (R) (47,-67,29)', 'networks.DefaultMode.PCC (1,-61,38)' };  % DMN
ROIs_group_4 = {'yeo.Ventral Attention'}; % VAN

% allROIs = {ROIs_group_1,ROIs_group_2, ROIs_group_3, ROIs_group_4};
%group_name = {'DMN','DAN', 'WMN', 'VAN'};
allROIs = {ROIs_group_1, ROIs_group_4};
group_name = {'DMN', 'VAN'};

% Initialize GCM container as a structure
GCMs = struct();
groups = {'old', 'young'};
sessions = {'sleep_deprived', 'full_sleep'};

% Initialize GCM structure dynamically
for roiGroupIdx = 1:length(allROIs)
    roiGroupName = group_name{roiGroupIdx};
    for groupIdx = 1:2
        for sessionIdx = 1:2
            groupName = groups{groupIdx};
            sessionName = sessions{sessionIdx};
            fieldName = sprintf('GCM_%s_%s_%s', groupName, sessionName, roiGroupName);
            GCMs.(fieldName) = {}; % Initialize as empty cell
        end
    end
end

%% Process Subjects and Estimate DCMs
for groupIdx = 1:2
    groupFolder = groups{groupIdx};
    groupPath = fullfile(rootDataPath, groupFolder);
    
    subjectFiles = dir(fullfile(groupPath, 'ROI_Subject*_Session*.mat'));

    for roiGroupIdx = 1:length(allROIs)
        desiredROIs = allROIs{roiGroupIdx};
        roiGroupName = group_name{roiGroupIdx};
        fprintf("Processing ROIs Group %d (%s)...\n", roiGroupIdx, roiGroupName);
        
        for i = 1:length(subjectFiles)
            [~, subjName, ~] = fileparts(subjectFiles(i).name);
            subjInfo = regexp(subjName, 'ROI_Subject(\d+)_Session(\d+)', 'tokens');
            if isempty(subjInfo)
                warning('File name format incorrect: %s', subjectFiles(i).name);
                continue;
            end
            subj = str2double(subjInfo{1}{1}); 
            session = str2double(subjInfo{1}{2});
            
            fprintf("Processing Subject %d, Session %d...\n", subj, session);

            fileName = fullfile(groupPath, subjectFiles(i).name);
            if ~isfile(fileName)
                warning("File %s not found.", fileName);
                continue;
            end
            load(fileName);
            if ~exist('names', 'var') || ~exist('data', 'var')
                warning('Required variables not found in file: %s', fileName);
                continue;
            end

            % Extract ROI time series
            TimeSeries = [];
            ROI_name = {};  

            for roiIdx = 1:length(desiredROIs)
                roiIndex = find(strcmp(names, desiredROIs{roiIdx}));
                if ~isempty(roiIndex)
                    TimeSeries = [TimeSeries, data{1, roiIndex}];
                    ROI_name{end + 1} = names{roiIndex};
                else
                    warning("ROI %s not found in subject data.", desiredROIs{roiIdx});
                end
            end

            % Create and configure DCM
            DCM = configureDCM(TimeSeries, ROI_name, TR, TE, length(desiredROIs));

            % Estimate DCM
            try
                DCM = spm_dcm_fmri_csd(DCM);
                fprintf("DCM estimation completed for Subject %d, Session %d.\n", subj, session);

                saveFileName = fullfile(groupPath, sprintf("DCM_Subject%03d_%s_%s_%s.mat", subj, roiGroupName, groupFolder, sessions{session}));
                save(saveFileName, 'DCM');

                % Store in appropriate GCM dynamically
                fieldName = sprintf('GCM_%s_%s_%s', groupFolder, sessions{session}, roiGroupName);
                GCMs.(fieldName){end+1, 1} = DCM;

            catch ME
                fprintf("Error in DCM estimation for Subject %d, Session %d: %s\n", subj, session, ME.message);
            end
        end
    end
end

%% Save GCMs
save(fullfile(rootDataPath, "GCMs.mat"), 'GCMs');
fprintf("DCM estimation and saving completed.\n");

%% Helper Function: Configure DCM
function DCM = configureDCM(TimeSeries, ROI_name, TR, TE, n)
    v = size(TimeSeries, 1); % Number of time points

    DCM.v = v;
    DCM.n = n;
    DCM.Y.dt = TR;
    DCM.Y.Q = spm_Ce(ones(1, n) * v); 
    DCM.Y.y = TimeSeries;
    DCM.Y.name = ROI_name;

    DCM.U.u = zeros(v, 1);
    DCM.U.name = {'null'};

    DCM.a = ones(n, n); 
    DCM.b = zeros(n, n, 0); 
    DCM.c = zeros(n, 0); 
    DCM.d = zeros(n, n, 0); 

    DCM.TE = TE; 
    DCM.delays = repmat(TR, n, 1); 

    DCM.model = 'spm_csd'; 
    DCM.analysis = 'CSD'; 
    DCM.options.analysis = 'CSD'; 
    DCM.options.spatial = 0;
end

% -------------------------------
%  Author: Vida Feizi
%  Affiliation: M.Sc. Student, Computer Engineering, University of Tabriz
%  Email: vFeizii75@gmail.com
% -------------------------------