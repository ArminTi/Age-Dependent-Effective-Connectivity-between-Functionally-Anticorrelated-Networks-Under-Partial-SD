clear;
clc;

%% Load Pre-Estimated GCM Structure
rootDataPath = "C:\Users\DELL\Desktop\Data_Analysis";
load(fullfile(rootDataPath, "GCMs.mat"), 'GCMs');

% Define groups, conditions, and ROI groups
groups = {'old', 'young'};
sessions = {'sleep_deprived', 'full_sleep'};
roiGroups = {'DMN'}; % Add more ROIs if needed

% Initialize a struct to dynamically store all GCMs
GCM_struct = struct();

for roiIdx = 1:length(roiGroups)
    roiGroup = roiGroups{roiIdx};
    
    for groupIdx = 1:length(groups)
        groupName = groups{groupIdx};
        
        % Initialize a variable to combine GCMs across sessions
        combined_GCM = [];
        
        for sessionIdx = 1:length(sessions)
            sessionName = sessions{sessionIdx};
            
            % Construct field name to extract from GCMs
            fieldName = sprintf('GCM_%s_%s_%s', groupName, sessionName, roiGroup);
            
            % Check if the field exists in the GCMs structure
            if isfield(GCMs, fieldName)
                % Concatenate GCMs
                combined_GCM = [combined_GCM; GCMs.(fieldName)];
            end
        end
        
        % Store in the new structure
        structField = sprintf('GCM_%s_%s', groupName, roiGroup);
        GCM_struct.(structField) = combined_GCM;
    end
end

%% Perform PEB Analysis for All ROIs

PEB_struct = struct();
BMA_struct = struct();

for roiIdx = 1:length(roiGroups)
    roiGroup = roiGroups{roiIdx};

    % **First-Level Analysis: Young Group**
    num_young = length(GCM_struct.(sprintf('GCM_young_%s', roiGroup)));
    
    % Design matrix for within-subject factor (sleep deprivation)
    X_young_first = [ones(num_young/2, 1), ones(num_young/2, 1);
                     ones(num_young/2, 1), -ones(num_young/2, 1)];
    X_labels_young_first = {'MeanEffect', 'SleepDeprivationEffect'};

    % Model settings
    M_young_first = struct();
    M_young_first.Q = 'all';
    M_young_first.X = X_young_first;
    M_young_first.Xnames = X_labels_young_first;

    field = {'A'}; % Specify which parameters to analyze
    PEB_young_first = spm_dcm_peb(GCM_struct.(sprintf('GCM_young_%s', roiGroup)), M_young_first, field);
    PEB_struct.(sprintf('PEB_young_%s', roiGroup)) = PEB_young_first;
    BMA_result_young = spm_dcm_peb_bmc(PEB_young_first);
    
    % Review results for Young Group and wait for Enter to continue
    spm_dcm_peb_review(BMA_result_young); 
    input('Press Enter to continue to the next result...');
    close(gcf); % Close review window automatically

    % **First-Level Analysis: Old Group**
    num_old = length(GCM_struct.(sprintf('GCM_old_%s', roiGroup)));
    
    % Design matrix for within-subject factor (sleep deprivation)
    X_old_first = [ones(num_old/2, 1), ones(num_old/2, 1);
                   ones(num_old/2, 1), -ones(num_old/2, 1)];
    X_labels_old_first = {'MeanEffect', 'SleepDeprivationEffect'};

    % Model settings
    M_old_first = struct();
    M_old_first.Q = 'all';
    M_old_first.X = X_old_first;
    M_old_first.Xnames = X_labels_old_first;

    PEB_old_first = spm_dcm_peb(GCM_struct.(sprintf('GCM_old_%s', roiGroup)), M_old_first, field);
    PEB_struct.(sprintf('PEB_old_%s', roiGroup)) = PEB_old_first;
    BMA_result_old = spm_dcm_peb_bmc(PEB_old_first);
    
    % Review results for Old Group and wait for Enter to continue
    spm_dcm_peb_review(BMA_result_old);
    input('Press Enter to continue to the next result...');
    close(gcf); % Close review window automatically

    % **Second-Level Analysis: Young vs. Old Comparison**
    PEB_combined = {PEB_young_first; PEB_old_first};

    % Design matrix for between-subject factor (age group)
    X_second = [1, 1;
                1, -1];
    X_labels_second = {'GroupMean', 'AgeDifference'};

    % Model settings
    M_second = struct();
    M_second.X = X_second;
    M_second.Xnames = X_labels_second;

    % Perform second-level PEB
    PEB_second = spm_dcm_peb(PEB_combined, M_second, field);
    PEB_struct.(sprintf('PEB_second_%s', roiGroup)) = PEB_second;

    % **Bayesian Model Averaging (BMA)**
    BMA_result = spm_dcm_peb_bmc(PEB_second);
    BMA_struct.(sprintf('BMA_%s', roiGroup)) = BMA_result;

    % Review results for Second-level analysis and wait for Enter to continue
    spm_dcm_peb_review(BMA_result);
    input('Press Enter to continue...');
    close(gcf); % Close review window automatically
end

% -------------------------------
%  Author: Vida Feizi
%  Affiliation: M.Sc. Student, Computer Engineering, University of Tabriz
%  Email: vFeizii75@gmail.com
% -------------------------------

