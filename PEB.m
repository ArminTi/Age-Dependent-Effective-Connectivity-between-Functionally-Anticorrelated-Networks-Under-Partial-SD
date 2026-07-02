clear
clc

%% 1. Load FD Data 
fd_path = 'E:\FMRI_Team\New folder';
load('GCMs_modified.mat') 
FD_young_data = readmatrix(fullfile(fd_path, 'young_FD2.xlsx'));
FD_old_data   = readmatrix(fullfile(fd_path, 'old_FD2.xlsx'));

%% 2. Setup GCMs
GCM_young = [GCMs.GCM_young_full_sleep_DMN_DAN_SN; ...
             GCMs.GCM_young_sleep_deprived_DMN_DAN_SN];

GCM_old   = [GCMs.GCM_old_full_sleep_DMN_DAN_SN; ...
             GCMs.GCM_old_sleep_deprived_DMN_DAN_SN];

numYoungSubjects = length(GCMs.GCM_young_sleep_deprived_DMN_DAN_SN);
numOldSubjects   = length(GCMs.GCM_old_sleep_deprived_DMN_DAN_SN);

field = {'A'};

%% 3. PEB Analysis for Young Group 
fd_young_cov = [FD_young_data(:, 2); FD_young_data(:, 1)];
 
X_young = [ones(numYoungSubjects,1), -ones(numYoungSubjects,1), fd_young_cov(1:numYoungSubjects); ...
           ones(numYoungSubjects,1),  ones(numYoungSubjects,1), fd_young_cov(numYoungSubjects+1:end)];
 
X_young(:,2) = X_young(:,2) - mean(X_young(:,2));
X_young(:,3) = X_young(:,3) - mean(X_young(:,3)); 

M_young = struct('X', X_young, ...
                 'Xnames', {{'MeanEffect', 'SleepDeprivationEffect', 'FDEffect'}});
PEB_young = spm_dcm_peb(GCM_young, M_young, field);
BMA_young = spm_dcm_peb_bmc(PEB_young);
save('BMA_young.mat', 'BMA_young');

%% 4. PEB Analysis for Old Group 
fd_old_cov = [FD_old_data(:, 2); FD_old_data(:, 1)];

X_old = [ones(numOldSubjects,1), -ones(numOldSubjects,1), fd_old_cov(1:numOldSubjects); ...
         ones(numOldSubjects,1),  ones(numOldSubjects,1), fd_old_cov(numOldSubjects+1:end)];

% Mean-centering  
X_old(:,2) = X_old(:,2) - mean(X_old(:,2));
X_old(:,3) = X_old(:,3) - mean(X_old(:,3)); 

M_old = struct('X', X_old, ...
               'Xnames', {{'MeanEffect', 'SleepDeprivationEffect', 'FDEffect'}});
PEB_old = spm_dcm_peb(GCM_old, M_old, field);
%%
BMA_old = spm_dcm_peb_bmc(PEB_old);
save('BMA_old.mat', 'BMA_old');

%% 5. 3rd Level Analysis (Hierarchical PEB)
PEBs = {PEB_young; PEB_old};

X3 = [1, -1; 
      1,  1];

M3 = struct('X', X3, ...
            'Xnames', {{'Overall Mean', 'Age Group Difference'}});

PEB_hier = spm_dcm_peb(PEBs, M3, field);
BMA_hier = spm_dcm_peb_bmc(PEB_hier);
 
spm_dcm_peb_review(BMA_hier);


% -------------------------------
%  Author: Vida Feizi
%  Affiliation: M.Sc. Student, Computer Engineering, University of Tabriz
%  Email: vFeizii75@gmail.com
% -------------------------------


