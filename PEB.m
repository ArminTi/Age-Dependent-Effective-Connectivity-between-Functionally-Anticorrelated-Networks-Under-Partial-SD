clear
clc
%% Loading

load('GCMs_modified.mat')
load("oldSleepDeprivedNorm.mat")
load("oldFullSleepNorm.mat")
load("youngSleepDeprivedNorm.mat")
load("youngFullSleepNorm.mat")

GCM_young = [GCMs.GCM_young_full_sleep_DMN_DAN_SN; GCMs.GCM_young_sleep_deprived_DMN_DAN_SN];
GCM_old = [GCMs.GCM_old_full_sleep_DMN_DAN_SN; GCMs.GCM_old_sleep_deprived_DMN_DAN_SN];
numYoungSubjects = length(GCMs.GCM_young_sleep_deprived_DMN_DAN_SN);
numOldSubjects = length(GCMs.GCM_old_sleep_deprived_DMN_DAN_SN);

%% Peb for Young

X_young = [ones(numYoungSubjects, 1), -ones(numYoungSubjects, 1), youngFullSleepNorm; ...
                     ones(numYoungSubjects, 1), ones(numYoungSubjects, 1), youngSleepDeprivedNorm];
X_labels_young_first = {'MeanEffect', 'SleepDeprivationEffect', 'KSSScore'};
M_young = struct();
M_young.X = X_young;
M_young.Xnames = X_labels_young_first;
field = {'A'};
PEB_young = spm_dcm_peb(GCM_young, M_young, field);
BMA_young = spm_dcm_peb_bmc(PEB_young);
save("BMA_young", "BMA_young")

%Peb for old

X_old = [ones(numOldSubjects, 1), -ones(numOldSubjects, 1), oldFullSleepNorm; ...
                     ones(numOldSubjects, 1), ones(numOldSubjects, 1), oldSleepDeprivedNorm];
X_labels_old_first = {'MeanEffect', 'SleepDeprivationEffect', 'KSSScore'};
M_old = struct();
M_old.X = X_old;
M_old.Xnames = X_labels_young_first;
field = {'A'};
PEB_old = spm_dcm_peb(GCM_old, M_old, field);
BMA_old = spm_dcm_peb_bmc(PEB_old);
save("BMA_old", "BMA_old")

%% Review for Young
spm_dcm_peb_review(BMA_young)

%% Review for Old
spm_dcm_peb_review(BMA_old)

% -------------------------------
%  Author: Vida Feizi
%  Affiliation: M.Sc. Student, Computer Engineering, University of Tabriz
%  Email: vFeizii75@gmail.com
% -------------------------------


