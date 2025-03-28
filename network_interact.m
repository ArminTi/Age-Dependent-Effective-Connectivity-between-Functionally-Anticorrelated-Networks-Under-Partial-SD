%% Loading GCMs
clear all;
clc;

load('GCM_old.mat');
load('GCM_young.mat');
ndcm_old = 10;
ndcm_young = 10;
%% PEB for young

M = struct();
M.X = [ones(ndcm_old,1),[repmat(0,[ndcm_old/2,1]); repmat(+1,[ndcm_old/2,1])]];
field = {'A'};
PEB_young = spm_dcm_peb(GCM_young, M, field);
BMA_young = spm_dcm_peb_bmc(PEB_young);
save('BMA_young.mat', 'BMA_young');
%% PEB for old
PEB_old = spm_dcm_peb(GCM_old, M, field);
BMA_old = spm_dcm_peb_bmc(PEB_old);
save('BMA_old.mat', 'BMA_old');

%% Load BMA
clear all
clc
load('BMA_old.mat');
load('BMA_young.mat');


%% for old

Ep = BMA_old.Ep;  
Pp = BMA_old.Pp;  
threshold = 0.95; 

M_cov1 = reshape(Ep(1:64), [8,8]);
M_cov2 = reshape(Ep(65:128), [8,8]);
Pp_cov1 = reshape(Pp(1:64), [8,8]);
Pp_cov2 = reshape(Pp(65:128), [8,8]);

%x = M_cov1;
M_cov1(Pp_cov1 <= threshold) = 0;
M_cov2(Pp_cov2 <= threshold) = 0;


figure;
subplot(1,2,1);
imagesc(M_cov1); colorbar;
title('Covariate 1 Connectivity (Filtered by Pp > 0.95)');
set(gca, 'XTick',1:8, 'YTick',1:8);
subplot(1,2,2);
imagesc(M_cov2); colorbar;
title('Covariate 2 Connectivity (Filtered by Pp > 0.95)');
set(gca, 'XTick',1:8, 'YTick',1:8);

%% Average for old 

DMN = 1:4;   
DAN = 5:8;   
valid_out_cov1 = M_cov1(DMN, DAN); % DAN →  DMN connections
valid_in_cov1  = M_cov1(DAN, DMN); % DMN →  DAN connections


avg_out_cov1 = sum(valid_out_cov1(:)) / nnz(valid_out_cov1);
avg_in_cov1  = sum(valid_in_cov1(:)) / nnz(valid_in_cov1);
hierarchy_DMN_cov1_old = abs(avg_in_cov1) - abs(avg_out_cov1);
hierarchy_DAN_cov1_old = abs(avg_out_cov1) - abs(avg_in_cov1);


valid_out_cov2 = M_cov2(DMN, DAN); % DMN → DAN connections
valid_in_cov2  = M_cov2(DAN, DMN); % DAN → DMN connections


avg_out_cov2 = sum(valid_out_cov2(:)) / nnz(valid_out_cov2);
avg_in_cov2  = sum(valid_in_cov2(:)) / nnz(valid_in_cov2);
hierarchy_DMN_cov2_old = abs(avg_in_cov2) - abs(avg_out_cov2);
hierarchy_DAN_cov2_old = abs(avg_out_cov2) - abs(avg_in_cov2);

%% for young

clear Ep Pp M_cov1 P_cov2 avg_out_cov1 avg_in_cov1 avg_out_cov2 avg_in_cov2

Ep = BMA_young.Ep;  
Pp = BMA_young.Pp;  
threshold = 0.95; 

M_cov1 = reshape(Ep(1:64), [8,8]);
M_cov2 = reshape(Ep(65:128), [8,8]);
Pp_cov1 = reshape(Pp(1:64), [8,8]);
Pp_cov2 = reshape(Pp(65:128), [8,8]);

%x = M_cov1;
M_cov1(Pp_cov1 <= threshold) = 0;
M_cov2(Pp_cov2 <= threshold) = 0;


figure;
subplot(1,2,1);
imagesc(M_cov1); colorbar;
title('Covariate 1 Connectivity (Filtered by Pp > 0.95)');
set(gca, 'XTick',1:8, 'YTick',1:8);
subplot(1,2,2);
imagesc(M_cov2); colorbar;
title('Covariate 2 Connectivity (Filtered by Pp > 0.95)');
set(gca, 'XTick',1:8, 'YTick',1:8);

%% Average for young

DMN = 1:4;   
DAN = 5:8;   
valid_out_cov1 = M_cov1(DMN, DAN); % DAN →  DMN connections
valid_in_cov1  = M_cov1(DAN, DMN); % DMN →  DAN connections


avg_out_cov1 = sum(valid_out_cov1(:)) / nnz(valid_out_cov1);
avg_in_cov1  = sum(valid_in_cov1(:)) / nnz(valid_in_cov1);
hierarchy_DMN_cov1_young = abs(avg_in_cov1) - abs(avg_out_cov1);
hierarchy_DAN_cov1_young = abs(avg_out_cov1) - abs(avg_in_cov1);


valid_out_cov2 = M_cov2(DMN, DAN); % DMN → DAN connections
valid_in_cov2  = M_cov2(DAN, DMN); % DAN → DMN connections

avg_out_cov2 = sum(valid_out_cov2(:)) / nnz(valid_out_cov2);
avg_in_cov2  = sum(valid_in_cov2(:)) / nnz(valid_in_cov2);
hierarchy_DMN_cov2_young = abs(avg_in_cov2) - abs(avg_out_cov2);
hierarchy_DAN_cov2_young = abs(avg_out_cov2) - abs(avg_in_cov2);


