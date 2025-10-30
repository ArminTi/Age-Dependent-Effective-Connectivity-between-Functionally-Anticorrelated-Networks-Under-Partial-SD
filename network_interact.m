clear all;
clc;
close all;

load('BMA_young.mat');
load('BMA_old.mat');

roiLabels = {'MPFC', 'LP(L)', 'LP(R)', 'PCC', 'FEF(L)', 'FEF(R)', 'IPS(L)', 'IPS(R)', ...
             'ACC', 'AInsula(L)', 'AInsula(R)', 'RPFC(L)', 'RPFC(R)', 'SMG(L)', 'SMG(R)'};

DMN = 1:4;   % (MPFC, LP(L), LP(R), PCC)
DAN = 5:8;   % (FEF(L), FEF(R), IPS(L), IPS(R))
SN = 9:15;   % (ACC, AInsula(L), AInsula(R), RPFC(L), RPFC(R), SMG(L), SMG(R))

threshold = 0.99;

%% Analysis for Old Group
Ep = BMA_old.Ep;
Pp = BMA_old.Pp;

M_cov1 = full(reshape(Ep(1:225), [15, 15]));        % MeanEffect
M_cov2 = full(reshape(Ep(226:450), [15, 15]));      % SleepDeprivationEffect
M_cov3 = full(reshape(Ep(451:675), [15, 15]));      % KSSScore
Pp_cov1 = full(reshape(Pp(1:225), [15, 15]));
Pp_cov2 = full(reshape(Pp(226:450), [15, 15]));
Pp_cov3 = full(reshape(Pp(451:675), [15, 15]));

M_cov1(Pp_cov1 <= threshold) = 0;
M_cov2(Pp_cov2 <= threshold) = 0;
M_cov3(Pp_cov3 <= threshold) = 0;

% Covariance 1 (MeanEffect)
% Reversed: M_cov1(DAN, DMN) for DMN -> DAN, meaning from DAN to DMN in standard convention
avg_out_DMN_DAN_cov1_old = safe_mean(M_cov1(DAN, DMN));
avg_out_DMN_SN_cov1_old = safe_mean(M_cov1(SN, DMN));
avg_out_DMN_cov1_old = safe_mean([avg_out_DMN_DAN_cov1_old; avg_out_DMN_SN_cov1_old]);
avg_in_DAN_DMN_cov1_old = safe_mean(M_cov1(DMN, DAN));
avg_in_SN_DMN_cov1_old = safe_mean(M_cov1(DMN, SN));
avg_in_DMN_cov1_old = safe_mean([avg_in_DAN_DMN_cov1_old; avg_in_SN_DMN_cov1_old]);
hierarchy_DMN_cov1_old = abs(avg_in_DMN_cov1_old) - abs(avg_out_DMN_cov1_old);

avg_out_DAN_DMN_cov1_old = safe_mean(M_cov1(DMN, DAN));
avg_out_DAN_SN_cov1_old = safe_mean(M_cov1(SN, DAN));
avg_out_DAN_cov1_old = safe_mean([avg_out_DAN_DMN_cov1_old; avg_out_DAN_SN_cov1_old]);
avg_in_DMN_DAN_cov1_old = safe_mean(M_cov1(DAN, DMN));
avg_in_SN_DAN_cov1_old = safe_mean(M_cov1(DAN, SN));
avg_in_DAN_cov1_old = safe_mean([avg_in_DMN_DAN_cov1_old; avg_in_SN_DAN_cov1_old]);
hierarchy_DAN_cov1_old = abs(avg_in_DAN_cov1_old) - abs(avg_out_DAN_cov1_old);

avg_out_SN_DMN_cov1_old = safe_mean(M_cov1(DMN, SN));
avg_out_SN_DAN_cov1_old = safe_mean(M_cov1(DAN, SN));
avg_out_SN_cov1_old = safe_mean([avg_out_SN_DMN_cov1_old; avg_out_SN_DAN_cov1_old]);
avg_in_DMN_SN_cov1_old = safe_mean(M_cov1(SN, DMN));
avg_in_DAN_SN_cov1_old = safe_mean(M_cov1(SN, DAN));
avg_in_SN_cov1_old = safe_mean([avg_in_DMN_SN_cov1_old; avg_in_DAN_SN_cov1_old]);
hierarchy_SN_cov1_old = abs(avg_in_SN_cov1_old) - abs(avg_out_SN_cov1_old);

% Display connectivity table for Covariance 1
disp('Old Group - Mean Effect Connectivity:');
T_old_cov1 = table({'DMN -> DAN'; 'DMN -> SN'; 'DAN -> DMN'; 'DAN -> SN'; 'SN -> DMN'; 'SN -> DAN'}, ...
    [avg_out_DMN_DAN_cov1_old; avg_out_DMN_SN_cov1_old; avg_out_DAN_DMN_cov1_old; ...
     avg_out_DAN_SN_cov1_old; avg_out_SN_DMN_cov1_old; avg_out_SN_DAN_cov1_old], ...
    'VariableNames', {'Connection', 'Average_Weight'});
disp(T_old_cov1);

% Covariance 2 (SleepDeprivationEffect)
avg_out_DMN_DAN_cov2_old = safe_mean(M_cov2(DAN, DMN));
avg_out_DMN_SN_cov2_old = safe_mean(M_cov2(SN, DMN));
avg_out_DMN_cov2_old = safe_mean([avg_out_DMN_DAN_cov2_old; avg_out_DMN_SN_cov2_old]);
avg_in_DAN_DMN_cov2_old = safe_mean(M_cov2(DMN, DAN));
avg_in_SN_DMN_cov2_old = safe_mean(M_cov2(DMN, SN));
avg_in_DMN_cov2_old = safe_mean([avg_in_DAN_DMN_cov2_old; avg_in_SN_DMN_cov2_old]);
hierarchy_DMN_cov2_old = abs(avg_in_DMN_cov2_old) - abs(avg_out_DMN_cov2_old);

avg_out_DAN_DMN_cov2_old = safe_mean(M_cov2(DMN, DAN));
avg_out_DAN_SN_cov2_old = safe_mean(M_cov2(SN, DAN));
avg_out_DAN_cov2_old = safe_mean([avg_out_DAN_DMN_cov2_old; avg_out_DAN_SN_cov2_old]);
avg_in_DMN_DAN_cov2_old = safe_mean(M_cov2(DAN, DMN));
avg_in_SN_DAN_cov2_old = safe_mean(M_cov2(DAN, SN));
avg_in_DAN_cov2_old = safe_mean([avg_in_DMN_DAN_cov2_old; avg_in_SN_DAN_cov2_old]);
hierarchy_DAN_cov2_old = abs(avg_in_DAN_cov2_old) - abs(avg_out_DAN_cov2_old);

avg_out_SN_DMN_cov2_old = safe_mean(M_cov2(DMN, SN));
avg_out_SN_DAN_cov2_old = safe_mean(M_cov2(DAN, SN));
avg_out_SN_cov2_old = safe_mean([avg_out_SN_DMN_cov2_old; avg_out_SN_DAN_cov2_old]);
avg_in_DMN_SN_cov2_old = safe_mean(M_cov2(SN, DMN));
avg_in_DAN_SN_cov2_old = safe_mean(M_cov2(SN, DAN));
avg_in_SN_cov2_old = safe_mean([avg_in_DMN_SN_cov2_old; avg_in_DAN_SN_cov2_old]);
hierarchy_SN_cov2_old = abs(avg_in_SN_cov2_old) - abs(avg_out_SN_cov2_old);

% Display connectivity table for Covariance 2
disp('Old Group - Sleep Deprivation Effect Connectivity:');
T_old_cov2 = table({'DMN -> DAN'; 'DMN -> SN'; 'DAN -> DMN'; 'DAN -> SN'; 'SN -> DMN'; 'SN -> DAN'}, ...
    [avg_out_DMN_DAN_cov2_old; avg_out_DMN_SN_cov2_old; avg_out_DAN_DMN_cov2_old; ...
     avg_out_DAN_SN_cov2_old; avg_out_SN_DMN_cov2_old; avg_out_SN_DAN_cov2_old], ...
    'VariableNames', {'Connection', 'Average_Weight'});
disp(T_old_cov2);

% Covariance 3 (KSSScore)
avg_out_DMN_DAN_cov3_old = safe_mean(M_cov3(DAN, DMN));
avg_out_DMN_SN_cov3_old = safe_mean(M_cov3(SN, DMN));
avg_out_DMN_cov3_old = safe_mean([avg_out_DMN_DAN_cov3_old; avg_out_DMN_SN_cov3_old]);
avg_in_DAN_DMN_cov3_old = safe_mean(M_cov3(DMN, DAN));
avg_in_SN_DMN_cov3_old = safe_mean(M_cov3(DMN, SN));
avg_in_DMN_cov3_old = safe_mean([avg_in_DAN_DMN_cov3_old; avg_in_SN_DMN_cov3_old]);
hierarchy_DMN_cov3_old = abs(avg_in_DMN_cov3_old) - abs(avg_out_DMN_cov3_old);

avg_out_DAN_DMN_cov3_old = safe_mean(M_cov3(DMN, DAN));
avg_out_DAN_SN_cov3_old = safe_mean(M_cov3(SN, DAN));
avg_out_DAN_cov3_old = safe_mean([avg_out_DAN_DMN_cov3_old; avg_out_DAN_SN_cov3_old]);
avg_in_DMN_DAN_cov3_old = safe_mean(M_cov3(DAN, DMN));
avg_in_SN_DAN_cov3_old = safe_mean(M_cov3(DAN, SN));
avg_in_DAN_cov3_old = safe_mean([avg_in_DMN_DAN_cov3_old; avg_in_SN_DAN_cov3_old]);
hierarchy_DAN_cov3_old = abs(avg_in_DAN_cov3_old) - abs(avg_out_DAN_cov3_old);

avg_out_SN_DMN_cov3_old = safe_mean(M_cov3(DMN, SN));
avg_out_SN_DAN_cov3_old = safe_mean(M_cov3(DAN, SN));
avg_out_SN_cov3_old = safe_mean([avg_out_SN_DMN_cov3_old; avg_out_SN_DAN_cov3_old]);
avg_in_DMN_SN_cov3_old = safe_mean(M_cov3(SN, DMN));
avg_in_DAN_SN_cov3_old = safe_mean(M_cov3(SN, DAN));
avg_in_SN_cov3_old = safe_mean([avg_in_DMN_SN_cov3_old; avg_in_DAN_SN_cov3_old]);
hierarchy_SN_cov3_old = abs(avg_in_SN_cov3_old) - abs(avg_out_SN_cov3_old);

% Display connectivity table for Covariance 3
disp('Old Group - KSS Score Connectivity:');
T_old_cov3 = table({'DMN -> DAN'; 'DMN -> SN'; 'DAN -> DMN'; 'DAN -> SN'; 'SN -> DMN'; 'SN -> DAN'}, ...
    [avg_out_DMN_DAN_cov3_old; avg_out_DMN_SN_cov3_old; avg_out_DAN_DMN_cov3_old; ...
     avg_out_DAN_SN_cov3_old; avg_out_SN_DMN_cov3_old; avg_out_SN_DAN_cov3_old], ...
    'VariableNames', {'Connection', 'Average_Weight'});
disp(T_old_cov3);

% Directed graphs for Old Group
weights_old_cov1 = [avg_out_DMN_DAN_cov1_old, avg_out_DMN_SN_cov1_old, avg_out_DAN_DMN_cov1_old, ...
                    avg_out_DAN_SN_cov1_old, avg_out_SN_DMN_cov1_old, avg_out_SN_DAN_cov1_old];
weights_old_cov1(isnan(weights_old_cov1)) = 0;

weights_old_cov2 = [avg_out_DMN_DAN_cov2_old, avg_out_DMN_SN_cov2_old, avg_out_DAN_DMN_cov2_old, ...
                    avg_out_DAN_SN_cov2_old, avg_out_SN_DMN_cov2_old, avg_out_SN_DAN_cov2_old];
weights_old_cov2(isnan(weights_old_cov2)) = 0;

weights_old_cov3 = [avg_out_DMN_DAN_cov3_old, avg_out_DMN_SN_cov3_old, avg_out_DAN_DMN_cov3_old, ...
                    avg_out_DAN_SN_cov3_old, avg_out_SN_DMN_cov3_old, avg_out_SN_DAN_cov3_old];
weights_old_cov3(isnan(weights_old_cov3)) = 0;

node_colors = [0.1, 0.7, 0.3; 0.9, 0.4, 0.2; 0.3, 0.5, 0.8]; % DMN, DAN, SN

figure('Name', 'Old Group Inter-Network Connectivity', 'Color', [0.9 0.9 0.9]);
subplot(1, 3, 1);
% Reversed: Edge [1,2] is now DMN -> DAN, which uses avg_out_DMN_DAN (M_cov1(DAN, DMN))
G_old_cov1 = digraph([1, 1, 2, 2, 3, 3], [2, 3, 1, 3, 1, 2], weights_old_cov1, {'DMN', 'DAN', 'SN'});
h = plot(G_old_cov1, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_old_cov1.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, ...
         'EdgeCData', weights_old_cov1, 'NodeLabel', {});
colormap([0 0 1; 1 0 0]);
h.LineWidth = abs(weights_old_cov1) * 5 + 1;
title('Mean Effect Connectivity (Old Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(1, 3, 2);
G_old_cov2 = digraph([1, 1, 2, 2, 3, 3], [2, 3, 1, 3, 1, 2], weights_old_cov2, {'DMN', 'DAN', 'SN'});
h = plot(G_old_cov2, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_old_cov2.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, ...
         'EdgeCData', weights_old_cov2, 'NodeLabel', {});
colormap([0 0 1; 1 0 0]);
h.LineWidth = abs(weights_old_cov2) * 5 + 1;
title('Sleep Deprivation Effect Connectivity (Old Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(1, 3, 3);
G_old_cov3 = digraph([1, 1, 2, 2, 3, 3], [2, 3, 1, 3, 1, 2], weights_old_cov3, {'DMN', 'DAN', 'SN'});
h = plot(G_old_cov3, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_old_cov3.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, ...
         'EdgeCData', weights_old_cov3, 'NodeLabel', {});
colormap([0 0 1; 1 0 0]);
h.LineWidth = abs(weights_old_cov3) * 5 + 1;
title('KSS Score Connectivity (Old Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

%% Analysis for Young Group
Ep = BMA_young.Ep;
Pp = BMA_young.Pp;

M_cov1 = full(reshape(Ep(1:225), [15, 15]));        % MeanEffect
M_cov2 = full(reshape(Ep(226:450), [15, 15]));      % SleepDeprivationEffect
M_cov3 = full(reshape(Ep(451:675), [15, 15]));      % KSSScore
Pp_cov1 = full(reshape(Pp(1:225), [15, 15]));
Pp_cov2 = full(reshape(Pp(226:450), [15, 15]));
Pp_cov3 = full(reshape(Pp(451:675), [15, 15]));

M_cov1(Pp_cov1 <= threshold) = 0;
M_cov2(Pp_cov2 <= threshold) = 0;
M_cov3(Pp_cov3 <= threshold) = 0;

% Covariance 1 (MeanEffect)
avg_out_DMN_DAN_cov1_young = safe_mean(M_cov1(DAN, DMN));
avg_out_DMN_SN_cov1_young = safe_mean(M_cov1(SN, DMN));
avg_out_DMN_cov1_young = safe_mean([avg_out_DMN_DAN_cov1_young; avg_out_DMN_SN_cov1_young]);
avg_in_DAN_DMN_cov1_young = safe_mean(M_cov1(DMN, DAN));
avg_in_SN_DMN_cov1_young = safe_mean(M_cov1(DMN, SN));
avg_in_DMN_cov1_young = safe_mean([avg_in_DAN_DMN_cov1_young; avg_in_SN_DMN_cov1_young]);
hierarchy_DMN_cov1_young = abs(avg_in_DMN_cov1_young) - abs(avg_out_DMN_cov1_young);

avg_out_DAN_DMN_cov1_young = safe_mean(M_cov1(DMN, DAN));
avg_out_DAN_SN_cov1_young = safe_mean(M_cov1(SN, DAN));
avg_out_DAN_cov1_young = safe_mean([avg_out_DAN_DMN_cov1_young; avg_out_DAN_SN_cov1_young]);
avg_in_DMN_DAN_cov1_young = safe_mean(M_cov1(DAN, DMN));
avg_in_SN_DAN_cov1_young = safe_mean(M_cov1(DAN, SN));
avg_in_DAN_cov1_young = safe_mean([avg_in_DMN_DAN_cov1_young; avg_in_SN_DAN_cov1_young]);
hierarchy_DAN_cov1_young = abs(avg_in_DAN_cov1_young) - abs(avg_out_DAN_cov1_young);

avg_out_SN_DMN_cov1_young = safe_mean(M_cov1(DMN, SN));
avg_out_SN_DAN_cov1_young = safe_mean(M_cov1(DAN, SN));
avg_out_SN_cov1_young = safe_mean([avg_out_SN_DMN_cov1_young; avg_out_SN_DAN_cov1_young]);
avg_in_DMN_SN_cov1_young = safe_mean(M_cov1(SN, DMN));
avg_in_DAN_SN_cov1_young = safe_mean(M_cov1(SN, DAN));
avg_in_SN_cov1_young = safe_mean([avg_in_DMN_SN_cov1_young; avg_in_DAN_SN_cov1_young]);
hierarchy_SN_cov1_young = abs(avg_in_SN_cov1_young) - abs(avg_out_SN_cov1_young);

% Display connectivity table for Covariance 1
disp('Young Group - Mean Effect Connectivity:');
T_young_cov1 = table({'DMN -> DAN'; 'DMN -> SN'; 'DAN -> DMN'; 'DAN -> SN'; 'SN -> DMN'; 'SN -> DAN'}, ...
    [avg_out_DMN_DAN_cov1_young; avg_out_DMN_SN_cov1_young; avg_out_DAN_DMN_cov1_young; ...
     avg_out_DAN_SN_cov1_young; avg_out_SN_DMN_cov1_young; avg_out_SN_DAN_cov1_young], ...
    'VariableNames', {'Connection', 'Average_Weight'});
disp(T_young_cov1);

% Covariance 2 (SleepDeprivationEffect)
avg_out_DMN_DAN_cov2_young = safe_mean(M_cov2(DAN, DMN));
avg_out_DMN_SN_cov2_young = safe_mean(M_cov2(SN, DMN));
avg_out_DMN_cov2_young = safe_mean([avg_out_DMN_DAN_cov2_young; avg_out_DMN_SN_cov2_young]);
avg_in_DAN_DMN_cov2_young = safe_mean(M_cov2(DMN, DAN));
avg_in_SN_DMN_cov2_young = safe_mean(M_cov2(DMN, SN));
avg_in_DMN_cov2_young = safe_mean([avg_in_DAN_DMN_cov2_young; avg_in_SN_DMN_cov2_young]);
hierarchy_DMN_cov2_young = abs(avg_in_DMN_cov2_young) - abs(avg_out_DMN_cov2_young);

avg_out_DAN_DMN_cov2_young = safe_mean(M_cov2(DMN, DAN));
avg_out_DAN_SN_cov2_young = safe_mean(M_cov2(SN, DAN));
avg_out_DAN_cov2_young = safe_mean([avg_out_DAN_DMN_cov2_young; avg_out_DAN_SN_cov2_young]);
avg_in_DMN_DAN_cov2_young = safe_mean(M_cov2(DAN, DMN));
avg_in_SN_DAN_cov2_young = safe_mean(M_cov2(DAN, SN));
avg_in_DAN_cov2_young = safe_mean([avg_in_DMN_DAN_cov2_young; avg_in_SN_DAN_cov2_young]);
hierarchy_DAN_cov2_young = abs(avg_in_DAN_cov2_young) - abs(avg_out_DAN_cov2_young);

avg_out_SN_DMN_cov2_young = safe_mean(M_cov2(DMN, SN));
avg_out_SN_DAN_cov2_young = safe_mean(M_cov2(DAN, SN));
avg_out_SN_cov2_young = safe_mean([avg_out_SN_DMN_cov2_young; avg_out_SN_DAN_cov2_young]);
avg_in_DMN_SN_cov2_young = safe_mean(M_cov2(SN, DMN));
avg_in_DAN_SN_cov2_young = safe_mean(M_cov2(SN, DAN));
avg_in_SN_cov2_young = safe_mean([avg_in_DMN_SN_cov2_young; avg_in_DAN_SN_cov2_young]);
hierarchy_SN_cov2_young = abs(avg_in_SN_cov2_young) - abs(avg_out_SN_cov2_young);

% Display connectivity table for Covariance 2
disp('Young Group - Sleep Deprivation Effect Connectivity:');
T_young_cov2 = table({'DMN -> DAN'; 'DMN -> SN'; 'DAN -> DMN'; 'DAN -> SN'; 'SN -> DMN'; 'SN -> DAN'}, ...
    [avg_out_DMN_DAN_cov2_young; avg_out_DMN_SN_cov2_young; avg_out_DAN_DMN_cov2_young; ...
     avg_out_DAN_SN_cov2_young; avg_out_SN_DMN_cov2_young; avg_out_SN_DAN_cov2_young], ...
    'VariableNames', {'Connection', 'Average_Weight'});
disp(T_young_cov2);

% Covariance 3 (KSSScore)
avg_out_DMN_DAN_cov3_young = safe_mean(M_cov3(DAN, DMN));
avg_out_DMN_SN_cov3_young = safe_mean(M_cov3(SN, DMN));
avg_out_DMN_cov3_young = safe_mean([avg_out_DMN_DAN_cov3_young; avg_out_DMN_SN_cov3_young]);
avg_in_DAN_DMN_cov3_young = safe_mean(M_cov3(DMN, DAN));
avg_in_SN_DMN_cov3_young = safe_mean(M_cov3(DMN, SN));
avg_in_DMN_cov3_young = safe_mean([avg_in_DAN_DMN_cov3_young; avg_in_SN_DMN_cov3_young]);
hierarchy_DMN_cov3_young = abs(avg_in_DMN_cov3_young) - abs(avg_out_DMN_cov3_young);

avg_out_DAN_DMN_cov3_young = safe_mean(M_cov3(DMN, DAN));
avg_out_DAN_SN_cov3_young = safe_mean(M_cov3(SN, DAN));
avg_out_DAN_cov3_young = safe_mean([avg_out_DAN_DMN_cov3_young; avg_out_DAN_SN_cov3_young]);
avg_in_DMN_DAN_cov3_young = safe_mean(M_cov3(DAN, DMN));
avg_in_SN_DAN_cov3_young = safe_mean(M_cov3(DAN, SN));
avg_in_DAN_cov3_young = safe_mean([avg_in_DMN_DAN_cov3_young; avg_in_SN_DAN_cov3_young]);
hierarchy_DAN_cov3_young = abs(avg_in_DAN_cov3_young) - abs(avg_out_DAN_cov3_young);

avg_out_SN_DMN_cov3_young = safe_mean(M_cov3(DMN, SN));
avg_out_SN_DAN_cov3_young = safe_mean(M_cov3(DAN, SN));
avg_out_SN_cov3_young = safe_mean([avg_out_SN_DMN_cov3_young; avg_out_SN_DAN_cov3_young]);
avg_in_DMN_SN_cov3_young = safe_mean(M_cov3(SN, DMN));
avg_in_DAN_SN_cov3_young = safe_mean(M_cov3(SN, DAN));
avg_in_SN_cov3_young = safe_mean([avg_in_DMN_SN_cov3_young; avg_in_DAN_SN_cov3_young]);
hierarchy_SN_cov3_young = abs(avg_in_SN_cov3_young) - abs(avg_out_SN_cov3_young);

% Display connectivity table for Covariance 3
disp('Young Group - KSS Score Connectivity:');
T_young_cov3 = table({'DMN -> DAN'; 'DMN -> SN'; 'DAN -> DMN'; 'DAN -> SN'; 'SN -> DMN'; 'SN -> DAN'}, ...
    [avg_out_DMN_DAN_cov3_young; avg_out_DMN_SN_cov3_young; avg_out_DAN_DMN_cov3_young; ...
     avg_out_DAN_SN_cov3_young; avg_out_SN_DMN_cov3_young; avg_out_SN_DAN_cov3_young], ...
    'VariableNames', {'Connection', 'Average_Weight'});
disp(T_young_cov3);

% Directed graphs for Young Group
weights_young_cov1 = [avg_out_DMN_DAN_cov1_young, avg_out_DMN_SN_cov1_young, avg_out_DAN_DMN_cov1_young, ...
                      avg_out_DAN_SN_cov1_young, avg_out_SN_DMN_cov1_young, avg_out_SN_DAN_cov1_young];
weights_young_cov1(isnan(weights_young_cov1)) = 0;

weights_young_cov2 = [avg_out_DMN_DAN_cov2_young, avg_out_DMN_SN_cov2_young, avg_out_DAN_DMN_cov2_young, ...
                      avg_out_DAN_SN_cov2_young, avg_out_SN_DMN_cov2_young, avg_out_SN_DAN_cov2_young];
weights_young_cov2(isnan(weights_young_cov2)) = 0;

weights_young_cov3 = [avg_out_DMN_DAN_cov3_young, avg_out_DMN_SN_cov3_young, avg_out_DAN_DMN_cov3_young, ...
                      avg_out_DAN_SN_cov3_young, avg_out_SN_DMN_cov3_young, avg_out_SN_DAN_cov3_young];
weights_young_cov3(isnan(weights_young_cov3)) = 0;

figure('Name', 'Young Group Inter-Network Connectivity', 'Color', [0.9 0.9 0.9]);
subplot(1, 3, 1);
G_young_cov1 = digraph([1, 1, 2, 2, 3, 3], [2, 3, 1, 3, 1, 2], weights_young_cov1, {'DMN', 'DAN', 'SN'});
h = plot(G_young_cov1, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_young_cov1.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, ...
         'EdgeCData', weights_young_cov1, 'NodeLabel', {});
colormap([0 0 1; 1 0 0]);
h.LineWidth = abs(weights_young_cov1) * 5 + 1;
title('Mean Effect Connectivity (Young Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(1, 3, 2);
G_young_cov2 = digraph([1, 1, 2, 2, 3, 3], [2, 3, 1, 3, 1, 2], weights_young_cov2, {'DMN', 'DAN', 'SN'});
h = plot(G_young_cov2, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_young_cov2.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, ...
         'EdgeCData', weights_young_cov2, 'NodeLabel', {});
colormap([0 0 1; 1 0 0]);
h.LineWidth = abs(weights_young_cov2) * 5 + 1;
title('Sleep Deprivation Effect Connectivity (Young Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(1, 3, 3);
G_young_cov3 = digraph([1, 1, 2, 2, 3, 3], [2, 3, 1, 3, 1, 2], weights_young_cov3, {'DMN', 'DAN', 'SN'});
h = plot(G_young_cov3, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_young_cov3.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, ...
         'EdgeCData', weights_young_cov3, 'NodeLabel', {});
colormap([0 0 1; 1 0 0]);
h.LineWidth = abs(weights_young_cov3) * 5 + 1;
title('KSS Score Connectivity (Young Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

% Hierarchy Directed Graphs
figure('Name', 'Hierarchy as Directed Graph', 'Color', [0.9 0.9 0.9]);
subplot(2, 3, 1);
[hierarchy_sources, hierarchy_targets, hierarchy_weights] = create_hierarchy_edges([hierarchy_DMN_cov1_old, hierarchy_DAN_cov1_old, hierarchy_SN_cov1_old]);
G_hierarchy_old_cov1 = digraph(hierarchy_sources, hierarchy_targets, hierarchy_weights, {'DMN', 'DAN', 'SN'});
h = plot(G_hierarchy_old_cov1, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_hierarchy_old_cov1.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, 'NodeLabel', {});
h.LineWidth = G_hierarchy_old_cov1.Edges.Weight * 5 + 1;
title('Hierarchy: Mean Effect (Old Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(2, 3, 2);
[hierarchy_sources, hierarchy_targets, hierarchy_weights] = create_hierarchy_edges([hierarchy_DMN_cov2_old, hierarchy_DAN_cov2_old, hierarchy_SN_cov2_old]);
G_hierarchy_old_cov2 = digraph(hierarchy_sources, hierarchy_targets, hierarchy_weights, {'DMN', 'DAN', 'SN'});
h = plot(G_hierarchy_old_cov2, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_hierarchy_old_cov2.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, 'NodeLabel', {});
h.LineWidth = G_hierarchy_old_cov2.Edges.Weight * 5 + 1;
title('Hierarchy: Sleep Deprivation (Old Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(2, 3, 3);
[hierarchy_sources, hierarchy_targets, hierarchy_weights] = create_hierarchy_edges([hierarchy_DMN_cov3_old, hierarchy_DAN_cov3_old, hierarchy_SN_cov3_old]);
G_hierarchy_old_cov3 = digraph(hierarchy_sources, hierarchy_targets, hierarchy_weights, {'DMN', 'DAN', 'SN'});
h = plot(G_hierarchy_old_cov3, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_hierarchy_old_cov3.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, 'NodeLabel', {});
h.LineWidth = G_hierarchy_old_cov3.Edges.Weight * 5 + 1;
title('Hierarchy: KSS Score (Old Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(2, 3, 4);
[hierarchy_sources, hierarchy_targets, hierarchy_weights] = create_hierarchy_edges([hierarchy_DMN_cov1_young, hierarchy_DAN_cov1_young, hierarchy_SN_cov1_young]);
G_hierarchy_young_cov1 = digraph(hierarchy_sources, hierarchy_targets, hierarchy_weights, {'DMN', 'DAN', 'SN'});
h = plot(G_hierarchy_young_cov1, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_hierarchy_young_cov1.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, 'NodeLabel', {});
h.LineWidth = G_hierarchy_young_cov1.Edges.Weight * 5 + 1;
title('Hierarchy: Mean Effect (Young Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(2, 3, 5);
[hierarchy_sources, hierarchy_targets, hierarchy_weights] = create_hierarchy_edges([hierarchy_DMN_cov2_young, hierarchy_DAN_cov2_young, hierarchy_SN_cov2_young]);
G_hierarchy_young_cov2 = digraph(hierarchy_sources, hierarchy_targets, hierarchy_weights, {'DMN', 'DAN', 'SN'});
h = plot(G_hierarchy_young_cov2, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_hierarchy_young_cov2.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, 'NodeLabel', {});
h.LineWidth = G_hierarchy_young_cov2.Edges.Weight * 5 + 1;
title('Hierarchy: Sleep Deprivation (Young Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

subplot(2, 3, 6);
[hierarchy_sources, hierarchy_targets, hierarchy_weights] = create_hierarchy_edges([hierarchy_DMN_cov3_young, hierarchy_DAN_cov3_young, hierarchy_SN_cov3_young]);
G_hierarchy_young_cov3 = digraph(hierarchy_sources, hierarchy_targets, hierarchy_weights, {'DMN', 'DAN', 'SN'});
h = plot(G_hierarchy_young_cov3, 'NodeColor', node_colors, 'MarkerSize', 20, 'Marker', 'o', ...
         'EdgeLabel', G_hierarchy_young_cov3.Edges.Weight, 'NodeFontSize', 14, 'EdgeFontSize', 12, 'NodeLabel', {});
h.LineWidth = G_hierarchy_young_cov3.Edges.Weight * 5 + 1;
title('Hierarchy: KSS Score (Young Group)', 'FontSize', 14);
set(gca, 'Color', 'none');

% Bar plots for Hierarchy Comparison
figure('Name', 'Hierarchy Comparison', 'Color', [0.9 0.9 0.9]);
subplot(1, 3, 1);
bar([hierarchy_DMN_cov1_old, hierarchy_DAN_cov1_old, hierarchy_SN_cov1_old; ...
     hierarchy_DMN_cov1_young, hierarchy_DAN_cov1_young, hierarchy_SN_cov1_young]');
set(gca, 'XTickLabel', {'DMN', 'DAN', 'SN'}, 'FontSize', 12);
legend('Old Group', 'Young Group', 'FontSize', 10);
title('Hierarchy (Mean Effect)', 'FontSize', 14);
ylabel('Hierarchy Index', 'FontSize', 12);

subplot(1, 3, 2);
bar([hierarchy_DMN_cov2_old, hierarchy_DAN_cov2_old, hierarchy_SN_cov2_old; ...
     hierarchy_DMN_cov2_young, hierarchy_DAN_cov2_young, hierarchy_SN_cov2_young]');
set(gca, 'XTickLabel', {'DMN', 'DAN', 'SN'}, 'FontSize', 12);
legend('Old Group', 'Young Group', 'FontSize', 10);
title('Hierarchy (Sleep Deprivation Effect)', 'FontSize', 14);
ylabel('Hierarchy Index', 'FontSize', 12);

subplot(1, 3, 3);
bar([hierarchy_DMN_cov3_old, hierarchy_DAN_cov3_old, hierarchy_SN_cov3_old; ...
     hierarchy_DMN_cov3_young, hierarchy_DAN_cov3_young, hierarchy_SN_cov3_young]');
set(gca, 'XTickLabel', {'DMN', 'DAN', 'SN'}, 'FontSize', 12);
legend('Old Group', 'Young Group', 'FontSize', 10);
title('Hierarchy (KSS Score)', 'FontSize', 14);
ylabel('Hierarchy Index', 'FontSize', 12);

% Display hierarchy results in console
disp('Hierarchy for Old Group:');
disp(['DMN (MeanEffect): ', num2str(hierarchy_DMN_cov1_old)]);
disp(['DAN (MeanEffect): ', num2str(hierarchy_DAN_cov1_old)]);
disp(['SN (MeanEffect): ', num2str(hierarchy_SN_cov1_old)]);
disp(['DMN (SleepDeprivationEffect): ', num2str(hierarchy_DMN_cov2_old)]);
disp(['DAN (SleepDeprivationEffect): ', num2str(hierarchy_DAN_cov2_old)]);
disp(['SN (SleepDeprivationEffect): ', num2str(hierarchy_SN_cov2_old)]);
disp(['DMN (KSSScore): ', num2str(hierarchy_DMN_cov3_old)]);
disp(['DAN (KSSScore): ', num2str(hierarchy_DAN_cov3_old)]);
disp(['SN (KSSScore): ', num2str(hierarchy_SN_cov3_old)]);

disp('Hierarchy for Young Group:');
disp(['DMN (MeanEffect): ', num2str(hierarchy_DMN_cov1_young)]);
disp(['DAN (MeanEffect): ', num2str(hierarchy_DAN_cov1_young)]);
disp(['SN (MeanEffect): ', num2str(hierarchy_SN_cov1_young)]);
disp(['DMN (SleepDeprivationEffect): ', num2str(hierarchy_DMN_cov2_young)]);
disp(['DAN (SleepDeprivationEffect): ', num2str(hierarchy_DAN_cov2_young)]);
disp(['SN (SleepDeprivationEffect): ', num2str(hierarchy_SN_cov2_young)]);
disp(['DMN (KSSScore): ', num2str(hierarchy_DMN_cov3_young)]);
disp(['DAN (KSSScore): ', num2str(hierarchy_DAN_cov3_young)]);
disp(['SN (KSSScore): ', num2str(hierarchy_SN_cov3_young)]);

%% Helper function
function avg = safe_mean(data)
    data = full(double(data));
    if nnz(data) == 0
        avg = 0;
    else
        avg = sum(data(:)) / nnz(data);
    end
    avg = double(avg);
end

function [sources, targets, edge_weights] = create_hierarchy_edges(hierarchy_values)
    sources = [];
    targets = [];
    edge_weights = [];
    for i = 1:3
        for j = 1:3
            if i ~= j && hierarchy_values(i) < hierarchy_values(j)
                sources = [sources; i];
                targets = [targets; j];
                edge_weights = [edge_weights; abs(hierarchy_values(j) - hierarchy_values(i))];
            end
        end
    end
end


% -------------------------------
%  Author: Vida Feizi
%  Affiliation: M.Sc. Student, Computer Engineering, University of Tabriz
%  Email: vFeizii75@gmail.com
% -------------------------------

