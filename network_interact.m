clear all;
clc;
close all; 

load('BMA_young.mat'); 
load('BMA_old.mat');
 
DMN = 1:4;
DAN = 5:8;
SN = 9:15;
network_names = {'DMN', 'DAN', 'SN'};
 
threshold = 0.99;
  
hierarchy_results_old = zeros(3, 3); % Rows: DMN, DAN, SN | Cols: Cov1, Cov2, Cov3
hierarchy_results_young = zeros(3, 3);

%% ====================  
 
Ep_old = BMA_old.Ep; %     
Pp_old = BMA_old.Pp;  
 
M_cov1_old = reshape(Ep_old(1:225), 15, 15); % Mean Effect
M_cov2_old = reshape(Ep_old(226:450), 15, 15); % Sleep Deprivation Effect
M_cov3_old = reshape(Ep_old(451:675), 15, 15); % FD Effect
 
Pp_cov1_old = reshape(Pp_old(1:225), 15, 15);
Pp_cov2_old = reshape(Pp_old(226:450), 15, 15);
Pp_cov3_old = reshape(Pp_old(451:675), 15, 15);
 
M_cov1_old(Pp_cov1_old <= threshold) = 0;
M_cov2_old(Pp_cov2_old <= threshold) = 0;
M_cov3_old(Pp_cov3_old <= threshold) = 0;

all_M_cov_old = {M_cov1_old, M_cov2_old, M_cov3_old};
cov_names = {'Mean Effect', 'Sleep Deprivation Effect', 'FD Effect'};

for c = 1:3 % Loop through each covariance matrix
    
    M_cov = all_M_cov_old{c};
    disp(['---------- ' cov_names{c} ' (Old) ----------']);
     
    out_DMN_to_DAN = M_cov(DAN, DMN);
    out_DMN_to_SN = M_cov(SN, DMN);
    avg_out_DMN_old = safe_mean([out_DMN_to_DAN(:); out_DMN_to_SN(:)]);
    
    
    in_DMN_from_DAN = M_cov(DMN, DAN);
    in_DMN_from_SN = M_cov(DMN, SN);
    avg_in_DMN_old = safe_mean([in_DMN_from_DAN(:); in_DMN_from_SN(:)]);
    
    
    out_DAN_to_DMN = M_cov(DMN, DAN);
    out_DAN_to_SN = M_cov(SN, DAN);
    avg_out_DAN_old = safe_mean([out_DAN_to_DMN(:); out_DAN_to_SN(:)]);

    % اتصالات ورودی (از دیگران به DAN)
    in_DAN_from_DMN = M_cov(DAN, DMN);
    in_DAN_from_SN = M_cov(DAN, SN);
    avg_in_DAN_old = safe_mean([in_DAN_from_DMN(:); in_DAN_from_SN(:)]);
    
 
    out_SN_to_DMN = M_cov(DMN, SN);
    out_SN_to_DAN = M_cov(DAN, SN);
    avg_out_SN_old = safe_mean([out_SN_to_DMN(:); out_SN_to_DAN(:)]);
 
    in_SN_from_DMN = M_cov(SN, DMN);
    in_SN_from_DAN = M_cov(SN, DAN);
    avg_in_SN_old = safe_mean([in_SN_from_DMN(:); in_SN_from_DAN(:)]);
 
    % Hierarchy = |mean(outgoing)| - |mean(incoming)|
    hierarchy_DMN = abs(avg_out_DMN_old) - abs(avg_in_DMN_old);
    hierarchy_DAN = abs(avg_out_DAN_old) - abs(avg_in_DAN_old);
    hierarchy_SN = abs(avg_out_SN_old) - abs(avg_in_SN_old);
    
    hierarchy_results_old(:, c) = [hierarchy_DMN; hierarchy_DAN; hierarchy_SN];
     
    From = {'DMN'; 'DMN'; 'DAN'; 'DAN'; 'SN'; 'SN'};
    To = {'DAN'; 'SN'; 'DMN'; 'SN'; 'DMN'; 'DAN'};
    Mean_Connectivity = [
        safe_mean(out_DMN_to_DAN(:)); safe_mean(out_DMN_to_SN(:));
        safe_mean(out_DAN_to_DMN(:)); safe_mean(out_DAN_to_SN(:));
        safe_mean(out_SN_to_DMN(:)); safe_mean(out_SN_to_DAN(:))
    ];
    tbl = table(From, To, Mean_Connectivity);
    disp(tbl);
     
    disp(['DMN: ', num2str(hierarchy_DMN)]);
    disp(['DAN: ', num2str(hierarchy_DAN)]);
    disp(['SN:  ', num2str(hierarchy_SN)]);
    disp('------------------------------------');
end


%% ======================= 
 
Ep_young = BMA_young.Ep;
Pp_young = BMA_young.Pp;
 
M_cov1_young = reshape(Ep_young(1:225), 15, 15);
M_cov2_young = reshape(Ep_young(226:450), 15, 15);
M_cov3_young = reshape(Ep_young(451:675), 15, 15);

Pp_cov1_young = reshape(Pp_young(1:225), 15, 15);
Pp_cov2_young = reshape(Pp_young(226:450), 15, 15);
Pp_cov3_young = reshape(Pp_young(451:675), 15, 15); 

M_cov1_young(Pp_cov1_young <= threshold) = 0;
M_cov2_young(Pp_cov2_young <= threshold) = 0;
M_cov3_young(Pp_cov3_young <= threshold) = 0;

all_M_cov_young = {M_cov1_young, M_cov2_young, M_cov3_young};

for c = 1:3  
    
    M_cov = all_M_cov_young{c};
    disp(['---------- ' cov_names{c} ' (Young) ----------']);
     
    out_DMN_to_DAN = M_cov(DAN, DMN);
    out_DMN_to_SN = M_cov(SN, DMN);
    avg_out_DMN_young = safe_mean([out_DMN_to_DAN(:); out_DMN_to_SN(:)]);
    
    in_DMN_from_DAN = M_cov(DMN, DAN);
    in_DMN_from_SN = M_cov(DMN, SN);
    avg_in_DMN_young = safe_mean([in_DMN_from_DAN(:); in_DMN_from_SN(:)]);
     
    out_DAN_to_DMN = M_cov(DMN, DAN);
    out_DAN_to_SN = M_cov(SN, DAN);
    avg_out_DAN_young = safe_mean([out_DAN_to_DMN(:); out_DAN_to_SN(:)]);

    in_DAN_from_DMN = M_cov(DAN, DMN);
    in_DAN_from_SN = M_cov(DAN, SN);
    avg_in_DAN_young = safe_mean([in_DAN_from_DMN(:); in_DAN_from_SN(:)]);
     
    out_SN_to_DMN = M_cov(DMN, SN);
    out_SN_to_DAN = M_cov(DAN, SN);
    avg_out_SN_young = safe_mean([out_SN_to_DMN(:); out_SN_to_DAN(:)]);

    in_SN_from_DMN = M_cov(SN, DMN);
    in_SN_from_DAN = M_cov(SN, DAN);
    avg_in_SN_young = safe_mean([in_SN_from_DMN(:); in_SN_from_DAN(:)]);
 
    hierarchy_DMN = abs(avg_out_DMN_young) - abs(avg_in_DMN_young);
    hierarchy_DAN = abs(avg_out_DAN_young) - abs(avg_in_DAN_young);
    hierarchy_SN = abs(avg_out_SN_young) - abs(avg_in_SN_young);
    
    hierarchy_results_young(:, c) = [hierarchy_DMN; hierarchy_DAN; hierarchy_SN];
     
    From = {'DMN'; 'DMN'; 'DAN'; 'DAN'; 'SN'; 'SN'};
    To = {'DAN'; 'SN'; 'DMN'; 'SN'; 'DMN'; 'DAN'};
    Mean_Connectivity = [
        safe_mean(out_DMN_to_DAN(:)); safe_mean(out_DMN_to_SN(:));
        safe_mean(out_DAN_to_DMN(:)); safe_mean(out_DAN_to_SN(:));
        safe_mean(out_SN_to_DMN(:)); safe_mean(out_SN_to_DAN(:))
    ];
    tbl = table(From, To, Mean_Connectivity);
    disp(tbl);
     
    disp(['DMN: ', num2str(hierarchy_DMN)]);
    disp(['DAN: ', num2str(hierarchy_DAN)]);
    disp(['SN:  ', num2str(hierarchy_SN)]);
    disp('------------------------------------');
end


%% =======   
figure('Name', '  (Mean Effect)', 'Position', [100, 100, 800, 600]);
bar_data = [hierarchy_results_old(:,1), hierarchy_results_young(:,1)];
b = bar(bar_data);
set(gca, 'xticklabel', network_names);
legend('Old Group', 'Young Group');
title('Hierarchy Strength Comparison (Mean Effect)');
ylabel('Hierarchy Value (|Out| - |In|)');
grid on;
 
figure('Name', 'Graph heir', 'Position', [950, 100, 1000, 500]);

% گراف گروه Old
subplot(1, 2, 1);
G_old = create_hierarchy_edges(hierarchy_results_old(:,1), network_names);
p_old = plot(G_old, 'Layout', 'force', 'EdgeLabel', G_old.Edges.Weight);
title('Hierarchy Graph - Old Group (Mean Effect)');
p_old.NodeColor = 'r';
 
subplot(1, 2, 2);
G_young = create_hierarchy_edges(hierarchy_results_young(:,1), network_names);
p_young = plot(G_young, 'Layout', 'force', 'EdgeLabel', G_young.Edges.Weight);
title('Hierarchy Graph - Young Group (Mean Effect)');
p_young.NodeColor = 'b';


%% ============================================== 

function m = safe_mean(data) 
    non_zero_data = data(data ~= 0);
    if isempty(non_zero_data)
        m = 0;
    else
        m = mean(non_zero_data);
    end
end

function G = create_hierarchy_edges(hierarchy_values, names) 
    num_nodes = length(hierarchy_values);
    s = []; % Source nodes
    t = []; % Target nodes
    weights = []; % Edge weights
    for i = 1:num_nodes
        for j = 1:num_nodes
            if i ~= j
                if hierarchy_values(i) < hierarchy_values(j)
                    s = [s, i];
                    t = [t, j];
                    weights = [weights, abs(hierarchy_values(j) - hierarchy_values(i))];
                end
            end
        end
    end
    G = digraph(s, t, weights, names);
end
 

% -------------------------------
%  Author: Vida Feizi
%  Affiliation: M.Sc. Student, Computer Engineering, University of Tabriz
%  Email: vFeizii75@gmail.com
% -------------------------------

