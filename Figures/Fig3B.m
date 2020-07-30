%% Plot lower uptake of each metal ion
% PCA
load('sLU_res.mat');

%% PCA plot
ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};
ion_clr_list = [228,26,28           % CA
                55,126,184          % CU
                255,127,0           % FE
                77,175,74           % K
                152,78,163          % MG
                166,86,40           % MN
                247,129,191         % NA
                153,153,153]/255;   % ZN
ref_clr = [0,0,0];

figure('Name','PCA');