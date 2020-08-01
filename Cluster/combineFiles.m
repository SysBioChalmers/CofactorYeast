%% combineFiles
cd tmp_results/;

%% simulationBiolog
% sB_res.fluxes = zeros(0,0);
% sB_res.labels = cell(0,0);
% k = 1:1:116;
% for i = 1:length(k)
%     display([num2str(i),'/',num2str(length(k))]);
%     file_flux = ['sB_fluxes_',num2str(k(i)),'.mat'];
%     file_label = ['sB_labels_',num2str(k(i)),'.mat'];
%     load(file_flux);
%     load(file_label);
% %     fluxes = fluxes(:,k(i):end);
% %     labels = labels(1,k(i):end);
%     sB_res.fluxes = [sB_res.fluxes fluxes];
%     sB_res.labels = [sB_res.labels labels];
% end
% cd ../../Results/;
% save('sB_res.mat','sB_res');
% cd ../;

%% simulationLowerUptake
% sLU_res.fluxes = zeros(0,0);
% sLU_res.labels = cell(0,0);
% ion_list = {'K' 'MG' 'FE' 'ZN' 'CA' 'MN' 'CU' 'NA'};
% for i = 1:length(ion_list)
%     file_flux = ['sLU_fluxes_',ion_list{i},'.mat'];
%     file_label = ['sLU_labels_',ion_list{i},'.mat'];
%     load(file_flux);
%     load(file_label);
%     sLU_res.fluxes = [sLU_res.fluxes fluxes];
%     sLU_res.labels = [sLU_res.labels labels];
% end
% cd ../../Results/;
% save('sLU_res.mat','sLU_res');
% cd ../;

%% simulationLowerUptakeOCNPS
sLUOCNPS_res.fluxes = zeros(0,0);
sLUOCNPS_res.labels = cell(0,0);
ion_list = {'O' 'C' 'N' 'P' 'S'};
for i = 1:length(ion_list)
    file_flux = ['sLUOCNPS_fluxes_',ion_list{i},'.mat'];
    file_label = ['sLUOCNPS_labels_',ion_list{i},'.mat'];
    load(file_flux);
    load(file_label);
    sLUOCNPS_res.fluxes = [sLUOCNPS_res.fluxes fluxes];
    sLUOCNPS_res.labels = [sLUOCNPS_res.labels labels];
end
cd ../../Results/;
save('sLUOCNPS_res.mat','sLUOCNPS_res');
cd ../;

%% simulationCNPS
% [~,raw,~] = xlsread('Exp_biolog.xlsx');
% metNameID = raw(2:end,2);
% kmax = sum(~ismember(metNameID,''));
% 
% sCNPS.fluxes = zeros(0,0);
% sCNPS.labels = cell(0,0);
% k = 1:1:kmax;
% for i = 1:length(k)
%     display([num2str(i),'/',num2str(length(k))]);
%     file_flux = ['sCNPS_fluxes_',num2str(k(i)),'.mat'];
%     file_label = ['sCNPS_labels_',num2str(k(i)),'.mat'];
%     load(file_flux);
%     load(file_label);
% %     fluxes = fluxes(:,k(i):end);
% %     labels = labels(1,k(i):end);
%     sCNPS.fluxes = [sCNPS.fluxes fluxes];
%     sCNPS.labels = [sCNPS.labels labels];
% end
% cd ../../Results/;
% save('sCNPS_res.mat','sCNPS_res');
% cd ../;

%% simulationAnaerobicCNPS
% [~,raw,~] = xlsread('Exp_biolog.xlsx');
% metNameID = raw(2:end,2);
% kmax = sum(~ismember(metNameID,''));
% 
% sAnaerobicCNPS.fluxes = zeros(0,0);
% sAnaerobicCNPS.labels = cell(0,0);
% k = 1:1:kmax;
% for i = 1:length(k)
%     display([num2str(i),'/',num2str(length(k))]);
%     file_flux = ['sAnaerobicCNPS_fluxes_',num2str(k(i)),'.mat'];
%     file_label = ['sAnaerobicCNPS_labels_',num2str(k(i)),'.mat'];
%     load(file_flux);
%     load(file_label);
% %     fluxes = fluxes(:,k(i):end);
% %     labels = labels(1,k(i):end);
%     sAnaerobicCNPS.fluxes = [sAnaerobicCNPS.fluxes fluxes];
%     sAnaerobicCNPS.labels = [sAnaerobicCNPS.labels labels];
% end
% cd ../../Results/;
% save('sAnaerobicCNPS_res.mat','sAnaerobicCNPS_res');
% cd ../;


