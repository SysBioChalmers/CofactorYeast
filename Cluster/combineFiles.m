%% combineFiles
cd tmp_results/;
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


sLU_res.fluxes = zeros(0,0);
sLU_res.labels = cell(0,0);
ion_list = {'K' 'MG' 'FE' 'ZN' 'CA' 'MN' 'CU' 'NA'};
for i = 1:length(ion_list)
    file_flux = ['sLU_fluxes_',ion_list{i},'.mat'];
    file_label = ['sLU_labels_',ion_list{i},'.mat'];
    load(file_flux);
    load(file_label);
    sLU_res.fluxes = [sLU_res.fluxes fluxes];
    sLU_res.labels = [sLU_res.labels labels];
end
cd ../../Results/;
save('sLU_res.mat','sLU_res');
cd ../;