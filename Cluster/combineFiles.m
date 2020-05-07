%% combineFiles
cd tmp_results/;
sB_res.fluxes = zeros(0,0);
sB_res.labels = cell(0,0);
k = 1:1:116;
for i = 1:length(k)
    display([num2str(i),'/',num2str(length(k))]);
    file_flux = ['sB_fluxes_',num2str(k(i)),'.mat'];
    file_label = ['sB_labels_',num2str(k(i)),'.mat'];
    load(file_flux);
    load(file_label);
%     fluxes = fluxes(:,k(i):end);
%     labels = labels(1,k(i):end);
    sB_res.fluxes = [sB_res.fluxes fluxes];
    sB_res.labels = [sB_res.labels labels];
end
cd ../../Results/;
save('sB_res.mat','sB_res');
cd ../;
