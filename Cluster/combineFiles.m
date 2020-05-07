%% combineFiles
cd Batch/;
sGK_res.fluxes = zeros(0,0);
sGK_res.genes = cell(0,0);
k = 1:6:1147;
for i = 1:length(k)
    display([num2str(i),'/',num2str(length(k))]);
    file_flux = ['sGK_fluxes_',num2str(k(i)),'.mat'];
    file_gene = ['sGK_genes_',num2str(k(i)),'.mat'];
    load(file_flux);
    load(file_gene);
    fluxes = fluxes(:,k(i):end);
    genes = genes(1,k(i):end);
    sGK_res.fluxes = [sGK_res.fluxes fluxes];
    sGK_res.genes = [sGK_res.genes genes];
end
cd ../Results/;
save('sGK_res.mat','sGK_res');
cd ../;
