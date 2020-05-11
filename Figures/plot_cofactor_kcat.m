load('enzymedata.mat');
load('CofactorDataset.mat');

kcat_conf_cutoff = 1;

idx = enzymedata.subunit_kcat_conf >= kcat_conf_cutoff;
protein_list = enzymedata.subunit(idx);
kcat_list = enzymedata.subunit_kcat(idx);

tmp_kcat_list = num2cell(kcat_list);
tmp_kcat_list = cellfun(@num2str,tmp_kcat_list,'UniformOutput',false);

tmp_list = strcat(protein_list,';',tmp_kcat_list);
tmp_list = unique(tmp_list);
tmp_list = split(tmp_list,';');

proteins = tmp_list(:,1);
kcats = tmp_list(:,2);
kcats = cell2mat(cellfun(@str2num,kcats,'UniformOutput',false))/3600;% /s

% cofactor_list = {'CA' 'CL' 'COA' 'FAD' 'FMN' 'GSH' 'K' 'MG' 'MN' 'NA'...
%                  'PLP' 'SAM' 'SO4' 'TDP' 'ZN' 'CU' 'FE'};
cofactor_list = {'CA' 'K' 'MG' 'MN' 'NA' 'ZN' 'CU' 'FE'};

figure();
for i = 1:length(cofactor_list)
    cofactorid = cofactor_list{i};
    
    if strcmp(cofactorid,'CU')
        searchid = {'CU_II','CU_I'};
    elseif strcmp(cofactorid,'FE')
        searchid = {'FE_III','FE_II','HEME_A','HEME_C','HEME_B','ISC_2FE2S','ISC_4FE4S','ISC_3FE4S'};
    else
        searchid = {cofactorid};
    end
    
    proteins_tmp = CofactorDataset.protein(ismember(CofactorDataset.cofactor,searchid));
    
    overlap_proteins = proteins(ismember(proteins,proteins_tmp));
    plot_kcats = log10(kcats(ismember(proteins,proteins_tmp)));
    
    subplot(length(cofactor_list),1,i);
    [f,x] = ecdf(plot_kcats);
    hold on;
    box on;
    plot(x,f);
    text(-3,0.8,['N = ',num2str(length(plot_kcats))],'Color','k','FontSize',6,'FontName','Helvetica');
    xlim([-4 6]);
    if i == length(cofactor_list)
        xlabel('log10(kcat)','FontSize',6,'FontName','Helvetica');
    end
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    title(cofactorid,'FontSize',7,'FontName','Helvetica','Color','k');
end
set(gcf,'position',[300 500 50 500]);

