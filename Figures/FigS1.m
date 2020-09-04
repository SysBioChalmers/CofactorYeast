%% Cofactor Dataset

load('CofactorDataset.mat');
load('CofactorYeast.mat');

ion_id_list = {'CA';'CU';'FE';'K';'MG';'MN';'NA';'ZN'};


%%

tot_num = zeros(length(ion_id_list),1);
modeled_num = zeros(length(ion_id_list),1);

for i = 1:length(ion_id_list)
    cft_tmp = ion_id_list{i};
    if strcmp(cft_tmp,'CU')
        cft_tmp = {'CU_I';'CU_II'};
    elseif strcmp(cft_tmp,'FE')
        cft_tmp = {'FE_III';'FE_II';'HEME_A';'HEME_C';'HEME_B';'ISC_2FE2S';'ISC_3FE4S';'ISC_4FE4S'};
    else
        cft_tmp = {cft_tmp};
    end
    tot_tmp = CofactorDataset.protein(ismember(CofactorDataset.cofactor,cft_tmp));
    tot_num(i,1) = length(tot_tmp);
    modeled_tmp = model.genes(ismember(model.genes,tot_tmp));
    modeled_num(i,1) = length(modeled_tmp);
end

figure();
clr1 = [77,77,77]/255;
[data1,b1] = sort(tot_num);
lbl1 = ion_id_list(b1);
subplot(1,2,1);
bar(data1,0.5,'FaceColor',clr1,'FaceAlpha',0.5,'EdgeColor',clr1,'LineWidth',0.8);
xticks(1:1:length(data1));
xticklabels(lbl1);
ylim([0 800]);
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Number of proteins','FontSize',7,'FontName','Helvetica','Color','k');
title('Total proteins (N = 6002)','FontSize',7,'FontName','Helvetica','Color','k');
% total number of proteins is from https://www.ncbi.nlm.nih.gov/genome/15?genome_assembly_id=22535
for i = 1:length(data1)
    if i > 5
        perc = round(data1(i)/6002*100);
        if i < 8
        text(i-0.35,data1(i)+60,num2str([num2str(perc) '%']),'FontSize',7,'FontName','Helvetica','Color',clr1);
        else
        text(i-0.5,data1(i)+60,num2str([num2str(perc) '%']),'FontSize',7,'FontName','Helvetica','Color',clr1);
        end
    end
end

[data2,b2] = sort(modeled_num);
clr2 = [178,24,43]/255;
lbl2 = ion_id_list(b2);
subplot(1,2,2);
bar(data2,0.5,'FaceColor',clr2,'FaceAlpha',0.5,'EdgeColor',clr2,'LineWidth',0.8);
xticks(1:1:length(data2));
xticklabels(lbl2);
ylim([0 200]);
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Number of proteins','FontSize',7,'FontName','Helvetica','Color','k');
title(['Metabolic proteins (N = ' num2str(length(model.genes)) ')'],'FontSize',7,'FontName','Helvetica','Color','k');
for i = 1:length(data2)
    if i > 5
        perc = round(data2(i)/length(model.genes)*100);
        if i < 8
        text(i-0.35,data2(i)+15,num2str([num2str(perc) '%']),'FontSize',7,'FontName','Helvetica','Color',clr2);
        else
        text(i-0.5,data2(i)+15,num2str([num2str(perc) '%']),'FontSize',7,'FontName','Helvetica','Color',clr2);
        end
    end
end
set(gcf,'position',[300 300 400 130]);









