load('CofactorDataset.mat');
load('CofactorYeast.mat');

unqcofactor = unique(CofactorDataset.cofactor);

tot_num = zeros(length(unqcofactor),1);
modeled_num = zeros(length(unqcofactor),1);

for i = 1:length(unqcofactor)
    cft_tmp = unqcofactor(i);
    tot_tmp = CofactorDataset.protein(ismember(CofactorDataset.cofactor,unqcofactor(i)));
    tot_num(i,1) = length(tot_tmp);
    modeled_tmp = model.genes(ismember(model.genes,tot_tmp));
    modeled_num(i,1) = length(modeled_tmp);
end

figure();
unqcofactorid = cellfun(@(x) strrep(x,'_','-'),unqcofactor,'UniformOutput',false);
clr = [8,81,156]/255;
[data1,b1] = sort(tot_num);
lbl1 = unqcofactorid(b1);
subplot(1,2,1);
barh(data1,0.5,'FaceColor',clr,'FaceAlpha',0.3,'EdgeColor',clr,'LineWidth',0.5);
yticks(1:1:length(data1));
yticklabels(lbl1);
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Number of proteins','FontSize',6,'FontName','Helvetica','Color','k');
title('Total proteins','FontSize',7,'FontName','Helvetica','Color','k');
for i = 1:length(data1)
    text(data1(i)+20,i,num2str(data1(i)),'FontSize',6,'FontName','Helvetica','Color',clr);
end

[data2,b2] = sort(modeled_num);
lbl2 = unqcofactorid(b2);
subplot(1,2,2);
barh(data2,0.5,'FaceColor',clr,'FaceAlpha',0.3,'EdgeColor',clr,'LineWidth',0.5);
yticks(1:1:length(data2));
yticklabels(lbl2);
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Number of proteins','FontSize',6,'FontName','Helvetica','Color','k');
title('Metabolic proteins','FontSize',7,'FontName','Helvetica','Color','k');
for i = 1:length(data2)
    text(data2(i)+5,i,num2str(data2(i)),'FontSize',6,'FontName','Helvetica','Color',clr);
end
set(gcf,'position',[300 500 300 380]);









