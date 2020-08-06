%% Plot iron
% Sensitivity of each enzyme
load('sI_res.mat');
load('sISP_res.mat');
load('CofactorYeast.mat');

fluxes = sISP_res.fluxes;
label = sISP_res.proteins';
flux_ref = sISP_res.flux_ref;
flux_50uptake = sI_res.fluxes(:,sI_res.k_cf == 0.5);

%% Exchange fluxes
mu = fluxes(strcmp(model.rxns,'r_2111'),:);
mu_ref = flux_ref(strcmp(model.rxns,'r_2111'),:);
mu_50uptake = flux_50uptake(strcmp(model.rxns,'r_2111'),:);

clr = [242,94,13]/255;
[data1,b1] = sort(mu);
lbl1 = label(b1);

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));
[~,b] = ismember(lbl1,gname_1);
lbl1 = gname_2(b);
clear b;

lbl1 = lower(lbl1);
lbl1 = cellfun(@(x) strcat(upper(x(1)),x(2:end)),lbl1,'UniformOutput',false);

figure();
data = [mu_50uptake data1 mu_ref];
lbl = ['50% uptake' lbl1' '100% uptake'];
lbl = strrep(lbl,'_enzyme','');
lbl = cellfun(@(x) strrep(x,'_','-'),lbl,'UniformOutput',false);
b1 = bar(data,0.7,'FaceColor','flat','LineWidth',0.5);
b1.CData(2:end-1,:) = repmat([242,94,13]/255,length(mu),1);
b1.CData(end,:) = [0 0 0];
b1.CData(1,:) = [0.5 0.5 0.5];
b1.EdgeColor = 'w';
xticks(1:1:length(data));
xticklabels(lbl);
xtickangle(90);
ylim([0.29 0.39]);
box off;
set(gca,'XColor','k');
set(gca,'YColor','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[300 300 310 120]);
set(gca,'position',[0.2 0.4 0.7 0.4]);

% figure();
% data = [mu_50uptake data1 mu_ref];
% lbl = ['50% uptake' lbl1' '100% uptake'];
% lbl = strrep(lbl,'_enzyme','');
% lbl = cellfun(@(x) strrep(x,'_','-'),lbl,'UniformOutput',false);
% b1 = barh(data,0.6,'FaceColor','flat','LineWidth',0.5);
% b1.CData(2:end-1,:) = repmat([242,94,13]/255,length(mu),1);
% b1.CData(end,:) = [0 0 0];
% b1.CData(1,:) = [0.5 0.5 0.5];
% b1.EdgeColor = 'w';
% yticks(1:1:length(data));
% yticklabels(lbl);
% xlim([0.29 0.39]);
% set(gca,'XColor','k');
% set(gca,'YColor','k');
% set(gca,'FontSize',6,'FontName','Helvetica');
% xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
% 
% set(gcf,'position',[300 300 120 310]);
% set(gca,'position',[0.4 0.1 0.5 0.6]);


