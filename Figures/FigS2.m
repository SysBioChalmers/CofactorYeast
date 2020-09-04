%% Max growth rates and Crabtree

load('CofactorYeast.mat');

%% plot max growth rates
load('sCS_res.mat');
[num,txt,~] = xlsread('Exp_carbon_sources.xlsx');
exp_cslist = txt(2:end,1);
exp_mulist = num;
clear num txt;
sim_mulist = zeros(length(exp_mulist),1);
for i = 1:length(exp_cslist)
    sim_mulist(i) = sCS_res.mulist(ismember(sCS_res.cslist,exp_cslist(i)));
end

unq_cs = unique(exp_cslist);
clr_cs = [  228,26,28;
            55,126,184;
            77,175,74;
            152,78,163;
            255,127,0;
            247,129,191;
            166,86,40;
            153,153,153]/255;
figure('Name','1');
hold on;
box on;
line([0 1],[0 1],'Color','k','LineWidth',0.5);
cslist_tmp = {'none'};
h_label_top = 0.45;
h_label_btm = 0.05;
for i = 1:length(exp_mulist)
    x_tmp = exp_mulist(i);
    y_tmp = sim_mulist(i);
    clr_tmp = clr_cs(ismember(unq_cs,exp_cslist{i}),:);
    scatter(x_tmp,y_tmp,20,'o','filled','LineWidth',1,'MarkerEdgeColor',clr_tmp,'MarkerFaceColor',clr_tmp,'MarkerFaceAlpha',0.5);
    if ~contains(cslist_tmp,exp_cslist(i))
        if length(cslist_tmp) < 5
            text(0.05,h_label_top,exp_cslist{i},'FontSize',6,'FontName','Helvetica','Color',clr_tmp);
            h_label_top = h_label_top - 0.05;
        else
            text(0.3,h_label_btm,exp_cslist{i},'FontSize',6,'FontName','Helvetica','Color',clr_tmp);
            h_label_btm = h_label_btm + 0.05;
        end
        cslist_tmp = [cslist_tmp;exp_cslist(i)];
    end
    
end
xlim([0.01 0.49]);
ylim([0.01 0.49]);

% rmse = sqrt(sum((exp_mulist-sim_mulist).^2)/numel(exp_mulist));
% 
% r2txt = ['R^2 = ',num2str(r2)];
% text(0.05,0.35,r2txt,'FontSize',6,'FontName','Helvetica','Color','k');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Simulated growth rate (/h)','FontSize',7,'FontName','Helvetica');
xlabel('Measured growth rate (/h)','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[200 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);


%% plot Crabtree
load('sC_fluxes.mat');

mu = fluxes(strcmp(model.rxns,'r_2111'),:);
glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
o2 = -1*fluxes(strcmp(model.rxns,'r_1992'),:);
co2 = fluxes(strcmp(model.rxns,'r_1672'),:);

% Yeast exp data (PMID: 9603825)
fluxes_exp_yeast = [0.1  0.15  0.2  0.25  0.28  0.3   0.32  0.35  0.36  0.38 ; % mu
                    1.1  1.67  2.15 2.83  3.24  3.7   5.44  8.09  8.33  10.23; % glucose
                    0    0     0    0     0     0.51  4.42  6.91  6.71  14.91; % ethanol
                    2.73 2.5   5.07 6.8   8.3   8.8   6.83  6.6   7.1   4.19]; % o2

figure('Name','2');
hold on;
box on;

plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(2,:),'o','LineWidth',0.75,'Color',[55,126,184]/255,'MarkerSize',5);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(3,:),'o','LineWidth',0.75,'Color',[255,127,0]/255,'MarkerSize',5);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(4,:),'o','LineWidth',0.75,'Color',[77,175,74]/255,'MarkerSize',5);
plot(mu,glc,'-','LineWidth',0.75,'Color',[55,126,184]/255);
plot(mu,etoh,'-','LineWidth',0.75,'Color',[255,127,0]/255);
plot(mu,o2,'-','LineWidth',0.75,'Color',[77,175,74]/255);

xlim([0 0.4]);

set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
ylabel('Flux (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica','Color','k');
legend({'Glucose uptake',...
        'Ethanol production',...
        'O2 uptake'},'FontSize',6,'FontName','Helvetica','location','nw');
set(gcf,'position',[0 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);