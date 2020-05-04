%% Plot ion uptake
load('sIU_res.mat');

load('CofactorYeast.mat');
load('enzymedata.mat');

% load reference
load('sC_fluxes.mat');
mu_ref_list = fluxes(strcmp(model.rxns,'r_2111'),:);
glc_ref_list = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
Ybio_ref_list = mu_ref_list./glc_ref_list;

%% Setting
decrease_value = 0.8:-0.2:0; %should be the same as "simulationIonUptake.m"
steps = length(decrease_value);
whole_value = [1 decrease_value];

%% Fluxes
figure('Name','1');
for i = 1:length(sIU_res.ionlist)
    ionid = sIU_res.ionlist{i};
    fluxes_changed = sIU_res.fluxes(:,(i-1)*steps+1:(i-1)*steps+steps);
    fluxes_whole = [sIU_res.fluxes_ref fluxes_changed];
    
    mu = fluxes_whole(strcmp(model.rxns,'r_2111'),:);
    glc = -1*fluxes_whole(strcmp(model.rxns,'r_1714'),:);
    etoh = fluxes_whole(strcmp(model.rxns,'r_1761'),:);
    o2 = -1*fluxes_whole(strcmp(model.rxns,'r_1992'),:);
    co2 = fluxes_whole(strcmp(model.rxns,'r_1672'),:);
    Ybio = mu./glc;
    
    subplot(4,length(sIU_res.ionlist),i);
    hold on;
    box on;
    plot(whole_value,mu,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[140,81,10]/255);
    xlim([0 1]);
    ylim([0 0.5]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Growth (/h)','FontSize',6,'FontName','Helvetica');
    end
    title(ionid,'FontSize',7,'FontName','Helvetica','Color','k');
    
    subplot(4,length(sIU_res.ionlist),i+length(sIU_res.ionlist));
    hold on;
    box on;
    plot(whole_value,glc,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[55,126,184]/255);
    plot(whole_value,etoh,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[255,127,0]/255);
    plot(whole_value,o2,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[77,175,74]/255);
    plot(whole_value,co2,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[152,78,163]/255);
    xlim([0 1]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Exchange fluxes','FontSize',6,'FontName','Helvetica');
    end
    
    subplot(4,length(sIU_res.ionlist),i+2*length(sIU_res.ionlist));
    hold on;
    box on;
    plot(whole_value,Ybio,'-o','MarkerSize',2,'LineWidth',0.75,'Color',[1,102,94]/255);
    xlim([0 1]);
    ylim([0 0.1]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Biomass yield','FontSize',6,'FontName','Helvetica');
    end
    xlabel('Relative uptake','FontSize',6,'FontName','Helvetica');
    
    subplot(4,length(sIU_res.ionlist),i+3*length(sIU_res.ionlist));
    hold on;
    box on;
    plot(mu,Ybio,'-o','MarkerSize',2,'LineWidth',0.75,'Color','r');
    plot(mu_ref_list,Ybio_ref_list,'-','LineWidth',0.75,'Color','k');
    xlim([0 0.5]);
    ylim([0 0.1]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    if i == 1
        ylabel('Biomass yield','FontSize',6,'FontName','Helvetica');
    end
    xlabel('Growth (/h)','FontSize',6,'FontName','Helvetica');
end
set(gcf,'position',[300 500 550 250]);


