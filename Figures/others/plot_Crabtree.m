% Plot Crabtree


%% Plot
load('sC_fluxes.mat');
load('CofactorYeast.mat');

mu = fluxes(strcmp(model.rxns,'r_2111'),:);
glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
o2 = -1*fluxes(strcmp(model.rxns,'r_1992'),:);
co2 = fluxes(strcmp(model.rxns,'r_1672'),:);
ac = fluxes(strcmp(model.rxns,'r_1634'),:);
aldh = fluxes(strcmp(model.rxns,'r_1631'),:);

figure('Name','orginial');
hold on;
box on;
plot(mu,glc,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[55,126,184]/255);
plot(mu,etoh,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[255,127,0]/255);
plot(mu,o2,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[77,175,74]/255);
plot(mu,co2,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[152,78,163]/255);
plot(mu,ac,'-o','LineWidth',0.75,'Color',[247,129,191]/255);
plot(mu,aldh,'-o','LineWidth',0.75,'Color',[247,129,191]/255);
xlim([0 0.4]);

set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel('Flux (mmol/gCDW/h)','FontSize',14,'FontName','Helvetica');
legend({'Glucose uptake',...
        'Ethanol production',...
        'O2 uptake',...
        'CO2 production'},'FontSize',12,'FontName','Helvetica','location','nw');
set(gcf,'position',[0 400 240 185]);
set(gca,'position',[0.17 0.2 0.76 0.75]);


dummy = fluxes(strcmp(model.rxns,'dilute_dummy'),:)./mu; %g/gCDW
figure('Name','dummy');
hold on;
box on;
plot(mu,dummy,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[82,82,82]/255);
xlim([0 0.4]);
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',14,'FontName','Helvetica');
ylabel('Dummy (g/gCDW)','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[0 0 240 140]);
set(gca,'position',[0.17 0.28 0.76 0.63]);

% Plot ions
ion_id_list = {'K';'MG';'FE';'ZN';'CA';'MN';'CU';'NA'};
ion_ex_list = {'r_2020'; ... % potassium exchange
               'r_4597'; ... % Mg(2+) exchange
               'r_1861'; ... % iron(2+) exchange
               'r_4596'; ... % Zn(2+) exchange
               'r_4600'; ... % Ca(2+) exchange
               'r_4595'; ... % Mn(2+) exchange
               'r_4594'; ... % Cu(2+) exchange
               'r_2049'}; ... % sodium exchange

figure('Name','ion');
for i = 1:length(ion_ex_list)
    ion_id = ion_id_list{i};
    ion_exR = ion_ex_list{i};
    ion_exR_flux = fluxes(ismember(model.rxns,ion_exR),:);
    ion_conc = -1 * ion_exR_flux ./ mu;
    ion_num = ion_conc*0.001*6.02e23*13e-12;
    
    subplot(length(ion_id_list)/2,2,i);
    hold on;
    box on;
    plot(mu,log10(ion_num));
    
    xlim([0 0.42]);
    ylim([5 8]);
    set(gca, 'XColor','k');
    set(gca, 'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    xlabel('Growth','FontSize',6,'FontName','Helvetica');
	ylabel('log10(atoms/cell)','FontSize',6,'FontName','Helvetica');
    title(ion_id,'FontSize',7,'FontName','Helvetica','Color','k');
end
set(gcf,'position',[800 500 180 350]);


% figure('Name','new');
% lower_values = glc/glc(end);
% yield = mu./(glc*180/1000);
% 
% subplot(3,1,1);
% hold on;
% box on;
% plot(lower_values,mu,'-o','MarkerSize',1,'LineWidth',0.75,'Color',[79,89,109]/255);
% xlim([0 1]);
% ylim([0 0.45]);
% set(gca, 'XColor','k');
% set(gca, 'YColor','k');
% set(gca,'FontSize',6,'FontName','Helvetica');
% if i == 1
%     ylabel('Growth (/h)','FontSize',6,'FontName','Helvetica');
% end
% title('Glcucose','FontSize',7,'FontName','Helvetica','Color','k');
% 
% subplot(3,1,2);
% hold on;
% box on;
% plot(lower_values,glc,'-o','MarkerSize',1,'LineWidth',0.75,'Color',[27,158,119]/255);
% plot(lower_values,etoh,'-o','MarkerSize',1,'LineWidth',0.75,'Color',[217,95,2]/255);
% 
% xlim([0 1]);
% ylim([0 40]);
% set(gca, 'XColor','k');
% set(gca, 'YColor','k');
% set(gca,'FontSize',6,'FontName','Helvetica');
% if i == 1
%     ylabel('Exchange rate','FontSize',6,'FontName','Helvetica');
% end
% 
% subplot(3,1,3);
% hold on;
% box on;
% plot(lower_values,yield,'-o','MarkerSize',1,'LineWidth',0.75,'Color',[79,89,109]/255);
% xlim([0 1]);
% ylim([0 0.55]);
% set(gca, 'XColor','k');
% set(gca, 'YColor','k');
% set(gca,'FontSize',6,'FontName','Helvetica');
% if i == 1
%     ylabel('Biomass yield','FontSize',6,'FontName','Helvetica');
% end
% xlabel('Relative uptake','FontSize',6,'FontName','Helvetica');
% set(gcf,'position',[300 500 80 180]);
