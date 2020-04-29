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
clear;

