%% simulationCrabtree
% Timing: ~ 500 s
load('CofactorYeast.mat');
load('enzymedata.mat');

tic;

%% Set model
model = changeRxnBounds(model,'r_1714',-1000,'l'); %

% block some reactions
model = changeRxnBounds(model,'r_0886_1',0,'b'); % iso-reaction of PFK
model = changeRxnBounds(model,'r_4262_fwd',0,'b'); % citrate hydroxymutase
model = changeRxnBounds(model,'r_4262_rvs',0,'b'); % citrate hydroxymutase

% block some reactions that done in PMID: 28779005.
model = changeRxnBounds(model,'r_2045_rvs',0,'b'); % serine transport from [m] to [c]
model = changeRxnBounds(model,'r_0659_fwd',0,'b'); % isocitrate dehydrogenase (NADP)
model = changeRxnBounds(model,'r_0659_rvs',0,'b'); % isocitrate dehydrogenase (NADP)

% model = changeRxnBounds(model,'r_0725_fwd',0,'b'); % methenyltetrahydrofolate cyclohydrolase
% model = changeRxnBounds(model,'r_0918',0,'b'); % phosphoserine transaminase


%% Set optimization
rxnID = 'r_1714'; %minimize glucose uptake rate
osenseStr = 'Maximize';

tot_protein = 0.46; %g/gCDW, estimated from the original GEM.
f_modeled_protein = extractModeledprotein(model,'r_4041','s_3717[c]'); %g/gProtein
% r_4041 is pseudo_biomass_rxn_id in the GEM
% s_3717[c] is protein id

f = tot_protein * f_modeled_protein;
clear tot_protein f_modeled_protein;

%% Solve LPs
mu_list = 0.02:0.02:0.38;

fluxes = zeros(length(model.rxns),length(mu_list));

for i = 1:length(mu_list)
    mu = mu_list(i);
    model_tmp = changeRxnBounds(model,'r_2111',mu,'b');
    disp(['mu = ' num2str(mu)]);
    fileName = writeLP(model_tmp,mu,f,osenseStr,rxnID,enzymedata,1);
    command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    [~,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    if strcmp(sol_status,'optimal')
        fluxes(:,i) = sol_full;
    end
end

cd Results/;
save('sC_fluxes.mat','fluxes');
cd ../;
clear;

%% Plot
load('sC_fluxes.mat');
load('CofactorYeast.mat');

mu = fluxes(strcmp(model.rxns,'r_2111'),:);
glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
o2 = -1*fluxes(strcmp(model.rxns,'r_1992'),:);
co2 = fluxes(strcmp(model.rxns,'r_1672'),:);
% ac = fluxes(strcmp(model.rxns,'r_1634'),:);
% aldh = fluxes(strcmp(model.rxns,'r_1631'),:);

figure('Name','orginial');
hold on;
box on;
plot(mu,glc,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[55,126,184]/255);
plot(mu,etoh,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[255,127,0]/255);
plot(mu,o2,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[77,175,74]/255);
plot(mu,co2,'-o','MarkerSize',5,'LineWidth',0.75,'Color',[152,78,163]/255);
% plot(mu,ac,'-o','LineWidth',0.75,'Color',[247,129,191]/255);
% plot(mu,aldh,'-o','LineWidth',0.75,'Color',[247,129,191]/255);
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


toc;

