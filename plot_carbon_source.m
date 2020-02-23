%% Plot carbon source
load('sCS_res.mat');
[num,txt,~] = xlsread('Exp_carbon_sources.xlsx');
exp_cslist = txt(2:end,1);
exp_mulist = num;
clear num txt;

load('CofactorYeast.mat');
load('enzymedata.mat');
model = addReaction(model,'exchange_cd','reactionFormula','s_3783[e] -> ','reversible',true);
model = changeRxnBounds(model,'exchange_cd',-1000,'l');
model = addReaction(model,'transport_cd','reactionFormula','s_3783[e] -> s_3782[c]','reversible',true);

%% plot exp mu v.s. sim mu
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
figure();
hold on;
box on;
line([0 1],[0 1],'Color','k','LineWidth',0.5);
for i = 1:length(exp_mulist)
    x_tmp = exp_mulist(i);
    y_tmp = sim_mulist(i);
    clr_tmp = clr_cs(ismember(unq_cs,exp_cslist{i}),:);
    scatter(x_tmp,y_tmp,20,'o','filled','LineWidth',1,'MarkerEdgeColor',clr_tmp,'MarkerFaceColor',clr_tmp,'MarkerFaceAlpha',0.5);
end
xlim([0.01 0.49]);
ylim([0.01 0.49]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Simulated growth rate (/h)','FontSize',7,'FontName','Helvetica');
xlabel('Measured growth rate (/h)','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[200 200 150 150]);
set(gca,'position',[0.2 0.2 0.6 0.6]);

%% plot change in cofactor

glc_idx = ismember(sCS_res.cslist,'Glucose');
flux_glc = sCS_res.fluxes(:,glc_idx);
[list_glc_tmp, concentration_glc_tmp] = simulatedCofactor(model,enzymedata,flux_glc);
[concentration_glc, I_tmp] = sort(concentration_glc_tmp);
list_glc = list_glc_tmp(I_tmp);

other_cs_tmp = sCS_res.cslist(~glc_idx);
other_mu_tmp = sCS_res.mulist(~glc_idx);
[other_mu, I_tmp] = sort(other_mu_tmp);
other_cs = other_cs_tmp(I_tmp);

data = zeros(length(list_glc),length(other_cs));

for i = 1:length(other_cs)
    flux_tmp = sCS_res.fluxes(:,ismember(sCS_res.cslist,other_cs(i)));
    [list_tmp, concentration_tmp] = simulatedCofactor(model,enzymedata,flux_tmp);
    for j = 1:length(list_glc)
        if ismember(list_glc(j),list_tmp)
            data(j,i) = concentration_tmp(ismember(list_tmp,list_glc(j)));
        else
            data(j,i) = nan;
        end
    end
end

data = log2(data./concentration_glc);
x = other_cs;
[~,Locb] = ismember(list_glc,model.mets);
y = model.metNames(Locb);
y = cellfun(@(x) strrep(x,' [cytoplasm]',''),y,'UniformOutput',false);

data = data([1:5,7:end],:);
y = y([1:5,7:end],:);

maxclr = [103,0,31]/255;
mdlclr = [255,255,255]/255;
minclr = [5,48,97]/255;
tmp11 = linspace(minclr(1),mdlclr(1),128)';
tmp12 = linspace(mdlclr(1),maxclr(1),128)';
tmp21 = linspace(minclr(2),mdlclr(2),128)';
tmp22 = linspace(mdlclr(2),maxclr(2),128)';
tmp31 = linspace(minclr(3),mdlclr(3),128)';
tmp32 = linspace(mdlclr(3),maxclr(3),128)';
tmp1 = [tmp11;tmp12(2:end)];
tmp2 = [tmp21;tmp22(2:end)];
tmp3 = [tmp31;tmp32(2:end)];
clrmap = [tmp1 tmp2 tmp3];

maxvalue = max(max(data));
minvalue = min(min(data));
minposition = 1;
mdlposition = find(any((clrmap == mdlclr)'));
maxposition = length(clrmap);
if abs(minvalue) >= abs(maxvalue)
    startposition = minposition;
    endposition = round(abs(maxvalue/minvalue)*(maxposition-mdlposition))+mdlposition;
else
    endposition = maxposition;
    startposition = mdlposition-floor(abs(minvalue/maxvalue)*mdlposition);
end

clrmap = clrmap(startposition:endposition,:);


figure();
h = heatmap(x,y,data,'Colormap',clrmap,'ColorMethod','count','CellLabelColor','none');
h.Title = 'log2FC';
set(h,'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[500 300 200 300]);
set(gca,'position',[0.45 0.3 0.35 0.5]);



