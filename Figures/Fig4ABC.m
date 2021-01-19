%% Plot iron

load('sI_res.mat');
load('CofactorYeast.mat');
load('enzymedata.mat');
load('CofactorDataset.mat');

selected_data = 1:1:10;
fluxes = [sI_res.fluxes(:,selected_data) sI_res.flux_ref];
label_tmp = num2str(sI_res.k_cf(1,selected_data));
label_tmp = strsplit(label_tmp);
label = [label_tmp 'Ref'];

%% Exchange fluxes
glc = -1*fluxes(strcmp(model.rxns,'r_1714'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
mu = fluxes(strcmp(model.rxns,'r_2111'),:);
fe = -1*fluxes(strcmp(model.rxns,'r_1861'),:);
glyc = fluxes(strcmp(model.rxns,'r_1808'),:);

figure('Name','1');
subplot(4,1,1);
b1 = bar(mu,0.7,'FaceColor','flat','LineWidth',0.1);
b1.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b1.CData(end,:) = [0 0 0];
b1.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([0.25 0.42]);
yticks(0.25:0.05:0.4);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica','Color','k');
subplot(4,1,2);
b2 = bar(glc,0.7,'FaceColor','flat','LineWidth',0.1);
b2.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b2.CData(end,:) = [0 0 0];
b2.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([5 25]);
yticks(5:10:25);
set(gca,'FontSize',6,'FontName','Helvetica');
title('Glucose uptake','FontSize',7,'FontName','Helvetica','Color','k');
subplot(4,1,3);
b3 = bar(etoh,0.7,'FaceColor','flat','LineWidth',0.1);
b3.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b3.CData(end,:) = [0 0 0];
b3.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([15 35]);
yticks(15:10:35);
set(gca,'FontSize',6,'FontName','Helvetica');
title('Ethanol production','FontSize',7,'FontName','Helvetica','Color','k');
ylabel('Flux (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica','Color','k');
subplot(4,1,4);
b4 = bar(glyc,0.7,'FaceColor','flat','LineWidth',0.1);
b4.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b4.CData(end,:) = [0 0 0];
b4.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([0 5]);
set(gca,'FontSize',6,'FontName','Helvetica');
title('Glycerol production','FontSize',7,'FontName','Helvetica','Color','k');
xlabel('theta','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[100 300 120 300]);

%% Protein levels
complexes = {'cytochrome c oxidase';'ubiquinol cytochrome-c reductase';'succinate dehydrogenase'};
complex_rxns = {'r_0438';'r_0439';'r_1021'};
proteins = {'ACO1';'HEM13';'HEM15';'ILV3';'LEU1';'LYS4';'MET5';'ERG1';'ERG11';'ERG25';'OLE1'};

protein_conc_ref = calculateProteinConc(model,model.genes,fluxes(:,end));
protein_conc_low = calculateProteinConc(model,model.genes,fluxes(:,1:end-1));
protein_conc = [protein_conc_low protein_conc_ref];
% remove low absolute protein level in reference
cutoff_low_abs = 0.05;
low_abs_value = quantile(protein_conc_ref(protein_conc_ref>0),cutoff_low_abs);
protein_list = model.genes(protein_conc_ref>low_abs_value);
data_abs = protein_conc(protein_conc_ref>low_abs_value,:);
data_rel = data_abs(:,1:end-1)./data_abs(:,end);

[~,txt,~] = xlsread('Yeast8_Modification.xlsx','SGDgeneNames');
gname_1 = txt(2:end,1);
gname_2 = txt(2:end,2);
clear txt;
gname_2(ismember(gname_2,'')) = gname_1(ismember(gname_2,''));
[~,b] = ismember(protein_list,gname_1);
protein_id_list = gname_2(b);
clear b;

% heatmap
heatmapdata = zeros(length(complexes)+length(proteins),length(label)-1);
heatmaplabel = cell(length(complexes)+length(proteins),1);
for i = 1:length(complexes)+length(proteins)
    if i <= length(complexes) % for complexes, expression change is approximated as metabolic flux change
        lbl = complexes(i);
        metrxn = complex_rxns(i);
        flux_tmp = sum(fluxes(startsWith(model.rxns,metrxn)&~contains(model.rxns,'enzyme'),:));
        heatmapdata(i,:) = flux_tmp(1,1:end-1)/flux_tmp(1,end);
    else
        lbl = proteins(i-length(complexes));
        heatmapdata(i,:) = data_rel(ismember(protein_id_list,lbl),:);
    end
    lbl = lower(lbl);
    lbl = cell2mat(lbl);
    lbl = strcat(upper(lbl(1)),lbl(2:end));
    heatmaplabel(i) = {lbl};
end

heatmapdata(heatmapdata == 0) = 1e-6;
heatmapdata = log2(heatmapdata);

heatmapdata(heatmapdata < -2) = -2;
heatmapdata(heatmapdata > 2) = 2;

minclr = [5,113,176]/255;
mdlclr = [255,255,255]/255;
maxclr = [202,0,32]/255;
tmp1 = [linspace(minclr(1),mdlclr(1),180) linspace(mdlclr(1),maxclr(1),180)]';
tmp2 = [linspace(minclr(2),mdlclr(2),180) linspace(mdlclr(2),maxclr(2),180)]';
tmp3 = [linspace(minclr(3),mdlclr(3),180) linspace(mdlclr(3),maxclr(3),180)]';
clrmap = [tmp1 tmp2 tmp3];

figure('Name','heatmap1');
for i = 1:length(heatmaplabel)
    subplot(length(heatmaplabel),1,i);
    data_tmp = heatmapdata(i,:);
    max_tmp = max(data_tmp);
    min_tmp = min(data_tmp);
    

    if max_tmp<0
        maxpos = round(((-2-max_tmp)/-2)*180,0);
        minpos = round(((-2-min_tmp)/-2)*180,0);
    elseif max_tmp>0 && min_tmp>0
        maxpos = round(max_tmp/2*180,0)+180;
        minpos = round(min_tmp/2*180,0)+180;
    elseif max_tmp>0 && min_tmp<0
        maxpos = round(max_tmp/2*180,0)+180;
        minpos = round(((-2-min_tmp)/-2)*180,0);
    end
    
    if maxpos==0
        maxpos = maxpos+1;
    end
    if minpos==0
        minpos = minpos+1;
    end
    
    clrmap_tmp = clrmap(minpos:maxpos,:);
    h_tmp = heatmap(label_tmp,heatmaplabel(i),data_tmp,'Colormap',clrmap_tmp,...
        'ColorMethod','count','CellLabelColor','none');
    h_tmp.XLabel = 'theta';
    set(h_tmp,'FontSize',6,'FontName','Helvetica');
    h_tmp.FontColor = 'k';
end
set(gcf,'position',[300 0 220 490]);


figure('Name','heatmap2');
h = heatmap(label_tmp,heatmaplabel,heatmapdata,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','none');
h.XLabel = 'theta';
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[600 300 250 300]);
set(gca,'position',[0.5 0.4 0.3 0.4]);


%% Usage for Erg25 and Ole1
cfusage_Erg25 = calculateCofactorUsage4protein(model,'FE',{'YGR060W'},fluxes);
cfusage_Ole1 = calculateCofactorUsage4protein(model,'FE',{'YGL055W'},fluxes);
tot = fe./ mu;
perc_Erg25 = cfusage_Erg25./tot*100;
perc_Ole1 = cfusage_Ole1./tot*100;

figure('Name','C');
subplot(2,1,1);
b1 = bar(perc_Erg25,0.7,'FaceColor','flat','LineWidth',0.1);
b1.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b1.CData(end,:) = [0 0 0];
b1.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([0 2.5]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Iron usage fraction (%)','FontSize',7,'FontName','Helvetica','Color','k');
subplot(2,1,2);
b1 = bar(perc_Ole1,0.7,'FaceColor','flat','LineWidth',0.1);
b1.CData(1:length(selected_data),:) = repmat([242,94,13]/255,length(selected_data),1);
b1.CData(end,:) = [0 0 0];
b1.EdgeColor = 'w';
xticks(1:1:length(label)+1);
xticklabels(label);
xtickangle(90);
ylim([0 19]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('theta','FontSize',7,'FontName','Helvetica','Color','k');

set(gcf,'position',[100 300 110 80]);





