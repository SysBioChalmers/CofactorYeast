load('sN3_res.mat');
load('modelNoscapine.mat');


minclr = [255,255,178]/255;
maxclr = [240,59,32]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];


figure('Name','1');
flux1 = sN3_res.fluxes_UnlimitedFe;
Noscapine1 = flux1(strcmp(model.rxns,'new_r_eNoscapine'),:);
data1 = (Noscapine1(2:end)-Noscapine1(1))/0.01;

flux2 = sN3_res.fluxes_LimitedFe;
Noscapine2 = flux2(strcmp(model.rxns,'new_r_eNoscapine'),:);
data2 = (Noscapine2(2:end)-Noscapine2(1))/0.01;

data = [data1;data2];
data(data > 0.2) = 0.2;

h = heatmap(sN3_res.lables(2:end),{'Unlimited iron uptake';'50% iron uptake'},data,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','none');
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[100 700 300 100]);
set(gca,'position',[0.25 0.4 0.45 0.16]);

