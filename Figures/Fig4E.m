load('sIRCAA_res.mat');
load('modelNoscapine.mat');
fluxes = sIRCAA_res.fluxes;
AA_Exchanges = {'r_1873' ... % L-alanine exchange
                'r_1879' ... % L-arginine exchange
                'r_1880' ... % L-asparagine exchange
                'r_1881' ... % L-aspartate exchange
                'r_1883' ... %	L-cysteine exchange
                'r_1889' ... %	L-glutamate exchange
                'r_1891' ... %	L-glutamine exchange
                'r_1810' ... %	L-glycine exchange
                'r_1893' ... %	L-histidine exchange
                'r_1897' ... %	L-isoleucine exchange
                'r_1899' ... %	L-leucine exchange
                'r_1900' ... %	L-lysine exchange
                'r_1902' ... %	L-methionine exchange
                'r_1903' ... %	L-phenylalanine exchange
                'r_1904' ... %	L-proline exchange
                'r_1906' ... %	L-serine exchange
                'r_1911' ... %	L-threonine exchange
                'r_1912' ... %	L-tryptophan exchange
                'r_1913' ... %	L-tyrosine exchange
                'r_1914'};    %	L-valine exchange
AA_IDs =   {'Ala' 'Arg' 'Asn' 'Asp' 'Cys' 'Glu' 'Gln' 'Gly' 'His' 'Ile' ...
            'Leu' 'Lys' 'Met' 'Phe' 'Pro' 'Ser' 'Thr' 'Trp' 'Tyr' 'Val'};


mu = fluxes(ismember(model.rxns,'r_2111'),:);
delta_mu = mu(2:end)-mu(1);
delta_aa = zeros(1,length(AA_IDs));
for i = 1:length(AA_IDs)
    fluxes_tmp = fluxes(:,ismember(sIRCAA_res.lables,AA_IDs(i)));
    delta_aa(1,i) = abs(fluxes_tmp(ismember(model.rxns,AA_Exchanges{i}),:));
end
data = delta_mu./delta_aa;

minclr = [255,255,178]/255;
maxclr = [242,94,13]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];


figure('Name','1');
h = heatmap(sIRCAA_res.lables(2:end),{'50% iron uptake'},data,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','none');
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[100 700 300 50]);
set(gca,'position',[0.25 0.4 0.45 0.15]);

