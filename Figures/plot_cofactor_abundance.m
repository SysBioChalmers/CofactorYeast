%% Plot cofactor abundance
load('CofactorDataset.mat');

[num,txt,~] = xlsread('protein_abundance.xlsx','paxdb');
proteinlist = txt(2:end,2);
abundancelist = num(:,3);
proteincopylist = abundancelist*50*1e6/1e6; %unit:molecule/cell
% 50*1e6 protein molecules per cell in yeast (PMID: 24114984).
clear num txt abundancelist;

%% Exp dataset
exp_cofactor = cell(0,1);
exp_atomcell = zeros(0,1);
exp_source = cell(0,1);

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet1');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet2');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet3');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet4');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet5');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet6');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

[~,~,tmp] = xlsread('Exp_cofactor_abundance.xlsx','Sheet7');
cftrid_tmp = tmp(2:end,ismember(tmp(1,:),'ID'));
atom_tmp = tmp(2:end,ismember(tmp(1,:),'atom/cell'));
atom_tmp = cell2mat(atom_tmp);
idx_tmp =  ~ismember(cftrid_tmp,'-') & ~isnan(atom_tmp);
exp_cofactor = [exp_cofactor; cftrid_tmp(idx_tmp)];
exp_atomcell = [exp_atomcell; atom_tmp(idx_tmp)];
exp_source = [exp_source; repmat(tmp(1,1),length(cftrid_tmp(idx_tmp)),1)];

clear atom_tmp idx_tmp cftrid_tmp tmp;

%% Estimate collected cofactors
cftrlist = unique(exp_cofactor);
atomlist = zeros(length(cftrlist),1);

for j = 1:length(cftrlist)
    
    if strcmp(cftrlist{j},'FE')
        idx_tmp = ismember(CofactorDataset.cofactor,{'FE_III';'FE_II';'HEME_A';'HEME_C';'HEME_B'});
        proteins = CofactorDataset.protein(idx_tmp);
        copies = CofactorDataset.copy(idx_tmp);
        abundance_tmp = 0; % atoms per cell
        for i = 1:length(proteins)
            proteinid = proteins(i);
            if ismember(proteinid,proteinlist)
                abund_tmp = proteincopylist(ismember(proteinlist,proteinid));
                abundance_tmp = abundance_tmp + abund_tmp * copies(i);
            end
        end
        idx_tmp = ismember(CofactorDataset.cofactor,{'ISC_2FE2S'});
        proteins = CofactorDataset.protein(idx_tmp);
        copies = CofactorDataset.copy(idx_tmp) * 2;% two FE per ISC
        for i = 1:length(proteins)
            proteinid = proteins(i);
            if ismember(proteinid,proteinlist)
                abund_tmp = proteincopylist(ismember(proteinlist,proteinid));
                abundance_tmp = abundance_tmp + abund_tmp * copies(i);
            end
        end
        idx_tmp = ismember(CofactorDataset.cofactor,{'ISC_3FE4S'});
        proteins = CofactorDataset.protein(idx_tmp);
        copies = CofactorDataset.copy(idx_tmp) * 3;% three FE per ISC
        for i = 1:length(proteins)
            proteinid = proteins(i);
            if ismember(proteinid,proteinlist)
                abund_tmp = proteincopylist(ismember(proteinlist,proteinid));
                abundance_tmp = abundance_tmp + abund_tmp * copies(i);
            end
        end
        idx_tmp = ismember(CofactorDataset.cofactor,{'ISC_4FE4S'});
        proteins = CofactorDataset.protein(idx_tmp);
        copies = CofactorDataset.copy(idx_tmp) * 4;% four FE per ISC
        for i = 1:length(proteins)
            proteinid = proteins(i);
            if ismember(proteinid,proteinlist)
                abund_tmp = proteincopylist(ismember(proteinlist,proteinid));
                abundance_tmp = abundance_tmp + abund_tmp * copies(i);
            end
        end
    elseif strcmp(cftrlist{j},'CU')
        idx_tmp = ismember(CofactorDataset.cofactor,{'CU_II';'CU_I'});
        proteins = CofactorDataset.protein(idx_tmp);
        copies = CofactorDataset.copy(idx_tmp);
        abundance_tmp = 0; % atoms per cell
        for i = 1:length(proteins)
            proteinid = proteins(i);
            if ismember(proteinid,proteinlist)
                abund_tmp = proteincopylist(ismember(proteinlist,proteinid));
                abundance_tmp = abundance_tmp + abund_tmp * copies(i);
            end
        end
    else
        idx_tmp = ismember(CofactorDataset.cofactor,cftrlist(j));
        proteins = CofactorDataset.protein(idx_tmp);
        copies = CofactorDataset.copy(idx_tmp);

        abundance_tmp = 0; % atoms per cell
        for i = 1:length(proteins)
            proteinid = proteins(i);
            if ismember(proteinid,proteinlist)
                abund_tmp = proteincopylist(ismember(proteinlist,proteinid));
                abundance_tmp = abundance_tmp + abund_tmp * copies(i);
            end
        end
    end
    atomlist(j) = abundance_tmp;
end
clear abund_tmp abundance_tmp copies i j idx_tmp proteinid proteins;

%% Plot

% color
unq_sources = unique(exp_source);
clr_sources = [ 228,26,28
                55,126,184
                77,175,74
                152,78,163
                255,127,0
                255,255,51
                166,86,40
                247,129,191
                153,153,153]/255;
% clr_sources = repmat(linspace(0,200,length(unq_sources))',1,3)/255;

figure();
hold on;
box on;
for i = 1:length(exp_cofactor)
    x_tmp = find(ismember(cftrlist,exp_cofactor{i}));
    y_tmp = exp_atomcell(i);
    clr_tmp = clr_sources(ismember(unq_sources,exp_source{i}),:);
    scatter(x_tmp,log10(y_tmp),30,'o','LineWidth',1,'MarkerEdgeColor',clr_tmp,'MarkerEdgeAlpha',0.75);
end
scatter(1:1:length(cftrlist),log10(atomlist),30,'kx','LineWidth',1);
ylim([4 11]);
xlim([0.1 length(cftrlist)+0.9]);
xticks(1:1:length(cftrlist));
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'XTickLabel',cftrlist);
set(gca,'FontSize',7,'FontName','Helvetica');
ylabel('log10(atoms/cell)','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[200 0 170 180]);
set(gca,'position',[0.13 0.2 0.85 0.6]);


%% Plot modeled and unmodeled cofactors
load('cofactor_info.mat');

figure();
subplot(2,1,1);
c = categorical(cofactor_info.element_id);
c = reordercats(c,cofactor_info.element_id);
b = bar(c,log10(cofactor_info.element_abund_total),'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
b.BarWidth = 0.5;
ylim([4 8]);
set(gca,'FontSize',6,'FontName','Helvetica');
% set(gca,'ycolor','k');
ylabel('log10(abundance/cell)','FontSize',6,'FontName','Helvetica','Color','k');


subplot(2,1,2);
modeled = cofactor_info.element_abund_modeled./cofactor_info.element_abund_total;
unmodeled = 1-modeled;
y = [modeled unmodeled];
bb = bar(c,y,'stacked');
bb(1).FaceColor = [239,138,98]/255;
% bb(1).FaceAlpha = 0.6;
bb(1).BarWidth = 0.5;
bb(2).FaceColor = [128,128,128]/255;
% bb(2).FaceAlpha = 0.6;
bb(2).BarWidth = 0.5;
ylim([0 1]);
set(gca,'FontSize',6,'FontName','Helvetica');
% set(gca,'ycolor','k');
ylabel('Fraction in total abundance','FontSize',6,'FontName','Helvetica','Color','k');

set(gcf,'position',[400 0 200 230]);



