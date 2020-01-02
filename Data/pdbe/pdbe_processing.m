%% pdbe_processer
%  Process data obtained from pdbe.

% The pdbe dataset was obtained 2019.12.30.

% 1. Different pdb ids even with the same corresponding gene id could have 
% different cofactor and protein stoichiometry information. Therefore, for
% the gene id with different pdb ids, the pdb id with the best blastp score
% would be selected. If multiple pdb ids have the same highest blastp
% score, then they would be all selected. If they have the same protein
% stoichiometry then the cofactor information could be integrated. If not,
% they should be manually checked.
% 2. Given that the collected pdbe information do not include all yeast
% genes, especially genes collected in the dataset but their homologue
% genes not. Therefore, yeast homologue genes information should be used to
% assign pdbe information for the genes that are not in the collected
% dataset but their homologues are.
% 3. The remaining genes are assumed to be monomer with no factors bound.

tic;

%% Load data

% load pdbe dataset
load('pdb.mat');%obtained 2019.12.30

% load blastp result
[num,txt,~] = xlsread('blastp_res.xlsx');%obtained 2020.01.02
% filter by identity >= 30% and coverage >= 80%
idx_filtered = num(:,1) >= 30 & num(:,2) >= 80;
num_new = num(idx_filtered,:);
txt_new = txt(idx_filtered,:);
clear idx_filtered;
% choose the highest bit score with lowest evalue
id_list = txt_new(:,1);
pdb_list = txt_new(:,2);
evalue_list = num_new(:,3);
bit_list = num_new(:,4);

%% Process data

% 1. filter gene id
% split multiple genes, remove 'NA' and non-yeast genes
pdb_processing_1 = struct();
pdb_processing_1.pdbid = cell(0,1);
pdb_processing_1.geneid = cell(0,1);
pdb_processing_1.name = cell(0,1);
pdb_processing_1.cofactor = cell(0,1);
pdb_processing_1.peptide_entityid = zeros(0,1);
pdb_processing_1.prot_stoic = zeros(0,1); % number of peptides per complex
pdb_processing_1.cofactor_stoic = cell(0,1); % number of cofactors in the complex
pdb_processing_1.pdb_chains = cell(0,1);
pdb_processing_1.all_pdb_chains = cell(0,1);

for i = 1:length(pdb.pdbid)
    disp(['Stage 1: ' num2str(i) '/' num2str(length(pdb.pdbid))]);
    gene_tmp = pdb.geneid{i};
    gene_list_tmp = strsplit(gene_tmp,'; ')';
    gene_list_tmp = unique(gene_list_tmp);
    for j = 1:length(gene_list_tmp)
        gene_tmp_tmp = gene_list_tmp{j};
        if ~strcmp('NA',gene_tmp_tmp) && (strcmp('Q',gene_tmp_tmp(1))||strcmp('R',gene_tmp_tmp(1))||strcmp('Y',gene_tmp_tmp(1)))
            pdb_processing_1.pdbid = [pdb_processing_1.pdbid;pdb.pdbid(i)];
            pdb_processing_1.geneid = [pdb_processing_1.geneid;{gene_tmp_tmp}];
            pdb_processing_1.name = [pdb_processing_1.name;pdb.name(i)];
            pdb_processing_1.cofactor = [pdb_processing_1.cofactor;pdb.cofactor(i)];
            pdb_processing_1.peptide_entityid = [pdb_processing_1.peptide_entityid;pdb.peptide_entityid(i)];
            pdb_processing_1.prot_stoic = [pdb_processing_1.prot_stoic;pdb.prot_stoic(i)];
            pdb_processing_1.cofactor_stoic = [pdb_processing_1.cofactor_stoic;pdb.cofactor_stoic(i)];
            pdb_processing_1.pdb_chains = [pdb_processing_1.pdb_chains;pdb.pdb_chains(i)];
            pdb_processing_1.all_pdb_chains = [pdb_processing_1.all_pdb_chains;pdb.all_pdb_chains(i)];
        end
    end
end

% 2. select gene-pdb pairs with highest blastp scores or highest coverage

genes_collected = unique(pdb_processing_1.geneid);

% % find the lowest identity and coverage
% str_collected = cell(0,1);
% for i = 1:length(genes_collected)
%     disp(['Define cutoff of blastp: ' num2str(i) '/' num2str(length(genes_collected))]);
%     gene_tmp = genes_collected(i);
%     idx_tmp = ismember(pdb_processing_1.geneid,gene_tmp);
%     pdbid_tmp = pdb_processing_1.pdbid(idx_tmp);
%     pdbchain_tmp = pdb_processing_1.all_pdb_chains(idx_tmp);
%     for k = 1:length(pdbid_tmp)
%         pdbid_tmp_tmp = pdbid_tmp(k);
%         pdbchain_tmp_tmp = split(pdbchain_tmp(k),'; ')';
%         for j = 1:length(pdbchain_tmp_tmp)
%             str_tmp = strcat(pdbid_tmp_tmp,'_',pdbchain_tmp_tmp(j),':',gene_tmp);
%             str_collected = [str_collected;str_tmp];
%         end
%     end
% end
% str_all_pdb = strcat(txt(:,2),':',txt(:,1));
% 
% identity = num(:,1);
% coverage = num(:,2);
% evalue = num(:,3);
% blastp_pdb = struct;
% blastp_pdb.str_collected_list = str_collected;
% blastp_pdb.identity_list = zeros(length(str_collected),1);
% blastp_pdb.coverage_list = zeros(length(str_collected),1);
% blastp_pdb.evalue_list = zeros(length(str_collected),1);
% for i = 1:length(str_collected)
%     disp(['Define cutoff of blastp: ' num2str(i) '/' num2str(length(str_collected))]);
%     str_tmp = str_collected(i);
%     if ismember(str_tmp,str_all_pdb)
%         idx_tmp = ismember(str_all_pdb,str_tmp);
%         identity_tmp = identity(idx_tmp);
%         coverage_tmp = coverage(idx_tmp);
%         evalue_tmp = evalue(idx_tmp);
%         idx_tmp_tmp = find(evalue_tmp == min(evalue_tmp),1);
%         blastp_pdb.identity_list(i) = identity_tmp(idx_tmp_tmp);
%         blastp_pdb.coverage_list(i) = coverage_tmp(idx_tmp_tmp);
%         blastp_pdb.evalue_list(i) = evalue_tmp(idx_tmp_tmp);
%     else
%         blastp_pdb.identity_list(i) = NaN;
%         blastp_pdb.coverage_list(i) = NaN;
%         blastp_pdb.evalue_list(i) = NaN;
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%
% min(blastp_pdb.identity_list);
% min(blastp_pdb.coverage_list);
% max(blastp_pdb.evalue_list);
% quantile(blastp_pdb.identity_list,0.01,1);
% quantile(blastp_pdb.coverage_list,0.01,1);
% quantile(blastp_pdb.evalue_list,0.99,1);
% %%%%%%%%%%%%%%%%%%%%%%

pdb_processing_2 = struct();
pdb_processing_2.pdbid = cell(0,1);
pdb_processing_2.geneid = cell(0,1);
pdb_processing_2.name = cell(0,1);
pdb_processing_2.cofactor = cell(0,1);
pdb_processing_2.peptide_entityid = zeros(0,1);
pdb_processing_2.prot_stoic = zeros(0,1); % number of peptides per complex
pdb_processing_2.cofactor_stoic = cell(0,1); % number of cofactors in the complex
pdb_processing_2.pdb_chains = cell(0,1);
pdb_processing_2.all_pdb_chains = cell(0,1);
for i = 1:length(genes_collected)
    disp(['Stage 2: ' num2str(i) '/' num2str(length(genes_collected))]);
    gene_tmp = genes_collected(i);
end









