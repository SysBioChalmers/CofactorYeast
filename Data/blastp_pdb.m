%% BlastP-PDB database
% To determine stoichiometry information of proteins.

% Timing: ~ 160 s

tic;
% "pdb_seqres.txt"
% PDB sequences should be first downloaded from the website:
% (https://www.rcsb.org/pages/general/summaries).

% "orf_trans_all_R64-2-1_20150113.fasta"
% Yeast sequences were downloaded (2019.11.18) from SGD website:
% (https://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/).

% Use the following command to perform blastp [Would take several hours]
% blastp -subject pdb_seqres.txt -query orf_trans_all_R64-2-1_20150113.fasta -out results.out -outfmt "7 qseqid sseqid pident qcovs evalue bitscore"
% The generated file "results.out" was converted to the excel file
% "blastp_res.xlsx" and moved to the current folder.

[num,txt,~] = xlsread('blastp_res.xlsx');

[~,pdb_raw,~] = xlsread('pdb_stoichiometry_20191118.xlsx');
pdb_id_list = pdb_raw(2:end,1);
sto_info_list = pdb_raw(2:end,2);

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

res = struct();
res.prot_id = cell(0,1);
res.pdb_id = cell(0,1);
res.peptide_n = zeros(0,1);
res.stoic_info = cell(0,1);

unq_id_list = unique(id_list);
for i = 1:length(unq_id_list)
    disp([num2str(i) '/' num2str(length(unq_id_list))]);
    id = unq_id_list(i);
    
        idx = ismember(id_list,id);
        pdb = pdb_list(idx);
        bit_score = bit_list(idx);
        evalue = evalue_list(idx);
        idx_tmp = bit_score == max(bit_score) & evalue == min(evalue);
        pdb_tmp1 = pdb(idx_tmp);
        pdb_tmp2 = cellfun(@(x) x(1:4),pdb_tmp1,'UniformOutput',false);
        pdb_tmp3 = unique(pdb_tmp2);
        
        pep_n_tmp = zeros(length(pdb_tmp3),1);
        stoic_tmp = cell(length(pdb_tmp3),1);
        for j = 1:length(pdb_tmp3)
            id_tmp = pdb_tmp3(j);
            pep_n_tmp(j,1) = sum(contains(pdb_tmp2,id_tmp));
            
            if sum(strcmpi(pdb_id_list,id_tmp)) > 0
                stoic_name = sto_info_list(strcmpi(pdb_id_list,id_tmp));
                stoic_tmp(j,1) = stoic_name;
            else
                stoic_tmp(j,1) = {'NA'};
            end
        end

        n = length(res.prot_id);
        res.pdb_id(n+1:n+length(pdb_tmp3),1) = pdb_tmp3;
        res.prot_id(n+1:n+length(pdb_tmp3),1) = repelem(id,length(pdb_tmp3));
        res.peptide_n(n+1:n+length(pdb_tmp3),1) = pep_n_tmp;
        res.stoic_info(n+1:n+length(pdb_tmp3),1) = stoic_tmp;

end
save('ProteinStoichiometry.mat','res');
toc;
clear;