%% pdbe_parser
%  Extract information of protein stoichiometry and cofactors from the PDBe
%  database.

%% Main code

% Download a list of pdb ids of Saccharomyces cerevisiae from PDBe, and
% save as "pdbe_yeast_ids.xlsx".

[~,txt,~] = xlsread('pdbe_yeast_ids.xlsx');

txt = txt(2:end,1);
txt = cellfun(@(x) x(1:4),txt,'UniformOutput',false);
pdbeids = unique(txt);
clear txt;

pdb = struct;
pdb.pdbid = cell(0,1);
pdb.geneid = cell(0,1);
pdb.name = cell(0,1);
pdb.cofactor = cell(0,1);
pdb.peptide_entityid = zeros(0,1);
pdb.prot_stoic = zeros(0,1); % number of peptides per complex
pdb.cofactor_stoic = cell(0,1); % number of cofactors per peptide
pdb.pdb_chain = cell(0,1);
pdb.all_pdb_chains = cell(0,1);

pdb_summary_prot_stoic = struct;
pdb_summary_prot_stoic.pdbid = cell(0,1);
pdb_summary_prot_stoic.form = cell(0,1);

pdb_summary_cofactor = struct;
pdb_summary_cofactor.pdbid = cell(0,1);
pdb_summary_cofactor.cofactor = cell(0,1);
pdb_summary_cofactor.cofactor_stoic = zeros(0,1); % number of cofactors per complex/assembly

retrieve_summary = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/';
retrieve_assembly = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/';
retrieve_ensembl = 'https://www.ebi.ac.uk/pdbe/api/mappings/ensembl/';
retrieve_molecules = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/';
retrieve_ligand = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/';


for i = 1:length(pdbeids)
    disp([num2str(i) '/' num2str(length(pdbeids))]);
    pdbid = pdbeids{i};
    data_summary = webread(strcat(retrieve_summary,pdbid));
    data_assembly = webread(strcat(retrieve_assembly,pdbid));
    try
        data_ensembl = webread(strcat(retrieve_ensembl,pdbid));
    catch
        data_ensembl = [];
    end
    data_molecules = webread(strcat(retrieve_molecules,pdbid));
    data_ligand = webread(strcat(retrieve_ligand,pdbid));
    
    name = strcat('x',pdbid);
    
    % check if there are multiple assemblies
    tmp_summary = getfield(data_summary,name);
	tmp = tmp_summary.assemblies;
	tmp = struct2table(tmp);
    if length(tmp_summary.assemblies) > 1
        assembly_id = cell2mat(tmp.assembly_id(tmp.preferred == true));
        form = strcat(tmp.form{tmp.preferred == true},tmp.name{tmp.preferred == true});
    else
        assembly_id = tmp.assembly_id; 
        form = strcat(tmp.form,tmp.name);
    end
    
    pdb_summary_prot_stoic.pdbid = [pdb_summary_prot_stoic.pdbid;{pdbid}];
    pdb_summary_prot_stoic.form = [pdb_summary_prot_stoic.form;{form}];
    
    % extract assembly information
    tmp_assembly = getfield(data_assembly,name);
    tmp_assembly = tmp_assembly(strcmp({tmp_assembly.assembly_id},assembly_id)).entities;
    tmp_assembly = struct2table(tmp_assembly);
    
    idx_bound = ismember(tmp_assembly.molecule_type,'Bound');
    bound_name = tmp_assembly.molecule_name(idx_bound);
    bound_copies = tmp_assembly.number_of_copies(idx_bound);
    bound_in_chains = tmp_assembly.in_chains(idx_bound);
    
    pdb_summary_cofactor.pdbid = [pdb_summary_cofactor.pdbid;repelem({pdbid},length(bound_name))'];
    pdb_summary_cofactor.cofactor = [pdb_summary_cofactor.cofactor;bound_name];
    pdb_summary_cofactor.cofactor_stoic = [pdb_summary_cofactor.cofactor_stoic;bound_copies];
    
    idx_peptd = ismember(tmp_assembly.molecule_type,'polypeptide(L)');
    peptd_name = tmp_assembly.molecule_name(idx_peptd);
    peptd_entityid = tmp_assembly.entity_id(idx_peptd);
    peptd_copies = tmp_assembly.number_of_copies(idx_peptd);
    peptd_pdb_chain = tmp_assembly.in_chains(idx_peptd);
    
    pdb.pdbid = [pdb.pdbid;repelem({pdbid},length(peptd_entityid))'];
    pdb.name = [pdb.name;peptd_name];
    pdb.peptide_entityid = [pdb.peptide_entityid;peptd_entityid];
    pdb.prot_stoic = [pdb.prot_stoic;peptd_copies]; % number of peptides per complex
    
    % map pdb id with gene id
    % collect chains
    chainlist_pdbgene = cell(0,1);
    tmp_molecules = getfield(data_molecules,name);
    % order molecule data
    for j = 1:length(tmp_molecules)
        tmp = tmp_molecules{j};
        entityname = strcat('entity_',mat2str(tmp.entity_id),' = tmp;');
        eval(entityname);
    end
%     tmp_molecules_peptd = tmp_molecules(idx_peptd);
%     for j = 1:length(tmp_molecules_peptd)
%         tmp = tmp_molecules_peptd{j};
%         chainlist_pdbgene = [chainlist_pdbgene;tmp.in_chains];
%     end

    for j = 1:length(peptd_entityid)
        read_struct = strcat('data_tmp = entity_',mat2str(peptd_entityid(j)),';');
        eval(read_struct);
        chainlist_pdbgene = [chainlist_pdbgene;data_tmp.in_chains];
    end

    if isempty(data_ensembl)
        geneidlist_pdbgene = repelem({'NA'},length(chainlist_pdbgene))';
    else
        tmp_ensembl = getfield(data_ensembl,name);
        tmp_ensembl = tmp_ensembl.Ensembl;
        tmp_ensembl = struct2cell(tmp_ensembl);
        
        tmp_relationship_gene = cell(0,1);
        tmp_relationship_chain = cell(0,1);
        for j = 1:length(tmp_ensembl)
            tmp = tmp_ensembl{j};
            tmp_table = tmp.mappings;
            if isfield(tmp_table,'end')
                tmp_table = rmfield(tmp_table,'end');
            end
            tmp_table = struct2table(tmp_table);
            tmp_relationship_gene = [tmp_relationship_gene;tmp_table.translation_id];
            tmp_relationship_chain = [tmp_relationship_chain;tmp_table.chain_id];
        end
        geneidlist_pdbgene = cell(length(chainlist_pdbgene),1);
        for j = 1:length(chainlist_pdbgene)
            tmp = chainlist_pdbgene(j);
            if ismember(tmp,tmp_relationship_chain)
                tmp_gene = tmp_relationship_gene(ismember(tmp_relationship_chain,tmp));
                tmp_gene = unique(tmp_gene);
                tmp_gene = {strjoin(tmp_gene,'; ')};
            else
                tmp_gene = {'NA'};
            end
            geneidlist_pdbgene(j) = tmp_gene;
        end
    end
    
    % map chain id with struct asym id
    tmp_ligand = getfield(data_ligand,name);
    tmp_ligand = struct2table(tmp_ligand);
    chainlist_chainasym = tmp_ligand.chain_id;
    asymlist_chainasym = tmp_ligand.struct_asym_id;
    for j = 1:length(peptd_entityid)
        tmp = peptd_entityid(j);
        read_struct = strcat('data_tmp = entity_',mat2str(peptd_entityid(j)),';');
        eval(read_struct);
        chainlist_chainasym = [chainlist_chainasym;data_tmp.in_chains];
        asymlist_chainasym = [asymlist_chainasym;data_tmp.in_struct_asyms];
    end
    
    % add gene id and chains
    for j = 1:length(peptd_entityid)
        tmp_asymchain = peptd_pdb_chain{j};
        [~, b] = ismember(tmp_asymchain,asymlist_chainasym);
        tmp_pdb_chain = chainlist_chainasym(b);
        pdb.pdb_chain = [pdb.pdb_chain;{strjoin(tmp_pdb_chain,'; ')}];
        
        [~, b] = ismember(tmp_pdb_chain,chainlist_pdbgene);
        tmp_gene = unique(geneidlist_pdbgene(b));
        pdb.geneid = [pdb.geneid;{strjoin(tmp_gene,'; ')}];
        
        read_struct = strcat('data_tmp = entity_',mat2str(peptd_entityid(j)),';');
        eval(read_struct);
        tmp_all_pdb_chains = data_tmp.in_chains;
        pdb.all_pdb_chains = [pdb.all_pdb_chains;{strjoin(tmp_all_pdb_chains,'; ')}];
        
        % add cofactors!!!
%         pdb.cofactor = cell(0,1);
%         pdb.cofactor_stoic = cell(0,1); % number of cofactors per peptide
    end
    
    

end

pdbid = '6oaa'
pdbid = '3fwc'
webread('https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/117e');





