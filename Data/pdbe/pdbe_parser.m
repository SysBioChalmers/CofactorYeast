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
pdb.pdb_chains = cell(0,1);
pdb.all_pdb_chains = cell(0,1);

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
    try
        data_ligand = webread(strcat(retrieve_ligand,pdbid));
    catch
        data_ligand = [];
    end
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
    
    % extract assembly information
    tmp_assembly = getfield(data_assembly,name);
    tmp_assembly = tmp_assembly(strcmp({tmp_assembly.assembly_id},assembly_id)).entities;
    if length(tmp_assembly) == 1
        if strcmp(tmp_assembly.molecule_type,'polypeptide(L)')
            check = 1;
            peptd_name = tmp_assembly.molecule_name;
            peptd_entityid = tmp_assembly.entity_id;
            peptd_copies = tmp_assembly.number_of_copies;
            peptd_pdb_chain = tmp_assembly.in_chains;
        else
            check = 0;
        end
    else
        tmp_assembly = struct2table(tmp_assembly);
        if ismember({'polypeptide(L)'},tmp_assembly.molecule_type)
            check = 1;
            idx_peptd = ismember(tmp_assembly.molecule_type,{'polypeptide(L)'});
            peptd_name = tmp_assembly.molecule_name(idx_peptd);
            peptd_entityid = tmp_assembly.entity_id(idx_peptd);
            peptd_copies = tmp_assembly.number_of_copies(idx_peptd);
            peptd_pdb_chain = tmp_assembly.in_chains(idx_peptd);
        else
            check = 0;
        end
    end
    
    if check == 1
        pdb.pdbid = [pdb.pdbid;repelem({pdbid},length(peptd_entityid))'];
        pdb.name = [pdb.name;peptd_name];
        pdb.peptide_entityid = [pdb.peptide_entityid;peptd_entityid];
        pdb.prot_stoic = [pdb.prot_stoic;peptd_copies]; % number of peptides per complex

        % map pdb id with gene id
        % collect chains
        chainlist_pdbgene = cell(0,1);
        tmp_molecules = getfield(data_molecules,name);
        % order molecule data
        if length(tmp_molecules) == 1
            entityname = strcat('entity_',mat2str(tmp_molecules.entity_id),' = tmp_molecules;');
            eval(entityname);
        else
            if iscell(tmp_molecules)
                for j = 1:length(tmp_molecules)
                    tmp = tmp_molecules{j};
                    entityname = strcat('entity_',mat2str(tmp.entity_id),' = tmp;');
                    eval(entityname);
                end
            elseif isstruct(tmp_molecules)
                for j = 1:length(tmp_molecules)
                    tmp = tmp_molecules(j);
                    entityname = strcat('entity_',mat2str(tmp.entity_id),' = tmp;');
                    eval(entityname);
                end
            end
        end

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
        if isempty(data_ligand)
            chainlist_chainasym = cell(0,1);
            asymlist_chainasym = cell(0,1);
            for j = 1:length(peptd_entityid)
                tmp = peptd_entityid(j);
                read_struct = strcat('data_tmp = entity_',mat2str(peptd_entityid(j)),';');
                eval(read_struct);
                chainlist_chainasym = [chainlist_chainasym;data_tmp.in_chains];
                asymlist_chainasym = [asymlist_chainasym;data_tmp.in_struct_asyms];
            end
        else
            tmp_ligand = getfield(data_ligand,name);
            if isfield(tmp_ligand,'author_insertion_code')
                tmp_ligand = rmfield(tmp_ligand,'author_insertion_code');
            end
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
        end

        if ismember('Bound',tmp_assembly.molecule_type)
            idx_bound = ismember(tmp_assembly.molecule_type,{'Bound'});
            bound_name = tmp_assembly.molecule_name(idx_bound);
            bound_copies = tmp_assembly.number_of_copies(idx_bound);
            bound_entity = tmp_assembly.entity_id(idx_bound);
        else
            bound_name = [];
            bound_copies = [];
            bound_entity = [];
        end

        % add gene id, chains and cofactors
        for j = 1:length(peptd_entityid)
            tmp_asymchain = peptd_pdb_chain{j};
            maxstr = max(cell2mat(cellfun(@(x) length(x),asymlist_chainasym,'UniformOutput',false)));
            if ~ischar(tmp_asymchain)
                maxstr_tmp = max(cell2mat(cellfun(@(x) length(x),tmp_asymchain,'UniformOutput',false)));
            else
                maxstr_tmp = length(tmp_asymchain);
            end
            if maxstr_tmp > maxstr
                if maxstr_tmp == 2
                    if ~ischar(tmp_asymchain)
                        tmp_asymchain = cellfun(@(x) x(1),tmp_asymchain,'UniformOutput',false);
                    else
                        tmp_asymchain = tmp_asymchain(1);
                    end
                elseif maxstr_tmp == 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
            [~, b] = ismember(tmp_asymchain,asymlist_chainasym);
            tmp_pdb_chain = chainlist_chainasym(b);
            pdb.pdb_chains = [pdb.pdb_chains;{strjoin(tmp_pdb_chain,'; ')}];

            [~, b] = ismember(tmp_pdb_chain,chainlist_pdbgene);
            tmp_gene = unique(geneidlist_pdbgene(b));
            pdb.geneid = [pdb.geneid;{strjoin(tmp_gene,'; ')}];

            read_struct = strcat('data_tmp = entity_',mat2str(peptd_entityid(j)),';');
            eval(read_struct);
            tmp_all_pdb_chains = data_tmp.in_chains;
            pdb.all_pdb_chains = [pdb.all_pdb_chains;{strjoin(tmp_all_pdb_chains,'; ')}];

            % add cofactors
            if ~isempty(bound_name)
                bounds_cell = cell(0,1);
                bounds_copies = zeros(0,1);
                for k = 1:length(bound_name)
                    tmp_b = bound_name(k);
                    [~, b] = ismember(tmp_b,tmp_assembly.molecule_name);
                    tmp_tmp = tmp_assembly.in_chains{b};
                    if ~ischar(tmp_tmp)
                        maxstr_tmp_tmp = max(cell2mat(cellfun(@(x) length(x),tmp_tmp,'UniformOutput',false)));
                    else
                        maxstr_tmp_tmp = length(tmp_tmp);
                    end
                    if maxstr_tmp_tmp > maxstr
                        if maxstr_tmp_tmp == 2
                            if ~ischar(tmp_tmp)
                                tmp_tmp = cellfun(@(x) x(1),tmp_tmp,'UniformOutput',false);
                            else
                                tmp_tmp = tmp_tmp(1);
                            end
                        elseif maxstr_tmp_tmp == 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end
                    end
                    tmp_c = 0;
                    for m = 1:length(tmp_tmp)
                        tmp_tmp_tmp = chainlist_chainasym(ismember(asymlist_chainasym,tmp_tmp(m)));
                        if ismember(tmp_tmp_tmp,tmp_pdb_chain)
                            tmp_c = tmp_c + 1;
                        end
                    end
                    if tmp_c > 0
                        bounds_cell = [bounds_cell;tmp_b];
                        bounds_copies = [bounds_copies;tmp_c];
                    end
                end
                if ~isempty(bounds_cell)
                    tmp_bound_copies = strjoin(cellstr(num2str(bounds_copies)),'; ');
                    pdb.cofactor = [pdb.cofactor;{strjoin(bounds_cell,'; ')}];
                    pdb.cofactor_stoic = [pdb.cofactor_stoic;tmp_bound_copies];
                else
                    pdb.cofactor = [pdb.cofactor;{'NA'}];
                    pdb.cofactor_stoic = [pdb.cofactor_stoic;{'NA'}];
                end
            else
                pdb.cofactor = [pdb.cofactor;{'NA'}];
                pdb.cofactor_stoic = [pdb.cofactor_stoic;{'NA'}];
            end
        end
    end
end

% pdbid = '6oaa'
% pdbid = '3fwc'
% pdbid = '6cp7'
% pdbid = '1a6r'
% pdbid = '6b8h'
% pdbid = '1ld4'
% pdbid = '1fnt'
% webread('https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/117e');





