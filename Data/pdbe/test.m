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

retrieve_summary = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/';
retrieve_molecules = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/';
retrieve_assembly = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/';

maxstring_molecule = 0;
id_molecule = cell(0,1);

for i = 1:length(pdbeids)
    disp([num2str(i) '/' num2str(length(pdbeids))]);
    pdbid = pdbeids{i};
    data_summary = webread(strcat(retrieve_summary,pdbid));
    data_molecules = webread(strcat(retrieve_molecules,pdbid));
    data_assembly = webread(strcat(retrieve_assembly,pdbid));
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
            peptd_pdb_chain = tmp_assembly.in_chains;
        else
            check = 0;
        end
    else
        tmp_assembly = struct2table(tmp_assembly);
        if ismember({'polypeptide(L)'},tmp_assembly.molecule_type)
            check = 1;
            idx_peptd = ismember(tmp_assembly.molecule_type,{'polypeptide(L)'});
            peptd_pdb_chain_tmp = tmp_assembly.in_chains(idx_peptd);
            peptd_pdb_chain = cell(0,1);
            for j = 1:length(peptd_pdb_chain_tmp)
                peptd_pdb_chain = [peptd_pdb_chain;peptd_pdb_chain_tmp{j}];
            end
        else
            check = 0;
        end
    end
    
    
    
    if check == 1
        if ~ischar(peptd_pdb_chain)
            maxstr_assembly = max(cell2mat(cellfun(@(x) length(x),peptd_pdb_chain,'UniformOutput',false)));
        else
            maxstr_assembly = length(peptd_pdb_chain);
        end
        % collect chains
        tmp_molecules = getfield(data_molecules,name);
        % order molecule data
        if length(tmp_molecules) == 1
            if strcmp(tmp_molecules.molecule_type,'polypeptide(L)')
                a = max(cell2mat(cellfun(@(x) length(x),tmp_molecules.in_struct_asyms,'UniformOutput',false)));
                if a > maxstring_molecule
                    maxstring_molecule = a;
                end
                if a > 1 && maxstr_assembly > a
                    id_molecule = [id_molecule;{pdbid}];
                end
            end
        else
            if iscell(tmp_molecules)
                for j = 1:length(tmp_molecules)
                    tmp = tmp_molecules{j};
                    if strcmp(tmp.molecule_type,'polypeptide(L)')
                        a = max(cell2mat(cellfun(@(x) length(x),tmp.in_struct_asyms,'UniformOutput',false)));
                        if a > maxstring_molecule
                            maxstring_molecule = a;
                        end
                        if a > 1 && maxstr_assembly > a
                            id_molecule = [id_molecule;{pdbid}];
                        end
                    end
                end
            elseif isstruct(tmp_molecules)
                for j = 1:length(tmp_molecules)
                    tmp = tmp_molecules(j);
                    if strcmp(tmp.molecule_type,'polypeptide(L)')
                        a = max(cell2mat(cellfun(@(x) length(x),tmp.in_struct_asyms,'UniformOutput',false)));
                        if a > maxstring_molecule
                            maxstring_molecule = a;
                        end
                        if a > 1 && maxstr_assembly > a
                            id_molecule = [id_molecule;{pdbid}];
                        end
                    end
                end
            end
        end
    end
end

id_molecule = unique(id_molecule);


% pdbid = '1fnt';
