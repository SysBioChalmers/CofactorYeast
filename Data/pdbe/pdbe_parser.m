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
pdb.peptide_entityid = cell(0,1);
pdb.prot_stoic = zeros(0,1); % number of peptides per complex
pdb.cofactor_stoic = zeros(0,1); % number of cofactors per peptide

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
    tmp_assembly = struct2cell(getfield(data_assembly,name));
    idx_tmp = cell2mat(tmp_assembly(6,:)) == assembly_id;
    tmp_assembly = struct2table(tmp_assembly{4,idx_tmp});
    
    idx_bound = ismember(tmp_assembly.molecule_type,'Bound');
    bound_name = tmp_assembly.molecule_name(idx_bound);
    bound_copies = tmp_assembly.number_of_copies(idx_bound);
    bound_entityid = tmp_assembly.entity_id(idx_bound);
    
    pdb_summary_cofactor.pdbid = [pdb_summary_cofactor.pdbid;repelem({pdbid},length(bound_name))'];
    pdb_summary_cofactor.cofactor = [pdb_summary_cofactor.cofactor;bound_name];
    pdb_summary_cofactor.cofactor_stoic = [pdb_summary_cofactor.cofactor_stoic;bound_copies];
    
    idx_peptd = ismember(tmp_assembly.molecule_type,'polypeptide(L)');
    peptd_name = tmp_assembly.molecule_name(idx_peptd);
    peptd_entityid = tmp_assembly.entity_id(idx_peptd);

    
%     1a3x
%     6bnp
%     ans1 = webread('https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/1a3x');
%     ans2 = webread('https://www.ebi.ac.uk/pdbe/api/mappings/ensembl/1a3x');
%     ans3 = webread('https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/1a3x');







end









