%% calculateMW 
function enzymedata = calculateMW(enzymedata,ProteinSequence)

enzymedata.enzyme_MW = zeros(length(enzymedata.enzyme),1);

for i = 1:length(enzymedata.enzyme)
    subunit_list = enzymedata.subunit(i,:);
    subunit_list = subunit_list(~ismember(subunit_list,''));
    MW_sum = 0;
    for j = 1:length(subunit_list) 
        subunit_tmp = subunit_list(j);
        subunit_seq_tmp = ProteinSequence.seq(ismember(ProteinSequence.id,subunit_tmp));
        subunit_MW_tmp = getProteinMW(cell2mat(subunit_seq_tmp));
        MW_sum = MW_sum + subunit_MW_tmp * enzymedata.subunit_stoichiometry(i,j);
    end
    enzymedata.enzyme_MW(i,1) = MW_sum;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MW = getProteinMW(seq)
% (From GECKO)
aa_codes = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N', ...
            'O','P','Q','R','S','T','U','V','W','X','Y','Z'};
aa_MWs   = [71.08 114.60 103.14 115.09 129.11 147.17 57.05 137.14 ...
            113.16 113.16 128.17 113.16 131.20 114.10 255.31 97.12 ...
            128.13 156.19 87.08 101.10 150.04 99.13 186.21 126.50 ...
            163.17 128.62];
MW = 18;
for i = 1:length(aa_codes)
    count = length(strfind(seq,aa_codes{i}));
    MW = MW + count*aa_MWs(i);
end
end