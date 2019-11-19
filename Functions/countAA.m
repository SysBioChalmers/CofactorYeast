%% countAA

function sum = countAA(seq,aa_list)

if iscell(seq)
    seq = cell2mat(seq);
end

sum = struct();
sum.subs = aa_list.subs;
sum.prod = aa_list.prod;
sum.abbr = aa_list.aa;
sum.num = [];

for i = 1:length(sum.abbr)
    AA_abbr = sum.abbr{i};
    sum.num(i,1) = length(strfind(seq,AA_abbr));
end