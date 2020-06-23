%% Estimate mitochondrial proteome
load('CofactorYeast.mat');

[num1,gene1,~] = xlsread('Exp_proteomics_data.xlsx','PMID28365149');
[num2,gene2,~] = xlsread('Exp_proteomics_data.xlsx','PMID32213592');
[num3,gene3,~] = xlsread('Exp_proteomics_data.xlsx','PMID32312967');
gene1 = gene1(2:end,1);
gene2 = gene2(2:end,1);
gene3 = gene3(2:end,1);

[u_mw,u_gene,~] = xlsread('uniprot_protein_mw.xlsx');
u_gene = u_gene(2:end,2);

% convert to mol/gCDW
num1 = num1/6.02e23/1e-12;
num2 = num2/6.02e23/1e-12;
num3 = num3/1e15/1e-3;
[~,b1] = size(num1);
[~,b2] = size(num2);
[~,b3] = size(num3);

[mito_list, reactions] = collectCompartment(model,{'m','mm'});
mito_data = zeros(length(mito_list),b1+b2+b3);
for i = 1:length(mito_list)
    if ismember(mito_list(i),u_gene)
        mw = u_mw(ismember(u_gene,mito_list(i)));
    else
        mw = median(u_mw);
    end
    if ismember(mito_list(i),gene1)
        mito_data(i,1:b1) = num1(ismember(gene1,mito_list(i)),:);
    end
    if ismember(mito_list(i),gene2)
        mito_data(i,b1+1:b1+b2) = num2(ismember(gene2,mito_list(i)),:);
    end
    if ismember(mito_list(i),gene3)
        mito_data(i,b1+b2+1:b1+b2+b3) = num3(ismember(gene3,mito_list(i)),:);
    end
    mito_data(i,:) = mito_data(i,:) * mw;
end
mito_data(isnan(mito_data)) = 0;
sum_mito_data = sum(mito_data);
scatter(1:1:length(sum_mito_data),sum_mito_data);
max(sum_mito_data)

