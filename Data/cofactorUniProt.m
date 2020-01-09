% cofactorUniProt

[~,txt,~] = xlsread('cofactor_UniProt.xlsx','Rawdata_20200107');
proteinlist = txt(2:end,1);
cofactorlist = txt(2:end,6);
clear txt;

CofactorUniProt = struct;
CofactorUniProt.cofactor = cell(0,1);
CofactorUniProt.protein = cell(0,1);
CofactorUniProt.copy = zeros(0,1);
CofactorUniProt.source = cell(0,1);

for i = 1:length(proteinlist)
    proteinid = proteinlist(i);
    rawcofactor = cofactorlist(i);
    tmp = split(rawcofactor,'COFACTOR: ');
    for j = 1:length(tmp)
        if contains(tmp(j),'Name=')
            tmptmp = split(tmp(j),';');
            nametmp = tmptmp{contains(tmptmp,'Name=')};
            nametmp = extractAfter(nametmp,'Name=');
            nametmp = strtrim(nametmp);
            if contains(tmp(j),'Note=Binds')
                sourcetmp = {'UniProt'};
                copytmp = tmptmp{contains(tmptmp,'Note=Binds')};
                copytmp = strtrim(copytmp);
                copytmp = extractAfter(copytmp,'Note=Binds');
                copytmp = strtrim(copytmp);
                copynumbertmp = extractBefore(copytmp,' ');
                if strcmp(copynumbertmp,'a')
                    copynumbertmp = 1;
                elseif strcmp(copynumbertmp,'two')
                    copynumbertmp = 2;
                else
                    copynumbertmp = str2double(copynumbertmp);
                end
                if contains(copytmp,'per dimer') || contains(copytmp,'per homodimer')
                    copynumbertmp = copynumbertmp/2;
                elseif contains(copytmp,'per heterotetramer')
                    copynumbertmp = copynumbertmp/4;
                end
            else
                copynumbertmp = 1;
                sourcetmp = {'UniProt, copy assumed to be 1'};
            end
            CofactorUniProt.protein = [CofactorUniProt.protein;proteinid];
            CofactorUniProt.cofactor = [CofactorUniProt.cofactor;{nametmp}];
            CofactorUniProt.copy = [CofactorUniProt.copy;copynumbertmp];
            CofactorUniProt.source = [CofactorUniProt.source;sourcetmp];
        end
    end
end
















