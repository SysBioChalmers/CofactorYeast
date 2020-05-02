%% setMedia
function model = setMedia(model,type)

% type  = 1: minimal media (Delft media) (default)
%       = 2: yeast nitrogen base without amino acids

exchangeRxns = findExcRxns(model);
model.lb(exchangeRxns) = 0;
model.ub(exchangeRxns) = 1000;

blockedExchanges = {'r_1663'; ... % bicarbonate exchange
                    'r_4062'; ... % lipid backbone exchange
                    'r_4064'};    % lipid chain exchange
% 1: minimal media (Delft media)                
desiredExchanges_1 = {'r_1654'; ... % ammonium exchange
                      'r_2005'; ... % phosphate exchange
                      'r_2060'; ... % sulphate exchange
                      'r_2020'; ... % potassium exchange
                      'r_1832'; ... % hydrogen exchange
                      'r_4597'; ... % Mg(2+) exchange
                      'r_2100'; ... % water exchange
                      'r_1861'; ... % iron(2+) exchange
                      'r_4596'; ... % Zn(2+) exchange
                      'r_4600'; ... % Ca(2+) exchange
                      'r_4593'; ... % chloride exchange
                      'r_4595'; ... % Mn(2+) exchange
                      'r_4594'; ... % Cu(2+) exchange
                      'r_2049'; ... % sodium exchange
                      'r_1671'; ... % biotin exchange
                      'r_2067'; ... % thiamine exchange
                      'r_2028'; ... % pyridoxine exchange >high uptake<
                      'r_1967'; ... % nicotinate exchange
                      'r_1947'; ... % myo-inositol exchange
                      'r_1604'; ... % 4-aminobenzoate exchange
                      'r_1548'};    % (R)-pantothenate exchange
% 2: yeast nitrogen base without amino acids
desiredExchanges_2 = {'r_1654'; ... % ammonium exchange
                      'r_2005'; ... % phosphate exchange
                      'r_2060'; ... % sulphate exchange
                      'r_2020'; ... % potassium exchange
                      'r_1832'; ... % hydrogen exchange
                      'r_4597'; ... % Mg(2+) exchange
                      'r_2100'; ... % water exchange
                      'r_1861'; ... % iron(2+) exchange
                      'r_4596'; ... % Zn(2+) exchange
                      'r_4600'; ... % Ca(2+) exchange
                      'r_4593'; ... % chloride exchange
                      'r_4595'; ... % Mn(2+) exchange
                      'r_4594'; ... % Cu(2+) exchange
                      'r_2049'; ... % sodium exchange
                      'r_1671'; ... % biotin exchange
                      'r_2067'; ... % thiamine exchange
                      'r_2028'; ... % pyridoxine exchange >high uptake<
                      'r_1967'; ... % nicotinate exchange
                      'r_1947'; ... % myo-inositol exchange
                      'r_1604'; ... % 4-aminobenzoate exchange
                      'r_1548'; ... % (R)-pantothenate exchange
                      'r_1792'; ... % folic acid exchange
                      'r_2038'};    % riboflavin exchange

if type == 1
    desiredExchanges = desiredExchanges_1;
elseif type == 2
    desiredExchanges = desiredExchanges_2;
    %remove iron(2+)
    desiredExchanges = desiredExchanges(~ismember(desiredExchanges,'r_1861'));
    %add iron(3+)
    desiredExchanges = [desiredExchanges;model.rxns(ismember(model.rxnNames,'iron(3+) exchange'))];
end

uptakeRxnIndexes     = findRxnIDs(model,desiredExchanges);
blockedRxnIndex      = findRxnIDs(model,blockedExchanges);

if length(find(uptakeRxnIndexes~= 0)) ~= length(desiredExchanges)
    warning('Not all exchange reactions were found.');
end

model.lb(uptakeRxnIndexes(uptakeRxnIndexes~=0)) = -1000;

model.lb(blockedRxnIndex) = 0;
model.ub(blockedRxnIndex) = 0;

