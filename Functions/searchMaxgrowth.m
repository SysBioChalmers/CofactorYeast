%% searchMaxgrowth
function [mu,fluxes] = searchMaxgrowth(model,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,precision,soplexpath,factor_k)

if exist('factor_k', 'var')
    if isempty(factor_k)
        factor_k = 1;
    end
else
    factor_k = 1;
end

if exist('precision', 'var')
    if isempty(precision)
        precision = 0.001;
    end
else
    precision = 0.001;
end

mu_low = 0;
mu_high = 1;

while mu_high-mu_low > precision
    mu_mid = (mu_low+mu_high)/2;
    model_tmp = changeRxnBounds(model,'r_2111',mu_mid,'b');
    disp(['mu = ' num2str(mu_mid)]);
    fileName = writeLP(model_tmp,mu_mid,f,f_mito,osenseStr,rxnID,enzymedata,factor_k_withoutcofator,factor_k);
    command = sprintf([soplexpath,' -s0 -g5 -t300 -f1e-20 -o1e-20 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 --real:fpfeastol=1e-9 --real:fpopttol=1e-9 %s > %s.out %s'],fileName,fileName);
    system(command,'-echo');
    [~,sol_status,sol_full] = readSoplexResult('Simulation.lp.out',model_tmp);
    disp(['solution status: ' sol_status]);
    
    if strcmp(sol_status,'optimal')
        mu_low = mu_mid;
        flux_tmp = sol_full;
    else
        mu_high = mu_mid;
    end
end

if mu_low > precision
    mu = mu_low;
    fluxes = flux_tmp;
else
    mu = 0;
    fluxes = zeros(length(model.rxns),1);
end


