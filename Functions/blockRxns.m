%% blockRxns 
function model = blockRxns(model)
% Some reactions should be block to avoid weird flux distributions.

% block Gcy1, an alternative glycerol dissimilation pathway that is active 
% under microaerobic conditions (PMID: 22979944)
model = changeRxnBounds(model,'r_0487',0,'b');
model.ub(ismember(model.rxns,'r_0487_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0487_withoutcofactor')) = 0;

% block newly added isozyme
model = changeRxnBounds(model,'r_0438_5',0,'b'); % ferrocytochrome-c:oxygen oxidoreductase
model.ub(ismember(model.rxns,'r_0438_5_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0438_5_withoutcofactor')) = 0;

% block newly added isozyme for enolase as Eno1 and Eno2 show higher protein levels
model = changeRxnBounds(model,'r_0366_1_fwd',0,'b'); % enolase
model = changeRxnBounds(model,'r_0366_1_rvs',0,'b'); % enolase
model.ub(ismember(model.rxns,'r_0366_1_fwd_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0366_1_fwd_withoutcofactor')) = 0;
model.ub(ismember(model.rxns,'r_0366_1_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0366_1_rvs_withoutcofactor')) = 0;
model = changeRxnBounds(model,'r_0366_4_fwd',0,'b'); % enolase
model = changeRxnBounds(model,'r_0366_4_rvs',0,'b'); % enolase
model.ub(ismember(model.rxns,'r_0366_4_fwd_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0366_4_fwd_withoutcofactor')) = 0;
model.ub(ismember(model.rxns,'r_0366_4_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0366_4_rvs_withoutcofactor')) = 0;
model = changeRxnBounds(model,'r_0366_5_fwd',0,'b'); % enolase
model = changeRxnBounds(model,'r_0366_5_rvs',0,'b'); % enolase
model.ub(ismember(model.rxns,'r_0366_5_fwd_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0366_5_fwd_withoutcofactor')) = 0;
model.ub(ismember(model.rxns,'r_0366_5_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0366_5_rvs_withoutcofactor')) = 0;

% Fe(3+) secretion should be blocked to avoid unlimited ATP and NADH
model.ub(ismember(model.rxnNames,'iron(3+) exchange')) = 0;

model = changeRxnBounds(model,'r_0886_1',0,'b'); % iso-reaction of PFK
model.ub(ismember(model.rxns,'r_0886_1_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0886_1_withoutcofactor')) = 0;
model = changeRxnBounds(model,'r_4262_fwd',0,'b'); % citrate hydroxymutase
model.ub(ismember(model.rxns,'r_4262_fwd_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4262_fwd_withoutcofactor')) = 0;
model = changeRxnBounds(model,'r_4262_rvs',0,'b'); % citrate hydroxymutase
model.ub(ismember(model.rxns,'r_4262_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4262_rvs_withoutcofactor')) = 0;

% block some reactions that done in PMID: 28779005.
model = changeRxnBounds(model,'r_2045_rvs',0,'b'); % serine transport from [m] to [c]
model.ub(ismember(model.rxns,'r_2045_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_2045_rvs_withoutcofactor')) = 0;
model = changeRxnBounds(model,'r_0659_fwd',0,'b'); % isocitrate dehydrogenase (NADP)
model.ub(ismember(model.rxns,'r_0659_fwd_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0659_fwd_withoutcofactor')) = 0;
model = changeRxnBounds(model,'r_0659_rvs',0,'b'); % isocitrate dehydrogenase (NADP)
model.ub(ismember(model.rxns,'r_0659_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0659_rvs_withoutcofactor')) = 0;

% model = changeRxnBounds(model,'r_0725_fwd',0,'b'); % methenyltetrahydrofolate cyclohydrolase
model.ub(ismember(model.rxns,'r_0725_fwd_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0725_fwd_withoutcofactor')) = 0;
% model = changeRxnBounds(model,'r_0918',0,'b'); % phosphoserine transaminase
model.ub(ismember(model.rxns,'r_0918_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_0918_withoutcofactor')) = 0;

model = changeRxnBounds(model,'r_4235',0,'b'); % weird reaction from glc to g6p
model.ub(ismember(model.rxns,'r_4235_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4235_withoutcofactor')) = 0;

model = changeRxnBounds(model,'r_4216_rvs',0,'b'); % block the reaction to produce FMN without ATP
model.ub(ismember(model.rxns,'r_4216_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4216_rvs_withoutcofactor')) = 0;

model = changeRxnBounds(model,'r_4264_rvs',0,'b'); % succinate:NAD+ oxidoreductase
model.ub(ismember(model.rxns,'r_4264_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4264_rvs_withoutcofactor')) = 0;

% The following two reactions account for the transport of pyruvate from
% [c] to [m] without enzyme cost, should be blocked.
model = changeRxnBounds(model,'r_1137',0,'b');
model.ub(ismember(model.rxns,'r_1137_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_1137_withoutcofactor')) = 0;
model = changeRxnBounds(model,'r_1138',0,'b');
model.ub(ismember(model.rxns,'r_1138_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_1138_withoutcofactor')) = 0;

% Block backward reactions for PLP synthesis.
plpbwrxns = {'r_4211_1_rvs' 'r_4211_2_rvs' 'r_4212_1_rvs' 'r_4212_2_rvs'};
model = changeRxnBounds(model,plpbwrxns,zeros(1,length(plpbwrxns)),'b');
model.ub(ismember(model.rxns,'r_4211_1_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4211_1_rvs_withoutcofactor')) = 0;
model.ub(ismember(model.rxns,'r_4211_2_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4211_2_rvs_withoutcofactor')) = 0;
model.ub(ismember(model.rxns,'r_4212_1_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4212_1_rvs_withoutcofactor')) = 0;
model.ub(ismember(model.rxns,'r_4212_2_rvs_withoutcofactor')) = 0;
model.lb(ismember(model.rxns,'r_4212_2_rvs_withoutcofactor')) = 0;


