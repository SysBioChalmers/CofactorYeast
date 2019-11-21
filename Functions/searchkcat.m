%% searchkcat
function [finalkcat, conf_score] = searchkcat(ec,substrate,org_name,allkcats)

% Input:  ec, EC number, e.g., 'EC1.1.1.1' or 'EC1.1.2.4; EC1.1.99.40';
%         substrate, metabolite name, e.g., 'acetyl-CoA; oxaloacetate';
%         org_name, organism name, e.g., 'saccharomyces cerevisiae';
%         allkcats, kcats from BRENDA database.
% Output: finalkcat, kcat value in the unit of "/h", max for multiple EC numbers;
%         conf_score, confidence score of the searched kcat value.
% CONFIDENCE SCORE
% 0:    No kcat reported for the EC number in BRENDA, then kcat = nan.
% 1:    Neither substrates nor organism matched, then kcat = median of all reported kcats.
% 2:    Substrates matched but organism not matched, then kcat = median of all substrates matched kcats.
% 3:    Both substrates and organism matched, then kcat = max of matched kcats.

% Note: if multiple EC numbers provided, then the kcat with highest
%       confidence score will be selected.

whole_ec = allkcats{1};
whole_substrate = allkcats{2};
whole_org = allkcats{3};
whole_kcat = allkcats{4};


eclist = split(ec,'; ');
kcatlist = zeros(length(eclist),1);
conflist = zeros(length(eclist),1);
substrate = lower(substrate);
sublist = split(substrate,'; ');
sublist_plus = cellfun(@(x) strcat(x,'+'),sublist,'UniformOutput',false);%nad/nadp and nad+/nadp+
sublist = [sublist;sublist_plus];

for i = 1:length(eclist)
    ec_tmp = eclist(i);
    if ~ismember(ec_tmp,whole_ec)
        kcat = nan;
        conf = 0;
    else
        idx_tmp = ismember(whole_ec,ec_tmp);
        sub_tmp = whole_substrate(idx_tmp);
        org_tmp = whole_org(idx_tmp);
        kcat_tmp = whole_kcat(idx_tmp);
        
        match_idx_sub = ismember(sub_tmp,sublist);
        match_idx_org = ismember(org_tmp,org_name);
        match_idx_combined = match_idx_sub & match_idx_org;
        
        if any(match_idx_combined) %both organism and substrate matched
            kcat = max(kcat_tmp(match_idx_combined)); %choose max
            conf = 3;
        else
            if any(match_idx_sub) %only substrate matched
                kcat = median(kcat_tmp(match_idx_sub)); %choose median
                conf = 2;
            else %if no substrate matched, then choose median of all kcats
                kcat = median(kcat_tmp);
                conf = 1;
            end
        end
    end
    kcatlist(i,1) = kcat;
    conflist(i,1) = conf;
end

finalkcat_tmp = kcatlist(conflist == max(conflist));

finalkcat = max(finalkcat_tmp);
conf_score = max(conflist);
end