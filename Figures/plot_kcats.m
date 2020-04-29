%% Compare kcats
% Two purposes
% 1. Determine which of substrate-matched and organism-matched is more
%    essential.
% 2. Calculate correlation of the kcat matching method.

org_id = 'saccharomyces cerevisiae';
% org_id = 'homo sapiens';
% org_id = 'escherichia coli';
% org_id = 'bacillus subtilis';

% Load all kcats from BRENDA database
allkcats = loadkcats;

whole_ec = allkcats{1};
whole_substrate = allkcats{2};
whole_org = allkcats{3};
whole_kcat = allkcats{4};

% Check if there exist duplicates of substrates
tmp = strcat(whole_ec,whole_substrate,whole_org);
if length(unique(tmp)) == length(tmp)
    disp('No duplicate substrates');
end
clear tmp;

% S. cerevisiae dataset
sc_idx = ismember(whole_org,org_id);

sc_data = struct();
sc_data.ec = whole_ec(sc_idx);
sc_data.sub = whole_substrate(sc_idx);
sc_data.kcat = whole_kcat(sc_idx);

% Non S. cerevisiae dataset
nonsc_data = struct();
nonsc_data.ec = whole_ec(~sc_idx);
nonsc_data.sub = whole_substrate(~sc_idx);
nonsc_data.ecsub = strcat(nonsc_data.ec,nonsc_data.sub);
nonsc_data.kcat = whole_kcat(~sc_idx);

clear sc_idx allkcats;

%% Purpose 1

% Reduce the dataset
tmp_sc = strcat(sc_data.ec,sc_data.sub);
tmp_nonsc = strcat(nonsc_data.ec,nonsc_data.sub);
tmp_idx = ismember(tmp_sc,tmp_nonsc);

sc_data_reduced = struct();
sc_data_reduced.ec = sc_data.ec(tmp_idx);
sc_data_reduced.sub = sc_data.sub(tmp_idx);
sc_data_reduced.kcat = sc_data.kcat(tmp_idx);
clear tmp_sc tmp_nonsc tmp_idx;

kcat_1 = zeros(length(sc_data_reduced.ec),1); % substrate-matched kcats
kcat_2 = zeros(length(sc_data_reduced.ec),1); % organism-matched kcats
kcat_3 = zeros(length(sc_data_reduced.ec),1); % just median kcats

for i = 1:length(sc_data_reduced.ec)
    ec_tmp = sc_data_reduced.ec(i);
    sub_tmp = sc_data_reduced.sub(i);
    ecsub_tmp = strcat(ec_tmp,sub_tmp);
    kcat1_tmp = median(nonsc_data.kcat(ismember(nonsc_data.ecsub,ecsub_tmp)));
    kcat_1(i) = kcat1_tmp;
    sub_list_tmp = sc_data_reduced.sub(ismember(sc_data_reduced.ec,ec_tmp));
    kcat_list_tmp = sc_data_reduced.kcat(ismember(sc_data_reduced.ec,ec_tmp));
    
    if length(kcat_list_tmp) > 1
        kcat2_tmp = median(kcat_list_tmp(~ismember(sub_list_tmp,sub_tmp)));
    else
        kcat2_tmp = 0;
    end
    
    kcat_2(i) = kcat2_tmp;
    
    kcat3_tmp = median(nonsc_data.kcat(ismember(nonsc_data.ec,ec_tmp)));
    kcat_3(i) = kcat3_tmp;
    
end

idx = (kcat_2 ~= 0);
kcatsc = sc_data_reduced.kcat(idx)/3600;
kcat1 = kcat_1(idx)/3600;
kcat2 = kcat_2(idx)/3600;
kcat3 = kcat_3(idx)/3600;
N_samples = length(kcatsc);

kcatsc = log10(kcatsc);
kcat1 = log10(kcat1);
kcat2 = log10(kcat2);
kcat3 = log10(kcat3);

rmse1 = sqrt(sum((kcatsc - kcat1).^2)/N_samples);
rmse2 = sqrt(sum((kcatsc - kcat2).^2)/N_samples);
rmse3 = sqrt(sum((kcatsc - kcat3).^2)/N_samples);
[R1,P1] = corrcoef(kcatsc,kcat1);
[R2,P2] = corrcoef(kcatsc,kcat2);
[R3,P3] = corrcoef(kcatsc,kcat3);

figure('Name','1');
hold on;
box on;
line([-7 7],[-7 7],'Color','k','LineStyle','-.');
scatter(kcatsc,kcat1,'filled','MarkerFaceAlpha',0.5);
xlim([-7 7]);
ylim([-7 7]);
set(gca,'FontSize',9,'FontName','Helvetica');
ylabel('log10 (median kcat)','FontSize',12,'FontName','Helvetica');
xlabel('log10 (matched kcat)','FontSize',12,'FontName','Helvetica');
title('matched substrate','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(rmse1,2))];
str_r = ['R = ',num2str(round(R1(1,2),2))];
str_n = ['N = ',num2str(N_samples)];
text(-6.5,6,str_rmse,'FontSize',14,'FontName','Helvetica');
text(-6.5,4.5,str_r,'FontSize',14,'FontName','Helvetica');
text(1.5,-6,str_n,'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[0 500 230 230]);
set(gca,'position',[0.15 0.15 0.7 0.7]);

figure('Name','2');
hold on;
box on;
line([-7 7],[-7 7],'Color','k','LineStyle','-.');
scatter(kcatsc,kcat2,'filled','MarkerFaceAlpha',0.5);
xlim([-7 7]);
ylim([-7 7]);
set(gca,'FontSize',9,'FontName','Helvetica');
ylabel('log10 (median kcat)','FontSize',12,'FontName','Helvetica');
xlabel('log10 (matched kcat)','FontSize',12,'FontName','Helvetica');
title('matched organism','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(rmse2,2))];
str_r = ['R = ',num2str(round(R2(1,2),2))];
str_n = ['N = ',num2str(N_samples)];
text(-6.5,6,str_rmse,'FontSize',14,'FontName','Helvetica');
text(-6.5,4.5,str_r,'FontSize',14,'FontName','Helvetica');
text(1.5,-6,str_n,'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[250 500 230 230]);
set(gca,'position',[0.15 0.15 0.7 0.7]);

figure('Name','3');
hold on;
box on;
line([-7 7],[-7 7],'Color','k','LineStyle','-.');
scatter(kcatsc,kcat3,'filled','MarkerFaceAlpha',0.5);
xlim([-7 7]);
ylim([-7 7]);
set(gca,'FontSize',9,'FontName','Helvetica');
ylabel('log10 (median kcat)','FontSize',12,'FontName','Helvetica');
xlabel('log10 (matched kcat)','FontSize',12,'FontName','Helvetica');
title('No match, just median','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(rmse3,2))];
str_r = ['R = ',num2str(round(R3(1,2),2))];
str_n = ['N = ',num2str(N_samples)];
text(-6.5,6,str_rmse,'FontSize',14,'FontName','Helvetica');
text(-6.5,4.5,str_r,'FontSize',14,'FontName','Helvetica');
text(1.5,-6,str_n,'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[500 500 230 230]);
set(gca,'position',[0.15 0.15 0.7 0.7]);



%% Purpose 2
combinedstring = strcat(whole_ec,whole_substrate,whole_org);

% Compare two methods
% Method A:
% 1. if both organism and substrate matched, then choose the kcat
% 2. if the match in step 1 is not found, then choose median of the other
%    orgnisms with the matched substrate.
% 3. if the match in step 2 is not found, then choose median of the
%    organism with all substrates.
% 4. if the match in step 3 is not found, then choose median of the EC
%    number.
kcat_a = zeros(length(sc_data.ec),1); % plan A kcats

for i = 1:length(sc_data.ec)
    ec_tmp = sc_data.ec(i);
    sub_tmp = sc_data.sub(i);
    cmbnstr_tmp = strcat(ec_tmp,sub_tmp,{org_id});
    idx_tmp = ismember(combinedstring,cmbnstr_tmp);
    whole_ec_tmp = whole_ec(~idx_tmp);
    whole_substrate_tmp = whole_substrate(~idx_tmp);
    whole_org_tmp = whole_org(~idx_tmp);
    whole_kcat_tmp = whole_kcat(~idx_tmp);
    
    if ~ismember(ec_tmp,whole_ec_tmp)
        kcatvalue = 0;
    else
        sublist = whole_substrate_tmp(ismember(whole_ec_tmp,ec_tmp));
        orglist = whole_org_tmp(ismember(whole_ec_tmp,ec_tmp));
        kcatlist = whole_kcat_tmp(ismember(whole_ec_tmp,ec_tmp));

        match_idx_sub = ismember(sublist,sub_tmp);
        match_idx_org = ismember(orglist,org_id);
        match_idx_combined = match_idx_sub & match_idx_org;

        if any(match_idx_combined) % step 1
            kcatvalue = kcatlist(match_idx_combined);
        else
            if any(match_idx_sub) % step 2
                kcatvalue = median(kcatlist(match_idx_sub));
            else
                if any(match_idx_org) % step 3
                    kcatvalue = median(kcatlist(match_idx_org));
                else % step 4
                    kcatvalue = median(kcatlist);
                end
            end
        end
        kcat_a(i) = kcatvalue;
    end
end

% Method B:
% 1. if both organism and substrate matched, then choose the kcat
% 2. if the match in step 1 is not found, then choose median of the
%    organism with all substrates.
% 3. if the match in step 2 is not found, then choose median of the other
%    orgnisms with the matched substrate.
% 4. if the match in step 3 is not found, then choose median of the EC
%    number.

kcat_b = zeros(length(sc_data.ec),1); % plan B kcats

for i = 1:length(sc_data.ec)
    ec_tmp = sc_data.ec(i);
    sub_tmp = sc_data.sub(i);
    cmbnstr_tmp = strcat(ec_tmp,sub_tmp,{org_id});
    idx_tmp = ismember(combinedstring,cmbnstr_tmp);
    whole_ec_tmp = whole_ec(~idx_tmp);
    whole_substrate_tmp = whole_substrate(~idx_tmp);
    whole_org_tmp = whole_org(~idx_tmp);
    whole_kcat_tmp = whole_kcat(~idx_tmp);
    
    if ~ismember(ec_tmp,whole_ec_tmp)
        kcatvalue = 0;
    else
        sublist = whole_substrate_tmp(ismember(whole_ec_tmp,ec_tmp));
        orglist = whole_org_tmp(ismember(whole_ec_tmp,ec_tmp));
        kcatlist = whole_kcat_tmp(ismember(whole_ec_tmp,ec_tmp));

        match_idx_sub = ismember(sublist,sub_tmp);
        match_idx_org = ismember(orglist,org_id);
        match_idx_combined = match_idx_sub & match_idx_org;

        if any(match_idx_combined) % step 1
            kcatvalue = kcatlist(match_idx_combined);
        else
            if any(match_idx_org) % step 2
                kcatvalue = median(kcatlist(match_idx_org));
            else
                if any(match_idx_sub) % step 3
                    kcatvalue = median(kcatlist(match_idx_sub));
                else % step 4
                    kcatvalue = median(kcatlist);
                end
            end
        end
        kcat_b(i) = kcatvalue;
    end
end

% Method C:
% 1. if both organism and substrate matched, then choose the kcat
% 2. if the match in step 1 is not found, then choose median of the
%    organism with all substrates.
% 3. if the match in step 2 is not found, then choose median of the EC
%    number.

kcat_c = zeros(length(sc_data.ec),1); % plan B kcats

for i = 1:length(sc_data.ec)
    ec_tmp = sc_data.ec(i);
    sub_tmp = sc_data.sub(i);
    cmbnstr_tmp = strcat(ec_tmp,sub_tmp,{org_id});
    idx_tmp = ismember(combinedstring,cmbnstr_tmp);
    whole_ec_tmp = whole_ec(~idx_tmp);
    whole_substrate_tmp = whole_substrate(~idx_tmp);
    whole_org_tmp = whole_org(~idx_tmp);
    whole_kcat_tmp = whole_kcat(~idx_tmp);
    
    if ~ismember(ec_tmp,whole_ec_tmp)
        kcatvalue = 0;
    else
        sublist = whole_substrate_tmp(ismember(whole_ec_tmp,ec_tmp));
        orglist = whole_org_tmp(ismember(whole_ec_tmp,ec_tmp));
        kcatlist = whole_kcat_tmp(ismember(whole_ec_tmp,ec_tmp));

        match_idx_sub = ismember(sublist,sub_tmp);
        match_idx_org = ismember(orglist,org_id);
        match_idx_combined = match_idx_sub & match_idx_org;

        if any(match_idx_combined) % step 1
            kcatvalue = kcatlist(match_idx_combined);
        else
            if any(match_idx_org) % step 2
                kcatvalue = median(kcatlist(match_idx_org));
            else % step 3
                kcatvalue = median(kcatlist);
            end
        end
        kcat_c(i) = kcatvalue;
    end
end


nonzero_idx = kcat_a ~= 0 & kcat_b ~= 0 & kcat_c ~= 0;
kcatsc = sc_data.kcat(nonzero_idx)/3600;
kcata = kcat_a(nonzero_idx)/3600;
kcatb = kcat_b(nonzero_idx)/3600;
kcatc = kcat_c(nonzero_idx)/3600;
N_samples = length(kcatsc);

kcatsc = log10(kcatsc);
kcata = log10(kcata);
kcatb = log10(kcatb);
kcatc = log10(kcatc);

rmsea = sqrt(sum((kcatsc - kcata).^2)/N_samples);
rmseb = sqrt(sum((kcatsc - kcatb).^2)/N_samples);
rmsec = sqrt(sum((kcatsc - kcatc).^2)/N_samples);
[Ra,Pa] = corrcoef(kcatsc,kcata);
[Rb,Pb] = corrcoef(kcatsc,kcatb);
[Rc,Pc] = corrcoef(kcatsc,kcatc);

figure('Name','4');
hold on;
box on;
line([-7 7],[-7 7],'Color','k','LineStyle','-.');
scatter(kcatsc,kcata,'filled','MarkerFaceAlpha',0.5);
xlim([-7 7]);
ylim([-7 7]);
set(gca,'FontSize',9,'FontName','Helvetica');
ylabel('log10 (collected kcat)','FontSize',12,'FontName','Helvetica');
xlabel('log10 (matched kcat)','FontSize',12,'FontName','Helvetica');
title('Method A','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(rmsea,2))];
str_r = ['R = ',num2str(round(Ra(1,2),2))];
str_n = ['N = ',num2str(N_samples)];
text(-6.5,6,str_rmse,'FontSize',14,'FontName','Helvetica');
text(-6.5,4.5,str_r,'FontSize',14,'FontName','Helvetica');
text(1.5,-6,str_n,'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[0 0 230 230]);
set(gca,'position',[0.15 0.15 0.7 0.7]);

figure('Name','5');
hold on;
box on;
line([-7 7],[-7 7],'Color','k','LineStyle','-.');
scatter(kcatsc,kcatb,'filled','MarkerFaceAlpha',0.5);
xlim([-7 7]);
ylim([-7 7]);
set(gca,'FontSize',9,'FontName','Helvetica');
ylabel('log10 (collected kcat)','FontSize',12,'FontName','Helvetica');
xlabel('log10 (matched kcat)','FontSize',12,'FontName','Helvetica');
title('Method B','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(rmseb,2))];
str_r = ['R = ',num2str(round(Rb(1,2),2))];
str_n = ['N = ',num2str(N_samples)];
text(-6.5,6,str_rmse,'FontSize',14,'FontName','Helvetica');
text(-6.5,4.5,str_r,'FontSize',14,'FontName','Helvetica');
text(1.5,-6,str_n,'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[250 0 230 230]);
set(gca,'position',[0.15 0.15 0.7 0.7]);

figure('Name','6');
hold on;
box on;
line([-7 7],[-7 7],'Color','k','LineStyle','-.');
scatter(kcatsc,kcatc,'filled','MarkerFaceAlpha',0.5);
xlim([-7 7]);
ylim([-7 7]);
set(gca,'FontSize',9,'FontName','Helvetica');
ylabel('log10 (collected kcat)','FontSize',12,'FontName','Helvetica');
xlabel('log10 (matched kcat)','FontSize',12,'FontName','Helvetica');
title('Method C','FontSize',14,'FontName','Helvetica');
str_rmse = ['RMSE = ',num2str(round(rmsec,2))];
str_r = ['R = ',num2str(round(Rc(1,2),2))];
str_n = ['N = ',num2str(N_samples)];
text(-6.5,6,str_rmse,'FontSize',14,'FontName','Helvetica');
text(-6.5,4.5,str_r,'FontSize',14,'FontName','Helvetica');
text(1.5,-6,str_n,'FontSize',14,'FontName','Helvetica');
set(gcf,'position',[500 0 230 230]);
set(gca,'position',[0.15 0.15 0.7 0.7]);
