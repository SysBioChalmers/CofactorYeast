
base = [1;2;3;nan;5;6;7];
scatterdev = [base-0.3;base;base+0.3];


[num, txt, ~] = xlsread('pHCA.xlsx','QL01');
label_BPS = txt(3:5,2);
label_BPS = strrep(label_BPS,' μM','');

OD_12 = num(1:3,1:3);
OD_20 = num(4:6,1:3);
pHCA_12 = num(1:3,4:6);
pHCA_20 = num(4:6,4:6);
rg_12 = num(1:3,7:9);
rg_20 = num(4:6,7:9);
PperCDW_12 = pHCA_12./(OD_12*0.65);
PperCDW_20 = pHCA_20./(OD_20*0.65);


OD_12_mean = mean(OD_12,2);
OD_20_mean = mean(OD_20,2);
pHCA_12_mean = mean(pHCA_12,2);
pHCA_20_mean = mean(pHCA_20,2);
rg_12_mean = mean(rg_12,2);
rg_20_mean = mean(rg_20,2);
PperCDW_12_mean = mean(PperCDW_12,2);
PperCDW_20_mean = mean(PperCDW_20,2);


OD_12_std = std(OD_12,1,2);
OD_20_std = std(OD_20,1,2);
pHCA_12_std = std(pHCA_12,1,2);
pHCA_20_std = std(pHCA_20,1,2);
rg_12_std = std(rg_12,1,2);
rg_20_std = std(rg_20,1,2);
PperCDW_12_std = std(PperCDW_12,1,2);
PperCDW_20_std = std(PperCDW_20,1,2);

figure('Name','main');
data_mean = [PperCDW_12_mean;nan;PperCDW_20_mean];
data_std = [PperCDW_12_std;nan;PperCDW_20_std];
data_raw = [PperCDW_12;[nan nan nan];PperCDW_20];
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
b.CData(5,:) = [255,255,255]/255;
b.CData(6,:) = [253,208,162]/255;
b.CData(7,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),8,'k');
s.LineWidth = 0.5;
xticks([2 6]);
xlim([0.25 7.75]);
xticklabels([{'12 h'};{'20 h'}]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylim([0 42]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
[~,p3] = ttest2(data_raw(5,:),data_raw(6,:),'Vartype','unequal');
[~,p4] = ttest2(data_raw(5,:),data_raw(7,:),'Vartype','unequal');
if p1 < 0.05 && p1 >= 0.01
    p1lable = '*';
    p1sz = 14;
elseif p1 < 0.01 && p1 >= 0.001
    p1lable = '**';
    p1sz = 14;
elseif p1 < 0.001
    p1lable = '***';
    p1sz = 14;
else
    p1lable = 'n.s.';
    p1sz = 7;
end
if p2 < 0.05 && p2 >= 0.01
    p2lable = '*';
    p2sz = 14;
elseif p2 < 0.01 && p2 >= 0.001
    p2lable = '**';
    p2sz = 14;
elseif p2 < 0.001
    p2lable = '***';
    p2sz = 14;
else
    p2lable = 'n.s.';
    p2sz = 7;
end
if p3 < 0.05 && p3 >= 0.01
    p3lable = '*';
    p3sz = 14;
elseif p3 < 0.01 && p3 >= 0.001
    p3lable = '**';
    p3sz = 14;
elseif p3 < 0.001
    p3lable = '***';
    p3sz = 14;
else
    p3lable = 'n.s.';
    p3sz = 7;
end
if p4 < 0.05 && p4 >= 0.01
    p4lable = '*';
    p4sz = 14;
elseif p4 < 0.01 && p4 >= 0.001
    p4lable = '**';
    p4sz = 14;
elseif p4 < 0.001
    p4lable = '***';
    p4sz = 14;
else
    p4lable = 'n.s.';
    p4sz = 7;
end
text(1.5,32,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,33,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(5.5,32,p3lable,'Color','black','FontSize',p3sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(6,33,p4lable,'Color','black','FontSize',p4sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],32*ones(1,2),'Color','black');
line([1 3],37*ones(1,2),'Color','black');
line([5 6],32*ones(1,2),'Color','black');
line([5 7],37*ones(1,2),'Color','black');
xlabel('Fermentation time','FontSize',7,'FontName','Helvetica');
ylabel('mg p-HCA/g CDW','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[500 400 110 110]);


figure('Name','SI');
subplot(1,3,1);
data_mean = [pHCA_12_mean;nan;pHCA_20_mean];
data_std = [pHCA_12_std;nan;pHCA_20_std];
data_raw = [pHCA_12;[nan nan nan];pHCA_20];
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
b.CData(5,:) = [255,255,255]/255;
b.CData(6,:) = [253,208,162]/255;
b.CData(7,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),8,'k');
s.LineWidth = 0.5;
xticks([2 6]);
xlim([0.25 7.75]);
xticklabels([{'12 h'};{'20 h'}]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylim([0 80]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
[~,p3] = ttest2(data_raw(5,:),data_raw(6,:),'Vartype','unequal');
[~,p4] = ttest2(data_raw(5,:),data_raw(7,:),'Vartype','unequal');
if p1 < 0.05 && p1 >= 0.01
    p1lable = '*';
    p1sz = 14;
elseif p1 < 0.01 && p1 >= 0.001
    p1lable = '**';
    p1sz = 14;
elseif p1 < 0.001
    p1lable = '***';
    p1sz = 14;
else
    p1lable = 'n.s.';
    p1sz = 7;
end
if p2 < 0.05 && p2 >= 0.01
    p2lable = '*';
    p2sz = 14;
elseif p2 < 0.01 && p2 >= 0.001
    p2lable = '**';
    p2sz = 14;
elseif p2 < 0.001
    p2lable = '***';
    p2sz = 14;
else
    p2lable = 'n.s.';
    p2sz = 7;
end
if p3 < 0.05 && p3 >= 0.01
    p3lable = '*';
    p3sz = 14;
elseif p3 < 0.01 && p3 >= 0.001
    p3lable = '**';
    p3sz = 14;
elseif p3 < 0.001
    p3lable = '***';
    p3sz = 14;
else
    p3lable = 'n.s.';
    p3sz = 7;
end
if p4 < 0.05 && p4 >= 0.01
    p4lable = '*';
    p4sz = 14;
elseif p4 < 0.01 && p4 >= 0.001
    p4lable = '**';
    p4sz = 14;
elseif p4 < 0.001
    p4lable = '***';
    p4sz = 14;
else
    p4lable = 'n.s.';
    p4sz = 7;
end
text(1.5,15,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,25,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(5.5,55,p3lable,'Color','black','FontSize',p3sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(6,65,p4lable,'Color','black','FontSize',p4sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],22*ones(1,2),'Color','black');
line([1 3],32*ones(1,2),'Color','black');
line([5 6],62*ones(1,2),'Color','black');
line([5 7],72*ones(1,2),'Color','black');
xlabel('Fermentation time','FontSize',7,'FontName','Helvetica');
ylabel('p-HCA titer (mg/L)','FontSize',7,'FontName','Helvetica');

subplot(1,3,2);
data_mean = [OD_12_mean;nan;OD_20_mean];
data_std = [OD_12_std;nan;OD_20_std];
data_raw = [OD_12;[nan nan nan];OD_20];
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
b.CData(5,:) = [255,255,255]/255;
b.CData(6,:) = [253,208,162]/255;
b.CData(7,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),8,'k');
s.LineWidth = 0.5;
xticks([2 6]);
xlim([0.25 7.75]);
xticklabels([{'12 h'};{'20 h'}]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylim([0 5]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
[~,p3] = ttest2(data_raw(5,:),data_raw(6,:),'Vartype','unequal');
[~,p4] = ttest2(data_raw(5,:),data_raw(7,:),'Vartype','unequal');
if p1 < 0.05 && p1 >= 0.01
    p1lable = '*';
    p1sz = 14;
elseif p1 < 0.01 && p1 >= 0.001
    p1lable = '**';
    p1sz = 14;
elseif p1 < 0.001
    p1lable = '***';
    p1sz = 14;
else
    p1lable = 'n.s.';
    p1sz = 7;
end
if p2 < 0.05 && p2 >= 0.01
    p2lable = '*';
    p2sz = 14;
elseif p2 < 0.01 && p2 >= 0.001
    p2lable = '**';
    p2sz = 14;
elseif p2 < 0.001
    p2lable = '***';
    p2sz = 14;
else
    p2lable = 'n.s.';
    p2sz = 7;
end
if p3 < 0.05 && p3 >= 0.01
    p3lable = '*';
    p3sz = 14;
elseif p3 < 0.01 && p3 >= 0.001
    p3lable = '**';
    p3sz = 14;
elseif p3 < 0.001
    p3lable = '***';
    p3sz = 14;
else
    p3lable = 'n.s.';
    p3sz = 7;
end
if p4 < 0.05 && p4 >= 0.01
    p4lable = '*';
    p4sz = 14;
elseif p4 < 0.01 && p4 >= 0.001
    p4lable = '**';
    p4sz = 14;
elseif p4 < 0.001
    p4lable = '***';
    p4sz = 14;
else
    p4lable = 'n.s.';
    p4sz = 7;
end
text(1.5,1.4,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,1.55,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(5.5,3.6,p3lable,'Color','black','FontSize',p3sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(6,3.75,p4lable,'Color','black','FontSize',p4sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],1.4*ones(1,2),'Color','black');
line([1 3],2*ones(1,2),'Color','black');
line([5 6],3.6*ones(1,2),'Color','black');
line([5 7],4.2*ones(1,2),'Color','black');
xlabel('Fermentation time','FontSize',7,'FontName','Helvetica');
ylabel('OD600','FontSize',7,'FontName','Helvetica');

subplot(1,3,3);
data_mean = [rg_12_mean;nan;rg_20_mean];
data_std = [rg_12_std;nan;rg_20_std];
data_raw = [rg_12;[nan nan nan];rg_20];
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
b.CData(5,:) = [255,255,255]/255;
b.CData(6,:) = [253,208,162]/255;
b.CData(7,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),8,'k');
s.LineWidth = 0.5;
xticks([2 6]);
xlim([0.25 7.75]);
xticklabels([{'12 h'};{'20 h'}]);
set(gca,'FontSize',6,'FontName','Helvetica');
ylim([0 20]);

xlabel('Fermentation time','FontSize',7,'FontName','Helvetica');
ylabel('Residual glucose (g/L)','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[500 100 360 110]);


