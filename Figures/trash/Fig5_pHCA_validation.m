scatterdev = repmat([1;2;3],3,1)+[repmat(-0.3,3,1);zeros(3,1);repmat(0.3,3,1)];

[num, txt, ~] = xlsread('pHCA.xlsx','QL01');
label_BPS = txt(3:5,2);
label_BPS = strrep(label_BPS,' μM','');

OD_12 = num(1:3,1:3);
OD_20 = num(4:6,1:3);
pHCA_12 = num(1:3,4:6);
pHCA_20 = num(4:6,4:6);

OD_12_mean = mean(OD_12,2);
OD_20_mean = mean(OD_20,2);
pHCA_12_mean = mean(pHCA_12,2);
pHCA_20_mean = mean(pHCA_20,2);

OD_12_std = std(OD_12,1,2);
OD_20_std = std(OD_20,1,2);
pHCA_12_std = std(pHCA_12,1,2);
pHCA_20_std = std(pHCA_20,1,2);

figure('Name','QL01');
subplot(2,2,1);
data_mean = OD_12_mean;
data_std = OD_12_std;
data_raw = OD_12;
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),15,'k');
s.LineWidth = 0.5;
xticks(1:1:3);
xlim([0.25 3.75]);
xticklabels(label_BPS);
set(gca,'FontSize',7,'FontName','Helvetica');
ylim([0 round(1.6*data_mean(1),2)]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
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
text(1.5,data_mean(1)*1.2,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,data_mean(1)*1.3,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],data_mean(1)*1.2*ones(1,2)*1,'Color','black');
line([1 3],data_mean(1)*1.3*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('OD600','FontSize',8,'FontName','Helvetica');
title('12 h','FontSize',8,'FontName','Helvetica');

subplot(2,2,2);
data_mean = OD_20_mean;
data_std = OD_20_std;
data_raw = OD_20;
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),15,'k');
s.LineWidth = 0.5;
set(gca,'FontSize',7,'FontName','Helvetica');
xticks(1:1:3);
xlim([0.25 3.75]);
xticklabels(label_BPS);
set(gca,'FontSize',7,'FontName','Helvetica');
ylim([0 round(1.6*data_mean(1),2)]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
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
text(1.5,data_mean(1)*1.2,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,data_mean(1)*1.3,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],data_mean(1)*1.2*ones(1,2)*1,'Color','black');
line([1 3],data_mean(1)*1.3*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('OD600','FontSize',8,'FontName','Helvetica');
title('20 h','FontSize',8,'FontName','Helvetica');

subplot(2,2,3);
data_mean = pHCA_12_mean;
data_std = pHCA_12_std;
data_raw = pHCA_12;
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),15,'k');
s.LineWidth = 0.5;
set(gca,'FontSize',7,'FontName','Helvetica');
xticks(1:1:3);
xlim([0.25 3.75]);
xticklabels(label_BPS);
set(gca,'FontSize',7,'FontName','Helvetica');
ylim([0 round(1.6*data_mean(1),2)]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
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
text(1.5,data_mean(1)*1.05,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,data_mean(1)*1.3,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],data_mean(1)*1.05*ones(1,2)*1.13,'Color','black');
line([1 3],data_mean(1)*1.3*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('p-HCA titer (mg/L)','FontSize',8,'FontName','Helvetica');
title('12 h','FontSize',8,'FontName','Helvetica');

subplot(2,2,4);
data_mean = pHCA_20_mean;
data_std = pHCA_20_std;
data_raw = pHCA_20;
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),15,'k');
s.LineWidth = 0.5;
set(gca,'FontSize',7,'FontName','Helvetica');
xticks(1:1:3);
xlim([0.25 3.75]);
xticklabels(label_BPS);
set(gca,'FontSize',7,'FontName','Helvetica');
ylim([0 round(1.6*data_mean(1),2)]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
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
text(1.5,data_mean(1)*1.05,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,data_mean(1)*1.3,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],data_mean(1)*1.05*ones(1,2)*1.13,'Color','black');
line([1 3],data_mean(1)*1.3*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('p-HCA titer (mg/L)','FontSize',8,'FontName','Helvetica');
title('20 h','FontSize',8,'FontName','Helvetica');

set(gcf,'position',[200 400 220 250]);


[num, txt, ~] = xlsread('pHCA.xlsx','IMX581');
label_BPS = txt(3:5,2);
label_BPS = strrep(label_BPS,' μM','');

OD_12 = num(1:3,1:3);
OD_20 = num(4:6,1:3);

OD_12_mean = mean(OD_12,2);
OD_20_mean = mean(OD_20,2);

OD_12_std = std(OD_12,1,2);
OD_20_std = std(OD_20,1,2);

figure('Name','IMX581');
subplot(1,2,1);
data_mean = OD_12_mean;
data_std = OD_12_std;
data_raw = OD_12;
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),15,'k');
s.LineWidth = 0.5;
set(gca,'FontSize',7,'FontName','Helvetica');
xticks(1:1:3);
xlim([0.25 3.75]);
xticklabels(label_BPS);
set(gca,'FontSize',7,'FontName','Helvetica');
ylim([0 round(1.6*data_mean(1),2)]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
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
text(1.5,data_mean(1)*1.2,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,data_mean(1)*1.3,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],data_mean(1)*1.2*ones(1,2)*1,'Color','black');
line([1 3],data_mean(1)*1.3*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('OD600','FontSize',8,'FontName','Helvetica');
title('12 h','FontSize',8,'FontName','Helvetica');

subplot(1,2,2);
data_mean = OD_20_mean;
data_std = OD_20_std;
data_raw = OD_20;
hold on;
b = bar(data_mean,0.7,'LineWidth',0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [255,255,255]/255;
b.CData(2,:) = [253,208,162]/255;
b.CData(3,:) = [242,94,13]/255;
errorbar(data_mean,data_std,'k','Marker','none','LineStyle','none','LineWidth',0.5,'CapSize',5);
s = scatter(scatterdev,data_raw(:),15,'k');
s.LineWidth = 0.5;
set(gca,'FontSize',7,'FontName','Helvetica');
xticks(1:1:3);
xlim([0.25 3.75]);
xticklabels(label_BPS);
set(gca,'FontSize',7,'FontName','Helvetica');
ylim([0 round(1.6*data_mean(1),2)]);
[~,p1] = ttest2(data_raw(1,:),data_raw(2,:),'Vartype','unequal');
[~,p2] = ttest2(data_raw(1,:),data_raw(3,:),'Vartype','unequal');
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
text(1.5,data_mean(1)*1.05,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,data_mean(1)*1.3,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],data_mean(1)*1.05*ones(1,2)*1.13,'Color','black');
line([1 3],data_mean(1)*1.3*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('OD600','FontSize',8,'FontName','Helvetica');
title('20 h','FontSize',8,'FontName','Helvetica');

set(gcf,'position',[200 100 220 125]);


