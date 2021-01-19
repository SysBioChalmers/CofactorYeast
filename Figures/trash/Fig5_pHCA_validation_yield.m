scatterdev = repmat([1;2;3],3,1)+[repmat(-0.3,3,1);zeros(3,1);repmat(0.3,3,1)];

[num, txt, ~] = xlsread('pHCA.xlsx','QL01');
label_BPS = txt(3:5,2);
label_BPS = strrep(label_BPS,' μM','');

OD_12 = num(1:3,1:3);
OD_20 = num(4:6,1:3);
pHCA_12 = num(1:3,4:6);
pHCA_20 = num(4:6,4:6);
rg_12 = num(1:3,7:9);
rg_20 = num(4:6,7:9);

yb_12 = OD_12*0.65./(20-rg_12);
yb_20 = OD_20*0.65./(20-rg_20);
yp_12 = pHCA_12./(20-rg_12)/1000;
yp_20 = pHCA_20./(20-rg_20)/1000;


yb_12_mean = mean(yb_12,2);
yb_20_mean = mean(yb_20,2);
yp_12_mean = mean(yp_12,2);
yp_20_mean = mean(yp_20,2);

yb_12_std = std(yb_12,1,2);
yb_20_std = std(yb_20,1,2);
yp_12_std = std(yp_12,1,2);
yp_20_std = std(yp_20,1,2);

figure('Name','QL01');
subplot(2,2,1);
data_mean = yb_12_mean;
data_std = yb_12_std;
data_raw = yb_12;
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
ylim([0 0.2]);
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
text(1.5,0.13,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,0.16,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],0.13*ones(1,2)*1,'Color','black');
line([1 3],0.16*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('Biomass yield (g/g)','FontSize',8,'FontName','Helvetica');
title('12 h','FontSize',8,'FontName','Helvetica');

subplot(2,2,2);
data_mean = yb_20_mean;
data_std = yb_20_std;
data_raw = yb_20;
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
ylim([0 0.2]);
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
text(1.5,0.13,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,0.16,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],0.13*ones(1,2)*1,'Color','black');
line([1 3],0.16*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('Biomass yield (g/g)','FontSize',8,'FontName','Helvetica');
title('20 h','FontSize',8,'FontName','Helvetica');

subplot(2,2,3);
data_mean = yp_12_mean;
data_std = yp_12_std;
data_raw = yp_12;
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
ylim([0 4.5e-3]);
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
text(1.5,3e-3,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,3.7e-3,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],3e-3*ones(1,2)*1.13,'Color','black');
line([1 3],3.7e-3*ones(1,2)*1.1,'Color','black');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('p-HCA yield (g/g)','FontSize',8,'FontName','Helvetica');
title('12 h','FontSize',8,'FontName','Helvetica');

subplot(2,2,4);
data_mean = yp_20_mean;
data_std = yp_20_std;
data_raw = yp_20;
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
ylim([0 4.5e-3]);
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
text(1.5,3e-3,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
text(2,3.7e-3,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
line([1 2],3e-3*ones(1,2)*1.13,'Color','black');
line([1 3],3.7e-3*ones(1,2)*1.1,'Color','black');
ylabel('p-HCA yield (g/g)','FontSize',8,'FontName','Helvetica');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
title('20 h','FontSize',8,'FontName','Helvetica');

set(gcf,'position',[200 400 220 250]);


[num, txt, ~] = xlsread('pHCA.xlsx','IMX581');
label_BPS = txt(3:5,2);
label_BPS = strrep(label_BPS,' μM','');

OD_12 = num(1:3,1:3);
OD_20 = num(4:6,1:3);

yb_12_mean = mean(OD_12,2);
yb_20_mean = mean(OD_20,2);

yb_12_std = std(OD_12,1,2);
yb_20_std = std(OD_20,1,2);

figure('Name','IMX581');
subplot(1,2,1);
data_mean = yb_12_mean;
data_std = yb_12_std;
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
text(1.5,data_mean(1)*1.3,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica');
text(2,data_mean(1)*1.5,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('OD600','FontSize',8,'FontName','Helvetica');
title('12 h','FontSize',8,'FontName','Helvetica');

subplot(1,2,2);
data_mean = yb_20_mean;
data_std = yb_20_std;
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
text(1.5,data_mean(1)*1.3,p1lable,'Color','black','FontSize',p1sz,'FontName','Helvetica');
text(2,data_mean(1)*1.5,p2lable,'Color','black','FontSize',p2sz,'FontName','Helvetica');
xlabel('BPS concentration (μM)','FontSize',8,'FontName','Helvetica');
ylabel('OD600','FontSize',8,'FontName','Helvetica');
title('20 h','FontSize',8,'FontName','Helvetica');

set(gcf,'position',[200 100 220 125]);


