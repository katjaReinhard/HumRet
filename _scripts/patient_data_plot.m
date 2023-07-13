pathcode = 'C:\Users\Lenovo2020\OneDrive - imec\HumRet\_scripts';
pathSave = 'C:\Users\Lenovo2020\OneDrive - imec\HumanRetina\PlosOne\resubmission';

addpath(genpath(pathcode))
%%
patient_data = [7	1	18	nan	2	15 19;7	1	20	nan	4	15 29;...
17	1	5	0	0	0 5;7	0	15	9	0	1 20;...
10	0	22	21	9	8 29;10	0	54	10	1	3 60;...
7	0	86	68	11	42 88;7	1	23	15	0	26 35;...
8	1	13	13	10	8 16;10	1	29	6	0	3 30]; %ischemia, disease, dg, chrp, bar flash, total
%%
cols = [27,158,119;217,95,2;117,112,179]; cols = cols./255;
figure
subplot(2,2,1)
dis = find(patient_data(:,2)==1);
nodis = find(patient_data(:,2)==0);
scatter(patient_data(nodis,1),patient_data(nodis,7),45,'k','linewidth',2);
hold on
scatter(patient_data(dis,1),patient_data(dis,7),35,[0.6 0.6 0.6],'filled');
axis([0 20 0 100])
xlabel('ischemia (min)')
ylabel('# responding units')

subplot(2,2,3)


P = polyfit(patient_data(:,1),patient_data(:,3),1);
yfit = P(1)*patient_data(:,1)+P(2);
hold on
plot(patient_data(:,1),yfit,'-','color',cols(1,:))

tmp = find(~isnan(patient_data(:,4)));
P = polyfit(patient_data(tmp,1),patient_data(tmp,4),1);
yfit = P(1)*patient_data(tmp,1)+P(2);
hold on
plot(patient_data(tmp,1),yfit,'-','color',cols(2,:))

P = polyfit(patient_data(:,1),patient_data(:,6),1);
yfit = P(1)*patient_data(:,1)+P(2);
hold on
plot(patient_data(:,1),yfit,'-','color',cols(3,:))

scatter(patient_data(:,1),patient_data(:,3),45,cols(1,:),'linewidth',1);
s1 = scatter(patient_data(:,1),patient_data(:,4),35,cols(2,:),'filled');
scatter(patient_data(:,1),patient_data(:,6),20,cols(3,:),'filled');
alpha(s1,0.5)
xlabel('ischemia (min)')
ylabel('# responding units')
axis([0 20 0 100])
legend({'drifting-grating','chirp','flash'})

subplot(2,2,4)
dis = find(patient_data(:,2)==1);
nodis = find(patient_data(:,2)==0);
data_dis = patient_data(dis,7);
data_nodis = patient_data(nodis,7);

scatter(repmat(1,1,length(data_nodis)),data_nodis,15,[0.2 0.2 0.2]);
% plotSpread(data_nodis,'distributionIdx',repmat(1,length(data_nodis),1),'distributionMarkers',{'.'},'distributionColors',{[0.2 0.2 0.2]})
 hold on
 p1 = plot([1-0.3 1+0.3],[median(data_nodis) median(data_nodis)],'-','color',[0.2 0.2 0.2],'linewidth',2)
scatter(repmat(7,1,length(data_dis)),data_dis,15,[0.6 0.6 0.6],'markerfacealpha',0.5);
%  plotSpread(data_dis,'distributionIdx',repmat(7,length(data_dis),1),'distributionMarkers',{'.'},'distributionColors',{[0.6 0.6 0.6]})
 hold on
 plot([7-0.3 7+0.3],[median(data_dis) median(data_dis)],'-','color',[0.6 0.6 0.6],'linewidth',2)

 data_dis = patient_data(dis,3);
data_nodis = patient_data(nodis,3);
scatter(repmat(2,1,length(data_nodis)),data_nodis,15,cols(1,:));
% plotSpread(data_nodis,'distributionIdx',repmat(1,length(data_nodis),1),'distributionMarkers',{'.'},'distributionColors',{[0.2 0.2 0.2]})
 hold on
 p1 = plot([2-0.3 2+0.3],[median(data_nodis) median(data_nodis)],'-','color',cols(1,:),'linewidth',2)
scatter(repmat(8,1,length(data_dis)),data_dis,15,[cols(1,:)],'markerfacealpha',0.5);
%  plotSpread(data_dis,'distributionIdx',repmat(7,length(data_dis),1),'distributionMarkers',{'.'},'distributionColors',{[0.6 0.6 0.6]})
 hold on
 plot([8-0.3 8+0.3],[median(data_dis) median(data_dis)],'-','color',[cols(1,:) 0.5],'linewidth',2)

 
 data_dis = patient_data(dis,4);
data_nodis = patient_data(nodis,4);
scatter(repmat(3,1,length(data_nodis)),data_nodis,15,cols(2,:));
% plotSpread(data_nodis,'distributionIdx',repmat(1,length(data_nodis),1),'distributionMarkers',{'.'},'distributionColors',{[0.2 0.2 0.2]})
 hold on
 p1 = plot([3-0.3 3+0.3],[median(data_nodis) median(data_nodis)],'-','color',cols(2,:),'linewidth',2)
scatter(repmat(9,1,length(data_dis)),data_dis,15,[cols(2,:)],'markerfacealpha',0.5);
%  plotSpread(data_dis,'distributionIdx',repmat(7,length(data_dis),1),'distributionMarkers',{'.'},'distributionColors',{[0.6 0.6 0.6]})
 hold on
 plot([9-0.3 9+0.3],[nanmedian(data_dis) nanmedian(data_dis)],'-','color',[cols(2,:) 0.5],'linewidth',2)

 
  data_dis = patient_data(dis,6);
data_nodis = patient_data(nodis,6);
scatter(repmat(4,1,length(data_nodis)),data_nodis,15,cols(3,:));
% plotSpread(data_nodis,'distributionIdx',repmat(1,length(data_nodis),1),'distributionMarkers',{'.'},'distributionColors',{[0.2 0.2 0.2]})
 hold on
 p1 = plot([4-0.3 4+0.3],[median(data_nodis) median(data_nodis)],'-','color',cols(3,:),'linewidth',2)
scatter(repmat(10,1,length(data_dis)),data_dis,15,[cols(3,:)],'markerfacealpha',0.5);
%  plotSpread(data_dis,'distributionIdx',repmat(7,length(data_dis),1),'distributionMarkers',{'.'},'distributionColors',{[0.6 0.6 0.6]})
 hold on
 plot([10-0.3 10+0.3],[median(data_dis) median(data_dis)],'-','color',[cols(3,:) 0.5],'linewidth',2)

 axis([0 11 0 100])
 set(gca,'xtick',[2.5 8.5],'xticklabel',{'healthy','with cond.'})
 ylabel('# responding units')
 
%  legend([p1 p2 p3 p4],{'total','drifting-grating','chirp','flash'},'fontname','Arial','fontsize',7)
 
 
subplot(2,2,2)
dis = find(patient_data(:,2)==1);
nodis = find(patient_data(:,2)==0);
data_dis = patient_data(dis,1);
data_nodis = patient_data(nodis,1);

plotSpread(data_nodis,'distributionIdx',repmat(1,length(data_nodis),1),'distributionMarkers',{'.'},'distributionColors',{[0.2 0.2 0.2]})
 hold on
 p1 = plot([1-0.3 1+0.3],[median(data_nodis) median(data_nodis)],'-','color',[0.2 0.2 0.2],'linewidth',2)
plotSpread(data_dis,'distributionIdx',repmat(2,length(data_dis),1),'distributionMarkers',{'.'},'distributionColors',{[0.6 0.6 0.6]})
 hold on
 plot([2-0.3 2+0.3],[median(data_dis) median(data_dis)],'-','color',[0.6 0.6 0.6],'linewidth',2)
 axis([0 3 0 20])
 set(gca,'xtick',[1 2],'xticklabel',{'healthy','with cond.'})
 ylabel('ischemia (min)')

% saveas(gcf,fullfile(pathSave,'data_per_donor.svg'))
