clear
close all
clc
%% pathes etc
path1='C:\Users\katja\OneDrive - imec\HumRet';
pathcode = 'C:\Users\katja\OneDrive - imec\HumRet\_scripts';
pathColor = 'F:\LabCode\Matlab\Plots\cmocean_v1.4\cmocean';
pathSave = 'C:\Users\katja\OneDrive - imec\HumanRetina\PlosOne';

%% load
addpath(genpath(pathcode))
addpath(genpath(pathColor))
load(fullfile(path1,'h_responding'))
load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);
load(fullfile(path1,'abs_ID_exampleCells'))
%% variables
spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];

fsize1 = 7; fsize2 = 9; fsize3 = 10;

ctr_edg = [-1 1];
sp_edg = [0.6 4.2];

spat=[4000 2000 1000 500 200 100];
freq=[1 2 4 8 ];
%% colors
CC = 10;
order = [4 7 8 3 5 6 2 9 1 10];
CLUSTER_INFO = zeros(CC,20,5);%mean,median,min max, std

clrs = [166,206,227;31,120,180;178,223,138;51,160,44;251,154,153;227,26,28;...
    253,191,111;255,127,0;202,178,214;106,61,154];
clrs = clrs./255;
%% start figure
FIG = figure('position',[635.6667 821 810.0000 1140],'color','w')

save_png = 1;
%% cluster heatmaps -> report contrast preference
axes('position',[0.08 0.95 0.8 0.04])


doclust = 0;
if doclust == 1
    CLUST = [];
    EVAL = [];
    for k = 2:30
        [idx,C] = kmeans(normP,k,'Replicates',1000);
        eva1 = evalclusters(normP,idx,'CalinskiHarabasz');
        eva2 = evalclusters(normP,idx,'DaviesBouldin');
        CLUST = [CLUST idx];
        evals = [eva1.CriterionValues; eva2.CriterionValues];
        EVAL = [EVAL evals];
    end
    save(fullfile(path1,'cluster_DGheatmaps'),'CLUST','EVAL','normP')
else
    load(fullfile(path1,'cluster_DGheatmaps'))
end

CL = CLUST(:,CC-1);

% --- resort
CL2 = CL;
for o = 1:length(order)
    tmp = find(CL==order(o));
    CL2(tmp) = o;
end
CL = CL2;

%----

HEATMAP = [];
for c = 1:CC
    tmp  =find(CL==c);
    
    htmp = 1-normP(tmp,:); htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    
    HEATMAP = [HEATMAP htmp];
    if c < CC
        HEATMAP = [HEATMAP ones(4,1)];
    end
end

imagesc(HEATMAP)
colormap(gca,'gray')
axis off
box off

hold on
strt = 0.5;
for x = 1:CC
    rectangle('position',[strt 0.5 6 4],'facecolor','none','edgecolor',clrs(x,:))
    strt=strt+7;
end

cnt = 0;
for x = 1:6
    cnt = cnt + 1;
    tx = text(x,4.7,int2str(spat(cnt)),'fontname','Arial','fontsize',fsize1,'horizontalalignment','right');
    set(tx,'rotation',90)
end
tx = text(6.5,5.1,'\mum','fontname','Arial','fontsize',fsize1);


cnt = 0;
for x = 4:-1:1
    cnt = cnt + 1;
    tx = text(0.4,x,int2str(freq(cnt)),'fontname','Arial','fontsize',fsize1,'horizontalalignment','right');
end
text(-1.7,2.5,'Hz','fontname','Arial','fontsize',fsize1)

axes('position',[0.08 0.9 0.8 0.02])
for x = 1:CC
    plot(x,length(find(CL==x)),'o','markersize',5,'markerfacecolor',clrs(x,:),'markeredgecolor',clrs(x,:))
    hold on
    CLUSTER_INFO(x,1,1)=length(find(CL==x));%mean,median,min max, std
end
yy = get(gca,'ylim');
axis([0.5 CC+0.5 0 yy(2)+10])
box off
ylabel('# cells','fontname','Arial','fontsize',fsize1)
%% TF vs SF
load(fullfile(path1,'h_normPeaks'))
togoDG = togo;

axes('position',[0.05 0.275 0.15 0.15])
hold off
normP = normPeaks(:,:,2);

[cX,cY,oX,oY,dX,dY] = get_GaussHeatmaps(normPeaks);

% scatter(7-cX,5-cY,10,[0.7 0.7 0.7])%spatial,temporal
hold on
for cc = 1:CC
    tmp = find(CL==cc);
    mx = median(cX(tmp)); my = median(cY(tmp));
    scatter(7-cX(tmp),5-cY(tmp),10,clrs(cc,:))
    hold on
    %     plot(7-mx,5-my,'o','markerfacecolor',clrs(cc,:),'markeredgecolor',clrs(cc,:),'markersize',5)
end


set(gca,'plotboxaspectratio',[1 1 1])
set(gca,'xtick',1:6,'xticklabel',[4000 2000 1000 500 200 100],...
    'ytick',1:4,'yticklabel',[1 2 4 8])
xlabel('\mum (4B)','fontname','Arial','fontsize',fsize2)
ylabel('Hz (4C)','fontname','Arial','fontsize',fsize2)
axis([sp_edg 0.8 4.2])
xtickangle(90)
%% contrast vs preferred SF
axes('position',[0.29 0.275 0.15 0.15])
[chirpFFTnorm,chirpAmp,freqSmoothChrp,togoChrp]=get_ChirpFFTnorm(path1,info);

rangeF = length(chirpAmp)-chirpAmp(1);
steps = round(1:length(chirpAmp)/3:length(chirpAmp));
chk1=steps(1):steps(2)-1;
chk3=steps(3):length(chirpAmp);

cntrstIDX = [];
for ce = 1:size(chirpAmp,1)
    low = sum(abs(chirpAmp(ce,chk1)));
    high = sum(abs(chirpAmp(ce,chk3)));
    
    idx = (high-low)./(low+high);
    cntrstIDX = [cntrstIDX;idx];
end

comp = intersect(togoChrp,togoDG);
idd = []; idc = [];
for ii = 1:length(comp)
    tmp  = find(togoDG==comp(ii));
    idd = [idd;tmp];
    tmp  = find(togoChrp==comp(ii));
    idc = [idc;tmp];
end

% scatter(7-cX(idd),cntrstIDX(idc),10,[0.7 0.7 0.7]);
hold on
for cc = 1:CC
    tmp = find(CL==cc);
    tdg = intersect(idd,tmp);
    mx = median(cX(tdg));
    absid = togoDG(tmp);
    tmp2 = [];
    for aa = 1:length(absid)
        tm = find(togoChrp==absid(aa));
        tmp2 = [tmp2 tm];
    end
    tmp2 = intersect(idc,tmp2);
    my = median(cntrstIDX(tmp2));
    scatter(7-cX(tdg),cntrstIDX(tmp2),10,clrs(cc,:))
    %     plot(7-mx,my,'o','markerfacecolor',clrs(cc,:),'markeredgecolor',clrs(cc,:),'markersize',5)
end

set(gca,'plotboxaspectratio',[1 1 1])
set(gca,'xtick',1:6,'xticklabel',[4000 2000 1000 500 200 100],...
    'ytick',-1:0.5:1)
xlabel('\mum (4B)','fontname','Arial','fontsize',fsize2)
ylabel('contrast pref. index (4E)','fontname','Arial','fontsize',fsize2)
xtickangle(90)
axis([sp_edg ctr_edg])
%% contrast vs TF

axes('position',[0.51 0.275 0.15 0.15])
% scatter(5-cY(idd),cntrstIDX(idc),10,[0.7 0.7 0.7]);
hold on
for cc = 1:CC
    tmp = find(CL==cc);
    tdg = intersect(idd,tmp);
    mx = median(cY(tdg));
    absid = togoDG(tmp);
    tmp2 = [];
    for aa = 1:length(absid)
        tm = find(togoChrp==absid(aa));
        tmp2 = [tmp2 tm];
    end
    tmp2 = intersect(idc,tmp2);
    my = median(cntrstIDX(tmp2));
    scatter(5-cY(tdg),cntrstIDX(tmp2),10,clrs(cc,:))
    %     plot(5-mx,my,'o','markerfacecolor',clrs(cc,:),'markeredgecolor',clrs(cc,:),'markersize',5)
end

set(gca,'plotboxaspectratio',[1 1 1])
set(gca,'xtick',1:4,'xticklabel',[1 2 4 8],...
    'ytick',-1:0.5:1)
xlabel('Hz (4C)','fontname','Arial','fontsize',fsize2)
ylabel('contrast pref. index (4E)','fontname','Arial','fontsize',fsize2)
axis([0.8 4.2 ctr_edg])
%% contrast vs broadness of SF/TF
axes('position',[0.75 0.275 0.15 0.15])
% scatter((dX(idd)+dY(idd))./2,cntrstIDX(idc),10,[0.7 0.7 0.7]);
hold on
for cc = 1:CC
    tmp = find(CL==cc);
    tdg = intersect(idd,tmp);
    mx = median((dX(tdg)+dY(tdg))./2);
    absid = togoDG(tmp);
    tmp2 = [];
    for aa = 1:length(absid)
        tm = find(togoChrp==absid(aa));
        tmp2 = [tmp2 tm];
    end
    tmp2 = intersect(idc,tmp2);
    my = median(cntrstIDX(tmp2));
    scatter((dX(tdg)+dY(tdg))./2,cntrstIDX(tmp2),10,clrs(cc,:))
    %     plot(mx,my,'o','markerfacecolor',clrs(cc,:),'markeredgecolor',clrs(cc,:),'markersize',5)
end

axis([0 4 ctr_edg])
set(gca,'plotboxaspectratio',[1 1 1])
set(gca,'xtick',1:4,...
    'ytick',[-1:0.5:1])
xlabel('width of S-T resp. (4D)','fontname','Arial','fontsize',fsize2)
ylabel('contrast pref. index (4E)','fontname','Arial','fontsize',fsize2)
%% cluster info
SUPER_COLL = zeros(length(togoDG),20);
SUPER_COLL(:,1) = togoDG;
SUPER_COLL(:,2) = CL;
% --- spatial ----
SP = unique(spat); FR =unique(freq);
ax1 = axes('position',[0.08 0.78 0.8 0.035])
for cc = 1:CC
    tmp = find(CL==cc);
    dat = cX(tmp);
    dat2 = [];
    for d = 1:length(dat)
        now = dat(d);
        in = floor(now);
        rest = now-in;
        if rest >0
            slp = (SP(in+1)-SP(in));
            um = slp*rest+SP(in);
            dat2 = [dat2;um];
        else
            dat2 = [dat2;SP(in)];
        end
    end
    SUPER_COLL(tmp,3) = dat2;
    CLUSTER_INFO(cc,2,1)=mean(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,2,2)=median(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,2,3)=min(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,2,4)=max(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,2,5)=std(dat2);%mean,median,min max, std
    plotSpread(dat2,'distributionIdx',repmat(cc,length(dat2),1),'distributionMarkers',{'.'},'distributionColors',{clrs(cc,:)})
    %     scatter(repmat(cc,1,length(dat2)),dat2,10,clrs(cc,:))
    hold on
    plot([cc-0.3 cc+0.3],[median(dat2) median(dat2)],'-','color',clrs(cc,:),'linewidth',1)
    
end
set(gca,'xtick',1:10,'xticklabel',1:10,'fontname','Arial','fontsize',fsize2)
axis([0.5 CC+0.5 0 4000])
box off
ylabel('\mum','fontname','Arial','fontsize',fsize2)
axes('position',[0.08 0.826 0.8 0.005])
text(0,0,'spatial frequency preference','fontname','Arial','fontsize',fsize2)
axis off
box off


% --- broadness ----
ax3 = axes('position',[0.08 0.65 0.8 0.035])
for cc = 1:CC
    tmp = find(CL==cc);
    dat2 = (dX(tmp)+dY(tmp))./2;
    SUPER_COLL(tmp,5) = dat2;
    %     scatter(repmat(cc,1,length(dat2)),dat2,10,clrs(cc,:))
    plotSpread(dat2,'distributionIdx',repmat(cc,length(dat2),1),'distributionMarkers',{'.'},'distributionColors',{clrs(cc,:)})
    hold on
    plot([cc-0.3 cc+0.3],[median(dat2) median(dat2)],'-','color',clrs(cc,:),'linewidth',1)
    CLUSTER_INFO(cc,3,1)=mean(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,3,2)=median(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,3,3)=min(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,3,4)=max(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,3,5)=std(dat2);%mean,median,min max, std
end
set(gca,'xtick',1:10,'xticklabel',1:10,'fontname','Arial','fontsize',fsize2)
axis([0.5 CC+0.5 0 4])
box off
ylabel('a.u.','fontname','Arial','fontsize',fsize2)
axes('position',[0.08 0.695 0.8 0.005])
text(0,0,'width of spatio-temporal response','fontname','Arial','fontsize',fsize2)
axis off
box off


% --- temporal ----
ax2 = axes('position',[0.08 0.715 0.8 0.035])

for cc = 1:CC
    tmp = find(CL==cc);
    dat = 5-cY(tmp);
    dat2 = [];
    for d = 1:length(dat)
        now = dat(d);
        in = floor(now);
        rest = now-in;
        if rest >0
            slp = (FR(in+1)-FR(in));
            um = slp*rest+FR(in);
            dat2 = [dat2;um];
        else
            dat2 = [dat2;FR(in)];
        end
    end
    SUPER_COLL(tmp,4) = dat2;
    plotSpread(dat2,'distributionIdx',repmat(cc,length(dat2),1),'distributionMarkers',{'.'},'distributionColors',{clrs(cc,:)})
    %     scatter(repmat(cc,1,length(dat2)),dat2,10,clrs(cc,:))
    hold on
    plot([cc-0.3 cc+0.3],[median(dat2) median(dat2)],'-','color',clrs(cc,:),'linewidth',1)
    CLUSTER_INFO(cc,4,1)=mean(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,4,2)=median(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,4,3)=min(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,4,4)=max(dat2);%mean,median,min max, std
    CLUSTER_INFO(cc,4,5)=std(dat2);%mean,median,min max, std
end
set(gca,'xtick',1:10,'xticklabel',1:10,'fontname','Arial','fontsize',fsize2)
axis([0.5 CC+0.5 0 8])
box off
ylabel('Hz','fontname','Arial','fontsize',fsize2)
axes('position',[0.08 0.76 0.8 0.005])
text(0,0,'temporal frequency preference','fontname','Arial','fontsize',fsize2)
axis off
box off

% --- contrast ----
ax4 = axes('position',[0.08 0.585 0.8 0.035])
for cc = 1:CC
    tmp = find(CL==cc);
    dat = togoDG(tmp);
    dat2 = []; tk = [];
    for d = 1:length(tmp)
        now = dat(d);
        tmp2 = find(togoChrp==now);
        if ~isempty(tmp2)
            tk = [tk;d];
            dat2 = [dat2;tmp2];
        end
    end
    SUPER_COLL(tmp(tk),6) = cntrstIDX(dat2);
    if ~isempty(dat2)
        %         scatter(repmat(cc,1,length(dat2)),cntrstIDX(dat2),10,clrs(cc,:))
        plotSpread(cntrstIDX(dat2),'distributionIdx',repmat(cc,length(dat2),1),'distributionMarkers',{'.'},'distributionColors',{clrs(cc,:)})
        
        hold on
        plot([cc-0.3 cc+0.3],[median(cntrstIDX(dat2)) median(cntrstIDX(dat2))],'-','color',clrs(cc,:),'linewidth',1)
        CLUSTER_INFO(cc,5,1)=mean(cntrstIDX(dat2));%mean,median,min max, std
        CLUSTER_INFO(cc,5,2)=median(cntrstIDX(dat2));%mean,median,min max, std
        CLUSTER_INFO(cc,5,3)=min(cntrstIDX(dat2));%mean,median,min max, std
        CLUSTER_INFO(cc,5,4)=max(cntrstIDX(dat2));%mean,median,min max, std
        CLUSTER_INFO(cc,5,5)=std(cntrstIDX(dat2));%mean,median,min max, std
    end
end
set(gca,'xtick',1:10,'xticklabel',1:10,'fontname','Arial','fontsize',fsize2)
axis([0.5 CC+0.5 -1 1])
box off
ylabel('a.u.','fontname','Arial','fontsize',fsize2)
axes('position',[0.08 0.63 0.8 0.005])
text(0,0,'contrast preference index','fontname','Arial','fontsize',fsize2)
axis off
box off



% ----ON-OFF----
ax5 = axes('position',[0.08 0.52 0.8 0.035])
for cc = 1:CC
    tmpo = find(CL==cc);
    ids = togoDG(tmpo);
    flsh = cell2mat(info(ids,8));
    SUPER_COLL(tmpo,7) = flsh;
    if ~isempty(flsh)
        tmp = find(flsh==-1);
        CLUSTER_INFO(cc,6,1)=length(tmp);%mean,median,min max, std
        rectangle('position',[cc-0.3 0 0.2 100./length(tmpo).*length(tmp)],'edgecolor',clrs(cc,:),'facecolor',clrs(cc,:))
        hold on
        tmp = find(flsh==1);
        CLUSTER_INFO(cc,6,2)=length(tmp);%mean,median,min max, std
        rectangle('position',[cc-0.1 0 0.2 100./length(tmpo).*length(tmp)],'edgecolor',clrs(cc,:),'facecolor','w')
        tmp = find(flsh==2);
        CLUSTER_INFO(cc,6,3)=length(tmp);%mean,median,min max, std
        rectangle('position',[cc+0.1 0 0.2 100./length(tmpo).*length(tmp)],'edgecolor',clrs(cc,:),'facecolor',[clrs(cc,:) 0.3])
        tmp = find(flsh==0);
        %       rectangle('position',[cc+0.075 0 0.075 100./length(tmpo).*length(tmp)],'edgecolor','k','facecolor','w')
        
    end
    
end
rectangle('position',[7.5 40 0.2 10],'edgecolor','k','facecolor','k')
text(7.75,45,'OFF','fontname','Arial','fontsize',fsize1)
rectangle('position',[8.5 40 0.2 10],'edgecolor','k','facecolor','w')
text(8.75,45,'ON','fontname','Arial','fontsize',fsize1)
rectangle('position',[9.5 40 0.2 10],'edgecolor','k','facecolor',[0 0 0 0.3])
text(9.75,45,'ON-OFF','fontname','Arial','fontsize',fsize1)
set(gca,'xtick',1:10,'xticklabel',1:10,'fontname','Arial','fontsize',fsize2)
axis([0.5 CC+0.5 0 52])
box off
ylabel('%','fontname','Arial','fontsize',fsize2)
axes('position',[0.08 0.565 0.8 0.005])
text(0,0,'polarity (OFF, ON, ON-OFF)','fontname','Arial','fontsize',fsize2)
axis off
box off

%----transiency----
ax6 = axes('position',[0.08 0.45 0.8 0.035])
load(fullfile(path1,'h_flashSpikes'))
for cc = 1:CC
    tmpo = find(CL==cc);
    ids = togoDG(tmpo);
    flsh = cell2mat(info(ids,8));
    trans = [];
    if ~isempty(flsh)
        for g=1:length(ids)
            if flsh(g)~=0
                
                numU = info{ids(g),2};
                dateU = info{ids(g),1};
                currUnit=numU;
                
                data=flash_raster(ids(g),:);
                
                allrate=[];
                for f=1:length(data)
                    spikes=data{f};
                    rate=MEA_spikerates(spikes,40,15000);
                    allrate=[allrate;rate];
                end
                meanrate=mean(allrate,1);
                bck = mean(meanrate(50:1500));
                stback = std( meanrate(50:1500));
                [m1,m2] = max(meanrate);
                tmp = find(meanrate(m2:end)<=bck+2*stback); tmp = tmp(1);
                %                 tmp = tmp/1000*40;
                if tmp>1000
                    tmp = 1000;
                end
                trans = [trans;tmp];
                SUPER_COLL(tmpo(g),8) = tmp;
            end
        end
    end
    plotSpread(trans,'distributionIdx',repmat(cc,length(trans),1),'distributionMarkers',{'.'},'distributionColors',{clrs(cc,:)})
    
    %     scatter(repmat(cc,1,length(trans)),trans,10,clrs(cc,:))
    hold on
    plot([cc-0.3 cc+0.3],[median(trans) median(trans)],'-','color',clrs(cc,:),'linewidth',1)
    
    CLUSTER_INFO(cc,7,1)=mean(trans);%mean,median,min max, std
    CLUSTER_INFO(cc,7,2)=median(trans);%mean,median,min max, std
    CLUSTER_INFO(cc,7,3)=min(trans);%mean,median,min max, std
    CLUSTER_INFO(cc,7,4)=max(trans);%mean,median,min max, std
    CLUSTER_INFO(cc,7,5)=std(trans);%mean,median,min max, std
end
set(gca,'xtick',1:10,'xticklabel',1:10,'fontname','Arial','fontsize',fsize2)
axis([0.5 CC+0.5 0 1020])
box off
ylabel('ms','fontname','Arial','fontsize',fsize2)
axes('position',[0.08 0.495 0.8 0.005])
text(0,0,'transiency of flash response','fontname','Arial','fontsize',fsize2)
axis off
box off
%% label example cells
% ex_names = {'A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','B6','B7','C1','C2','C3'};
ex_names = {'A_{1}','A_{2}','A_{3}','A_{4}','A_{5}',...
    'B_{1}','B_{2}','B_{3}','B_{4}','B_{5}','B_{6}','B_{7}','C_{1}','C_{2}','C_{3}'};

counter = zeros(15,1);
for ii = 1:length(COLL_ID)
    inow = COLL_ID(ii);
    pos = find(SUPER_COLL(:,1)==inow);
    cl = SUPER_COLL(pos,2);
    counter(cl) = counter(cl)+1;
    for aa = [1:4 6]
        eval(['currAX=ax',int2str(aa),';'])
        if cl == 6
            if counter(cl) == 1
                plot(currAX,cl-0.2,SUPER_COLL(pos,aa+2),'ok','markerfacecolor','k','markersize',4);
            elseif counter(cl) == 2
                plot(currAX,cl-0.1,SUPER_COLL(pos,aa+2),'^k','markerfacecolor','k','markersize',4);
            elseif counter(cl) == 3
                plot(currAX,cl,SUPER_COLL(pos,aa+2),'sk','markerfacecolor','k','markersize',4);
            elseif counter(cl) == 4
                plot(currAX,cl+0.1,SUPER_COLL(pos,aa+2),'o','markerfacecolor','none','markeredgecolor','k','markersize',4,'linewidth',1);
            elseif counter(cl) == 5
                plot(currAX,cl+0.2,SUPER_COLL(pos,aa+2),'^','markerfacecolor','none','markeredgecolor','k','markersize',4,'linewidth',1);
            end
        elseif  counter(cl) == 1
            plot(currAX,cl,SUPER_COLL(pos,aa+2),'ok','markerfacecolor','k','markersize',4);
        elseif counter(cl) == 2
            plot(currAX,cl,SUPER_COLL(pos,aa+2),'ok','markerfacecolor','none','markersize',4);
        end
    end
    %     text(ax1,cl+0.1,SUPER_COLL(pos,3),)
end
save(fullfile(path1,'EXAMPLE_INFO'),'SUPER_COLL','CLUSTER_INFO')
%% label example cell names
axes('position',[0.08 0.84 0.8 0.04])
list_label = {'o','^','s','o','^'};
list_color = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0];
list_size = [4 4 4 4 4];
counter = zeros(15,1);
for ii = 1:length(COLL_ID)
    inow = COLL_ID(ii);
    pos = find(SUPER_COLL(:,1)==inow);
    cl = SUPER_COLL(pos,2);
    counter(cl) = counter(cl)+1;
    
    if cl == 6
        if counter<4
        plot(cl,0.95-(0.25*counter(cl)-1),list_label{counter(cl)},'markeredgecolor',list_color(counter(cl),:),...
            'markerfacecolor','k','markersize',list_size(counter(cl)))
        else
            plot(cl,0.95-(0.25*counter(cl)-1),list_label{counter(cl)},'markeredgecolor',list_color(counter(cl),:),...
            'markerfacecolor','none','markersize',list_size(counter(cl)))
        end
        hold on
        text(cl+0.1,0.95-(0.25*counter(cl)-1),ex_names{ii},'fontname','Arial','fontsize',fsize1)
    elseif counter(cl) == 1
        plot(cl,0.95-(0.25*counter(cl)-1),'ok','markerfacecolor','k','markersize',4)
        hold on
        text(cl+0.1,0.95-(0.25*counter(cl)-1),ex_names{ii},'fontname','Arial','fontsize',fsize1)
    elseif counter(cl) == 2
        plot(cl,0.95-(0.25*counter(cl)-1),'ok','markerfacecolor','none','markersize',4)
        hold on
        text(cl+0.1,0.95-(0.25*counter(cl)-1),ex_names{ii},'fontname','Arial','fontsize',fsize1)
    end
end
axis([0.5 CC+0.5 -inf inf])
axis off
box off
%% letters
axes('position',[ 0.01 0.99 0.01 0.01])
text(0,0,'A','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off


axes('position',[ 0.01 0.828 0.01 0.01])
text(0,0,'B','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.01 0.76 0.01 0.01])
text(0,0,'C','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.01 0.695 0.01 0.01])
text(0,0,'D','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.01 0.63 0.01 0.01])
text(0,0,'E','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.01 0.565 0.01 0.01])
text(0,0,'F','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.01 0.497 0.01 0.01])
text(0,0,'G','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.01 0.42 0.01 0.01])
text(0,0,'H','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.24 0.42 0.01 0.01])
text(0,0,'I','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.45 0.42 0.01 0.01])
text(0,0,'J','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off

axes('position',[ 0.7 0.42 0.01 0.01])
text(0,0,'K','fontsize',fsize3,'fontweight','bold','fontname','Arial')
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
box off
%% find example cells
% clc
% midget = [10 6 4];%B5, B1, A4
% parasol = [1 3];
% others = [2 5 7 8 9 11:15];
% disp('----midget-----')
% for ii = 1:length(midget)
%     inow = midget(ii);
%     inow2 = COLL_ID(inow);
%     iddg = find(togoDG==inow2);
%     idch = find(togoChrp==inow2);
%     [~,mx] = max(normP(iddg,:));
%     disp([int2str(inow),' // ','contrastIDX: ',num2str(cntrstIDX(idch)),' | ',int2str(spat(mx)),'um, ',int2str(freq(mx)),'Hz',...
%         ' | ',num2str(7-cX(iddg)),'um, ',num2str(5-cY(iddg)),'Hz, | clust: ',int2str(CL(iddg))])
% end
% disp('----parasol-----')
% for ii = 1:length(parasol)
%     inow = parasol(ii);
%     inow2 = COLL_ID(inow);
%     iddg = find(togoDG==inow2);
%     idch = find(togoChrp==inow2);
%     [~,mx] = max(normP(iddg,:));
%     disp([int2str(inow),' // ','contrastIDX: ',num2str(cntrstIDX(idch)),' | ',int2str(spat(mx)),'um, ',int2str(freq(mx)),'Hz',...
%         ' | ',num2str(7-cX(iddg)),'um, ',num2str(5-cY(iddg)),'Hz, | clust: ',int2str(CL(iddg))])
% end
%
% disp('----others-----')
% for ii = 1:length(others)
%     inow = others(ii);
%     inow2 = COLL_ID(inow);
%     iddg = find(togoDG==inow2);
%     idch = find(togoChrp==inow2);
%     [~,mx] = max(normP(iddg,:));
%     disp([int2str(inow),' // ','contrastIDX: ',num2str(cntrstIDX(idch)),' | ',int2str(spat(mx)),'um, ',int2str(freq(mx)),'Hz',...
%         ' | ',num2str(7-cX(iddg)),'um, ',num2str(5-cY(iddg)),'Hz, | clust: ',int2str(CL(iddg))])
% end
%%
for ii = 1:length(COLL_ID)
    inow = COLL_ID(ii);
    iddg = find(togoDG==inow);
    cln = CL(iddg);
    disp(['example ',int2str(ii),': cluster ',int2str(cln)])
end
%%
disp('number of cells per cluster')
for c = 1:CC
    tmp = find(CL ==c);
    disp(['cluster ',int2str(c),': ',int2str(length(tmp)),' cells'])
end
%% save
exportgraphics(gcf,fullfile(pathSave,'NEW_Fig5_types.png'),'Resolution',600)
exportgraphics(gcf,fullfile(pathSave,'NEW_Fig5_types.eps'),'Resolution',600)

% saveas(gcf,fullfile(pathSave,'NEW_Fig5_types.tif'))