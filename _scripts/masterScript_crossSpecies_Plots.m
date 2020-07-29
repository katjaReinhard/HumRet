%% speed preference from DG
clear
close all
superpath='e:\all_species_temp';
cols='bgr';
Dtype='wph';
spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
speed=spat.*freq/1000; speed=speed';
% uspeed=unique(speed);
% uspeed=[0.5 1 2 4 8 16 32];
uspeed=[1 2 4 8 16];
% uspeed=[1 2 4 8];

ids=[];
for u=1:length(uspeed)
    tmp=find(speed==uspeed(u));
    ids=[ids;tmp];
end

speed=speed(ids);
figure
colldata=cell(length(uspeed),3);
for d=1:3
    
    load(fullfile(superpath,'newAttempt','DG',[Dtype(d),'_FFTpeaks']))
    data=reshape(FFTpeaks(:,:,2),size(FFTpeaks,1),24);
    data=data./repmat(max(data,[],2),1,24);
    data=data(:,ids);
    
    
    
    
    newdata=[];
    for s=1:length(uspeed)
        tmp=find(speed==uspeed(s));
        newdata=[newdata max(data(:,tmp),[],2)];
    end
    
    
    DGspeed=[];
    for c=1:size(newdata,1)
        su=sum(newdata(c,:));
        cumdata=cumsum(newdata(c,:));
        tmp1=find(cumdata<su/2);
        if isempty(tmp1)
            optDG=uspeed(1);
        else
            tmp1=tmp1(end);
            tmp2=find(cumdata>su/2);
            if isempty(tmp2)
                optDG=uspeed(end);
            else
                tmp2=tmp2(1);
                slope=(cumdata(tmp2)-cumdata(tmp1))/(uspeed(tmp2)-uspeed(tmp1));
                optDG=(su/2-cumdata(tmp1))/slope+uspeed(tmp1);
            end
        end
        DGspeed=[DGspeed;optDG];
    end
    if d==1
        mall=DGspeed;
    elseif d==2
        pall=DGspeed;
    else
        hall=DGspeed;
    end
    DGspeed=[DGspeed;-1;max(uspeed)+1];
    figure(1)
    subplot(3,2,d*2-1)
    hist(DGspeed,38)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cols(d),'EdgeColor',cols(d))
    med=median(DGspeed(1:end-2));
    st=std(DGspeed(1:end-2));
    conf = 1.96*st/sqrt(length(DGspeed)-2);
    title([Dtype(d), ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(DGspeed)-2),')'])
    hold on
    plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
    plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
    yy=get(gca,'ylim');
    plot([med-conf med-conf],[0 yy(2)+1],'--k')
    plot([med+conf med+conf],[0 yy(2)+1],'--k')
    %     set(gca,'xtick',speeds)
    % set(gca,'xticklabel',speeds)
    xlabel('speed (mm/s)')
    axis([0 uspeed(end) 0 yy(2)+1])
    
    subplot(3,2,[2 4])
    c=cdfplot(DGspeed(1:end-2))
    set(c,'color',cols(d),'linestyle','-','linewidth',2)
    hold on
    xlabel('speed (mm/s)')
    ylabel('probability')
    title('DG speed tuning')
    
    figure(2)
    collme=[]; collst=[]; collci=[];
    for s=1:length(uspeed)
        tmp=find(speed==uspeed(s));
        %         tmp=find(speed(13:16)==uspeed(s)); tmp=tmp+12;
        currdata=max(data(:,tmp),[],2);
        colldata{s,d}=currdata;
        %         currdata=mean(data(:,tmp),2);
        me=mean(currdata);
        st=std(currdata);
        ci=1.96*st/sqrt(size(currdata,1));
        collme=[collme;me];
        collst=[collst;st];
        collci=[collci;ci];
        subplot(2,2,d)
        plot([uspeed(s) uspeed(s)],[me-st me+st],'-','color',cols(d))
        hold on
        plot([uspeed(s) uspeed(s)],[me-ci me+ci],'-','color',cols(d),'linewidth',4)
        subplot(2,2,4)
        plot([uspeed(s) uspeed(s)],[me-st me+st],'-','color',cols(d))
        hold on
        plot([uspeed(s) uspeed(s)],[me-ci me+ci],'-','color',cols(d),'linewidth',4)
    end
    
    subplot(2,2,d)
    plot(uspeed,collme,'o-','color',cols(d),'markerfacecolor',cols(d),'markeredgecolor',cols(d),'linewidth',2)
    subplot(2,2,4)
    plot(uspeed,collme,'o-','color',cols(d),'markerfacecolor',cols(d),'markeredgecolor',cols(d),'linewidth',2)
end
p1=ranksum(mall, pall);
p2=ranksum(mall, hall);
p3=ranksum(pall, hall);

pvals=zeros(3,length(uspeed));
cnt=0;
for c=1:2
    for cc=c+1:3
        cnt=cnt+1;
        for s=1:length(uspeed)
            data1=colldata{s,c};
            data2=colldata{s,cc};
            p=ranksum(data1,data2);
            pvals(cnt,s)=p;
        end
    end
end
figure(2)
subplot(2,2,4)
for c=1:3
    for s=1:length(uspeed)
        text(uspeed(s),1.2+c*0.1,num2str(round(pvals(c,s)*100)/100))
    end
end
axis([0 uspeed(end)+1 0 1.6])
%% speed preference 6vel for max(black,white)
clear
close all
superpath='f:\all_species_temp';
savepath=fullfile(superpath,'overview_pics');


names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[0.2; 1; 2; 4; 8; 16];
types='wph';
cols='bgr';
titles={'mouse','pig','human'};

for t=1:length(types)
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    allval=zeros(size(info,1),2);
    black=cell2mat(info(goodB,24));
    white=cell2mat(info(goodW,25));
    allval(goodB,1)=black;
    allval(goodW,2)=white;
    allval=max(allval,[],2);
    tmp=find(allval>0);
    allval=allval(tmp);
    if t==1
        mall=allval;
    elseif t==2
        pall=allval;
    else
        hall=allval;
    end
    
    allval=[allval;-1;17];
    
    
    
    subplot(3,2,t*2-1)
    hist(allval,38)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
    med=median(allval(1:end-2));
    st=std(allval(1:end-2));
    conf = 1.96*st/sqrt(length(allval)-2);
    title([titles{t}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(allval)-2),')'])
    hold on
    plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
    plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
    yy=get(gca,'ylim');
    plot([med-conf med-conf],[0 yy(2)+1],'--k')
    plot([med+conf med+conf],[0 yy(2)+1],'--k')
    %     set(gca,'xtick',speeds)
    % set(gca,'xticklabel',speeds)
    xlabel('speed (mm/s)')
    axis([0 16 0 yy(2)+1])
    
    subplot(3,2,[2 4])
    c=cdfplot(allval(1:end-2))
    set(c,'color',cols(t),'linestyle','-','linewidth',2)
    hold on
    xlabel('speed (mm/s)')
    ylabel('probability')
    title('better bar')
    
    
end

set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(savepath,'speed_all_noSlow.bmp'))
% saveas(gcf,fullfile(savepath,'speed_all_noSlow.fig'))
% close

p1=ranksum(mall, pall);
p2=ranksum(mall, hall);
p3=ranksum(pall, hall);
%% latencies and polarity
% calculated by finding first "real" peak and then 75% of its hight
% manual correction for random burst induced or very unstable peaks
% based on sigma=60 plots

% for pig, exclude 20140723d

close all
clear

badIDs=202:250;

path1='F:\all_species_temp';
list='wph';
cols='bgr';
allLat=[]; type=[]; flashID=[];
for l=1:length(list);
    load(fullfile(path1,[list(l),'_latency2']))
    load(fullfile(path1,[list(l),'_info']))
    
    if l==2
        goodies=[];
        for f=1:length(flash)
            if isempty(find(badIDs==flash(f)))
                goodies=[goodies;f];
            end
        end
        flash=flash(goodies);
        lats=lats(goodies,:);
        peaks=peaks(goodies,:);
    end
    
    
    for ii=1:length(info)
        if ~isempty(find(flash==ii))
            info{ii,26}=1;
        else
            info{ii,26}=0;
        end
    end
    save(fullfile(path1,[list(l),'_info']),'info','goodones')
    allLat=[allLat;lats];
    type=[type;zeros(length(lats),1)+l];
    flashID=[flashID;flash];
end

save(fullfile('F:\all_species_temp','latency'),'allLat','type','flashID')

path2='F:\all_species_temp\final_pics';
figure
for i=1:4
    tmp=find(allLat(:,i)>0);
    data=allLat(tmp,i);
    me=median(data);
    st=std(data);
    data=[0;data;2000];
    
    subplot(2,2,i)
    [h int]=hist(data,40);
    bar(int(2:end-1),100/length(data)*h(2:end-1))
    hold on
    plot(me,max(100/length(data)*h(2:end-1))+2,'or')
    plot([me-st me+st],[max(100/length(data)*h(2:end-1))+2 max(100/length(data)*h(2:end-1))+2],'-r')
    axis([-100 2000 0 max(100/length(data)*h(2:end-1))+4])
    title(['latency all species for step ',int2str(i),'; median=',num2str(round(me))])
    ylabel('% of cells with flash-responses')
    xlabel('latency [ms]')
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'latency_allSpecies_perCond.emf'))
close

close all
figure(1)
compl=reshape(allLat,length(allLat)*4,1);
tmp=find(compl>0);
data=compl(tmp);
me=median(data);
st=std(data);
data=[0;data;2000];
subplot(2,2,1)
figure
h=histfit(data,40,'kernel');
figure(1)
% [h int]=hist(data,40);
int=get(h(1),'xdata');
int=mean(int(2:3,:));
hi=get(h(1),'ydata');
hi=hi(2,:);
xfit=get(h(2),'xdata');
yfit=get(h(2),'ydata');
bar(int(2:end-1),100/length(data-2)*hi(2:end-1))
hold on
plot(xfit,100/length(data-2)*yfit,'-')
plot(me,max(100/length(data-2)*hi(2:end-1))+2,'or')
plot([me-st me+st],[max(100/length(data-2)*hi(2:end-1))+2 max(100/length(data-2)*hi(2:end-1))+2],'-r')
axis([-100 2000 0 max(100/length(data-2)*hi(2:end-1))+4])
title(['all','; median=',num2str(round(me))])
ylabel('% of cells with flash-responses')
xlabel('latency [ms]')

for t=1:length(list)
    tst=find(type==t);
    compl=reshape(allLat(tst,:),length(tst)*4,1);
    tmp=find(compl>0);
    data=compl(tmp);
    me=median(data);
    st=std(data);
    
    data=[0;data;2000];
    figure
    h=histfit(data,40,'kernel');
    figure(1)
    
    subplot(2,2,t+1)
    
    int=get(h(1),'xdata');
    keepint=reshape(int(1:2,:),length(int)*2,1);
    keepint=sort(keepint);
    keepint=[keepint;keepint(end)+diff(keepint(end-1:end));keepint(end)+diff(keepint(end-1:end))];
    int=mean(int(2:3,:));
    hi=get(h(1),'ydata');
    keephi=[];
    for x=1:length(hi)
        keephi=[keephi;hi(2:3,x)];
    end
    keephi(1:2)=0;
    keephi(end-1:end)=0;
    keephi=[0;keephi;0];
    keephi=100/length(data-2)*keephi;
    hi=hi(2,:);
    xfit=get(h(2),'xdata');
    yfit=get(h(2),'ydata');
    bar(int(2:end-1),100/length(data-2)*hi(2:end-1))
    hold on
    plot(xfit,100/length(data-2)*yfit,'-')
    plot(me,max(100/length(data-2)*hi(2:end-1))+2,'or')
    plot([me-st me+st],[max(100/length(data-2)*hi(2:end-1))+2 max(100/length(data-2)*hi(2:end-1))+2],'-r')
    axis([-100 2000 0 max(max(100/length(data-2)*yfit),max(100/length(data-2)*hi(2:end-1)))+4])
    title([list(t),' (all conditions)','; median=',num2str(round(me))])
    ylabel('% of cells with flash-responses')
    xlabel('latency [ms]')
    figure(100)
    plot(keepint,keephi,cols(t))
    hold on
end
figure(1)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'latency_acrossCond.emf'))
close
figure(100)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'latency_acrossCond2.emf'))
close

for t=1:3
    figure
    for i=1:4
        tmp=find(type==t);
        la=allLat(tmp,:);
        tmp=find(la(:,i)>0);
        data=la(tmp,i);
        me=median(data);
        st=std(data);
        data=[0;data;2000];
        subplot(2,2,i)
        [h int]=hist(data,50);
        bar(int(2:end-1),100/length(data)*h(2:end-1))
        hold on
        plot(me,max(100/length(data)*h(2:end-1))+2,'or')
        plot([me-st me+st],[max(100/length(data)*h(2:end-1))+2 max(100/length(data)*h(2:end-1))+2],'-r')
        axis([-100 2000 0 max(100/length(data)*h(2:end-1))+4])
        title([list(t),', condition ',int2str(i),'; median=',num2str(round(me))])
        ylabel('% of cells with flash-responses')
        xlabel('latency [ms]')
    end
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(path2,['latency_',list(t),'_perCond.emf']))
    close
end

early=zeros(size(allLat,1),size(allLat,2));
all=zeros(size(allLat,1),size(allLat,2));
for i=1:4
    tmp=find(allLat(:,i)<=550);
    tmp2=find(allLat(tmp,i)>0);
    early(tmp(tmp2),i)=1;
    tmp=find(allLat(:,i)>0);
    all(tmp,i)=1;
end

onoff=zeros(size(allLat,1),1);
anyresp=zeros(size(allLat,1),1);
for i=1:length(early)
    if isempty(find(early(i,:)>0));
        onoff(i)=0;
    elseif ~isempty(find(early(i,[1 4])==1)) && isempty(find(early(i,[2 3])==1))
        onoff(i)=-1;
    elseif isempty(find(early(i,[1 4])==1)) && ~isempty(find(early(i,[2 3])==1))
        onoff(i)=1;
    else
        onoff(i)=2;
    end
    
    if isempty(find(all(i,:)>0));
        anyresp(i)=0;
    elseif ~isempty(find(all(i,[1 4])==1)) && isempty(find(all(i,[2 3])==1))
        anyresp(i)=-1;
    elseif isempty(find(all(i,[1 4])==1)) && ~isempty(find(all(i,[2 3])==1))
        anyresp(i)=1;
    else
        anyresp(i)=2;
    end
end

save(fullfile('F:\all_species_temp','polarity'),'type','onoff','anyresp','flashID')

figure
re=[-1 1 2 0];
for t=1:length(list)
    tst=find(type==t);
    num=[];
    for r=1:length(re)
        tmp=find(onoff(tst)==re(r));
        num=[num length(tmp)];
    end
    subplot(2,2,t)
    h=pie(num)
    hText = findobj(h,'Type','text'); % text handles
    percentValues = get(hText,'String'); % percent values
    str = {'off ';'on ';'onoff ';'none '}; % text
    tmp=find(num>0);
    combinedstrings = strcat(str(tmp),percentValues);
    set(hText,{'String'},combinedstrings);
    title(['polartiy based on early responses (<=550ms) of ',list(t),' (total: ',int2str(length(tst)),' cells)'])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'polarity_early.emf'))
close

figure
re=[-1 1 2];
for t=1:length(list)
    tst=find(type==t);
    num=[];
    for r=1:length(re)
        tmp=find(anyresp(tst)==re(r));
        num=[num length(tmp)];
    end
    subplot(2,2,t)
    h=pie(num)
    hText = findobj(h,'Type','text'); % text handles
    percentValues = get(hText,'String'); % percent values
    str = {'off ';'on ';'onoff ';'none '}; % text
    tmp=find(num>0);
    combinedstrings = strcat(str(tmp),percentValues);
    set(hText,{'String'},combinedstrings);
    title(['polartiy based on all responses of ',list(t),' (total: ',int2str(length(tst)),' cells)'])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'polarity_all.emf'))
close

%% chirp temporal tuning: p-values
clear
close all
nfft=4;
nlist=[2000 4000 6000 8000 10000 15000];
superpath='f:\all_species_temp';
load(fullfile(superpath,'newAttempt2','chirp',['chirpRatios_',int2str(nlist(nfft))]))

step=5;
data=colldata2;
fr=frs2{1};
dataset='smooth'; %data2
% dataset='noSmooth_noNorm'; %data1a
% dataset='noSmooth_Norm'; %data1b

pvals=zeros(ceil(size(data{1},2)/step),3); %M-P, M-H, P-H
cnt=0;
for t=1:step:size(data{1},2)
    cnt=cnt+1;
    now=t:t+step-1;
    if t+step-1 > size(data{1},2)
        now=t:size(data{1},2);
    end
    p=ranksum(mean(data{1}(:,now),2),mean(data{2}(:,now),2));
    pvals(cnt,1)=p;
    p=ranksum(mean(data{1}(:,now),2),mean(data{3}(:,now),2));
    pvals(cnt,2)=p;
    p=ranksum(mean(data{2}(:,now),2),mean(data{3}(:,now),2));
    pvals(cnt,3)=p;
end

figure
startpoint=floor(step/2);
if startpoint==0
    startpoint=1;
end
if length(fr(1,startpoint:step:end))<size(pvals,1)
    frnow=fr(1,startpoint:step:end);
    frnow=[frnow frnow(end)+diff([frnow(1) frnow(2)])];
else
    frnow=fr(1,startpoint:step:end);
end
subplot(2,2,1)
% patch([fr(1,startpoint:step:end)' fliplr(fr(1,startpoint:step:end)')],[pvals(:,1) repmat(0,length(fr(1,startpoint:step:end)'),1)],[0 0.5 0.5],'edgecolor',[0 0.5 0.5])
plot(frnow,pvals(:,1),'color',[0 0.5 0.5],'linewidth',3)
title('mouse - pig')
subplot(2,2,2)
% patch([fr(1,startpoint:step:end)' fliplr(fr(1,startpoint:step:end)')],[pvals(:,2) repmat(0,length(fr(1,startpoint:step:end)'),1)],[0.5 0 0.5],'edgecolor',[0.5 0 0.5])
plot(frnow,pvals(:,2),'color',[0.5 0 0.5],'linewidth',3)
title('mouse - human')
subplot(2,2,3)
% patch([fr(1,startpoint:step:end)' fliplr(fr(1,startpoint:step:end)')],[repmat(-0.001,length(fr(1,startpoint:step:end)'),1) fliplr(pvals(:,3)) ],[0.5 0.5 0],'edgecolor',[0.5 0.5 0])
plot(frnow,pvals(:,3),'color',[0.5 0.5 0],'linewidth',3)
title('pig - human')
for s=1:3
    subplot(2,2,s)
    hold on
    %     plot([fr(1) fr(end)],[0.05 0.05],'-k')
    plot([fr(1) fr(end)],[0.01 0.01],'-k')
    axis([fr(1) fr(end) 0 0.05])
    xlabel('fr [Hz]')
    ylabel('p-value')
end
% saveas(gcf,fullfile(superpath,'newAttempt','stat_figs',['chirpPval_',dataset,'_',int2str(step),'step',int2str(nlist(nfft)),'.bmp']))
% saveas(gcf,fullfile(superpath,'newAttempt','stat_figs',['chirpPval_',dataset,'_',int2str(step),'step',int2str(nlist(nfft)),'.fig']))
% close
% save(fullfile(superpath,'newAttempt','chirp',[Dtype,'_chirpFFT']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates','varInfo')
%% chirp temporal tuning (binary, NFFT=8000)
clear
close all

nfft=4; %2000,4000,6000,8000,10000,15000
nlist=[2000 4000 6000 8000 10000 15000];
superpath='f:\all_species_temp';

Dtype='wph';
cols='bgr';

colldata2=cell(3,1); colldata1a=cell(3,1); colldata1b=cell(3,1);colldata3=cell(3,1);
frs1=cell(3,1); frs2=cell(3,1); togos=cell(3,1);
figure(1)
for d=1:3
    load(fullfile(superpath,'newAttempt2','chirp',[Dtype(d),'_chirpFFT']))
    
    resp=cell2mat(reshape(chirpFFTabs(:,6,nfft),size(chirpFFTabs,1),1)')'; %non-smooth binary
    realgood=[];
    for ts=1:size(chirpFFTabs,1)
        if ~isempty(chirpFFTabs{ts,1,1})
            realgood=[realgood;ts];
        end
    end
    togo=togo(ts);
    stim=cell2mat(reshape(chirpFFTabs(:,7,nfft),size(chirpFFTabs,1),1)')';
    data2=cell2mat(reshape(chirpFFTsmooth(:,6,nfft),size(chirpFFTabs,1),1)')'; %smooth binary
    fr1=cell2mat(reshape(chirpFFTabs(:,1,nfft),size(chirpFFTabs,1),1));
    fr2=cell2mat(reshape(chirpFFTsmooth(:,1,nfft),size(chirpFFTabs,1),1)')';
    frs1{d}=fr1;
    frs2{d}=fr2;
    
    clear chirpFFTabs chirpFFTsmooth
    
    data1a=resp./stim;
    resp2=[]; stim2=[];
    for i=1:size(resp,1)
        resp2=[resp2;resp(i,:)/max(resp(i,:))];
        stim2=[stim2;stim(i,:)/max(stim(i,:))];
    end
    data1b=resp2./stim2;
    
    data3=[];
    for i=1:size(data2,1)
        data3=[data3;data2(i,:)/max(data2(i,:))];
    end
    
    subplot(4,4,d)
    plot(fr1(1,:),data1a,'color',[0.6 0.6 0.6])
    hold on
    subplot(4,4,d+4)
    plot(fr1(1,:),data1b,'color',[0.6 0.6 0.6])
    hold on
    subplot(4,4,d+8)
    plot(fr2(1,:),data2,'color',[0.6 0.6 0.6])
    hold on
    subplot(4,4,d+12)
    plot(fr2(1,:),data3,'color',[0.6 0.6 0.6])
    hold on
    
    subplot(4,4,d)
    plot(fr1,median(data1a),'color',cols(d),'linewidth',2)
    patch([fr1 fliplr(fr1)],[median(data1a)-std(data1a) fliplr(median(data1a)+std(data1a))],cols(d),'facealpha',0.4)
    %     plot(fr1,median(data1a)+std(data1a),'color',cols(d),'linewidth',1)
    %     plot(fr1,median(data1a)-std(data1a),'color',cols(d),'linewidth',1)
    title('non-smoothed FFTresp/FFTstim')
    axis tight
    subplot(4,4,4)
    plot(fr1,median(data1a),'color',cols(d),'linewidth',2)
    hold on
    patch([fr1 fliplr(fr1)],[median(data1a)-std(data1a) fliplr(median(data1a)+std(data1a))],cols(d),'facealpha',0.2)
    %     plot(fr1,median(data1a)+std(data1a),'color',cols(d),'linewidth',1)
    %     plot(fr1,median(data1a)-std(data1a),'color',cols(d),'linewidth',1)
    axis tight
    
    subplot(4,4,d+4)
    plot(fr1,median(data1b),'color',cols(d),'linewidth',2)
    patch([fr1 fliplr(fr1)],[median(data1b)-std(data1b) fliplr(median(data1b)+std(data1b))],cols(d),'facealpha',0.4)
    %     plot(fr1,median(data1b)+std(data1b),'color',cols(d),'linewidth',1)
    %     plot(fr1,median(data1b)-std(data1b),'color',cols(d),'linewidth',1)
    title('non-smoothed norm(FFTresp)/norm(FFTstim)')
    axis tight
    subplot(4,4,8)
    plot(fr1,median(data1b),'color',cols(d),'linewidth',2)
    hold on
    patch([fr1 fliplr(fr1)],[median(data1b)-std(data1b) fliplr(median(data1b)+std(data1b))],cols(d),'facealpha',0.2)
    %     plot(fr1,median(data1b)+std(data1b),'color',cols(d),'linewidth',1)
    %     plot(fr1,median(data1b)-std(data1b),'color',cols(d),'linewidth',1)
    axis tight
    
    subplot(4,4,d+8)
    plot(fr2,median(data2),'color',cols(d),'linewidth',2)
    patch([fr2 fliplr(fr2)],[median(data2)-std(data2) fliplr(median(data2)+std(data2))],cols(d),'facealpha',0.4)
    %     plot(fr2,median(data2)+std(data2),'color',cols(d),'linewidth',1)
    %     plot(fr2,median(data2)-std(data2),'color',cols(d),'linewidth',1)
    title('smoothed norm(FFTresp)/norm(FFTstim)')
    axis tight
    
    subplot(4,4,12)
    plot(fr2,median(data2),'color',cols(d),'linewidth',2)
    hold on
    patch([fr2 fliplr(fr2)],[median(data2)-std(data2) fliplr(median(data2)+std(data2))],cols(d),'facealpha',0.2)
    
    %     plot(fr2,median(data2)+std(data2),'color',cols(d),'linewidth',1)
    %     plot(fr2,median(data2)-std(data2),'color',cols(d),'linewidth',1)
    axis tight
    
    figure(2)
    
    plot(fr2,median(data2),'color',cols(d),'linewidth',2)
    hold on
    %     patch([fr2 fliplr(fr2)],[median(data2)-std(data2) fliplr(median(data2)+std(data2))],cols(d),'facealpha',0.2)
    axis tight
    
    figure(1)
    subplot(4,4,d+12)
    plot(fr2,median(data3),'color',cols(d),'linewidth',2)
    patch([fr2 fliplr(fr2)],[median(data3)-std(data3) fliplr(median(data3)+std(data3))],cols(d),'facealpha',0.4)
    %     plot(fr2,median(data2)+std(data2),'color',cols(d),'linewidth',1)
    %     plot(fr2,median(data2)-std(data2),'color',cols(d),'linewidth',1)
    title('norm(smoothed nnorm(FFTresp)/norm(FFTstim))')
    axis tight
    
    subplot(4,4,16)
    plot(fr2,median(data3),'color',cols(d),'linewidth',2)
    hold on
    patch([fr2 fliplr(fr2)],[median(data3)-std(data3) fliplr(median(data3)+std(data3))],cols(d),'facealpha',0.2)
    
    %     plot(fr2,median(data2)+std(data2),'color',cols(d),'linewidth',1)
    %     plot(fr2,median(data2)-std(data2),'color',cols(d),'linewidth',1)
    axis tight
    
    colldata2{d}=data2; %smooth resp/stim (mouse/pig/human)
    colldata3{d}=data3; %norm-smooth resp/stim (mouse/pig/human)
    colldata1a{d}=data1a; %non-smooth rsp/stim
    colldata1b{d}=data1b; %non-smooth norm(resp)/norm(stim)
    togos{d}=togo;
end
set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(superpath,'newAttempt2',['chirp_comparison',int2str(nlist(nfft)),'.bmp']))
% saveas(gcf,fullfile(superpath,'newAttempt2',['chirp_comparison',int2str(nlist(nfft)),'.fig']))
% save(fullfile(superpath,'newAttempt2','chirp',['chirpRatios_',int2str(nlist(nfft))]),'colldata2','colldata3','colldata1a','colldata1b','frs1','frs2','togos')

figure(2)
title(['w: ',int2str(length(togos{1})),', p: ',int2str(length(togos{2})),' ,h: ',int2str(length(togos{3}))])
set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(superpath,'newAttempt2',['chirp_comparison',int2str(nlist(nfft)),'.eps']))
%% spatial and temporal tuning (normalized F1)
% plot for each temporal/spatial frequency the spatial/temporal tuning
% curve based on normalized F1 peaks (normalized for each cell to its max.
% F1)
% blue=mouse, green=pig, red=human
% including STD (thin line), and 95%CI (thick line)
% version1: take always all cells which have been classified as
% "DG-positive" (at least 1 response to one DG combination)
% version2: take for each single compared value only the cells which
% respond
clear
close all
superpath='f:\all_species_temp';
version=1;

col1={'bbcc','ggkk','rrmm'};
col2={'bbbccc','gggkkk','rrrmmm'};
style1={'-','--','-','--'};
style2={'-','--','-.','-','--','-.'};
cols='bgr';
Dtype='wph';
spatials=[100 200 500 1000 2000 4000];
freqs=[1 2 4 8];
harmlist=[0.5 1 2 3];

collNorm=cell(3,1);
NofCT=zeros(3,6,4);
NofCS=zeros(3,4,6);
goodT=cell(3,6,4);
goodS=cell(3,4,6);
% NofC=zeros(3,24);
meanPeakT=zeros(3,6,4);
meanPeakS=zeros(3,4,6);
stT=zeros(3,6,4);
stS=zeros(3,4,6);
ciT=zeros(3,6,4);
ciS=zeros(3,4,6);
for d=1:3
    fcnt=0;
    load(fullfile(superpath,'newAttempt2','DG',[Dtype(d),'_normPeaks']))
    collNorm{d,1}=reshape(normPeaks(:,:,2),size(normPeaks,1),24);
    scnt=0;
    for s=1:4:24
        scnt=scnt+1;
        
        fcnt=fcnt+1;
        figure(fcnt)
        
        if version==2
            meanPeak=[]; st=[]; ci=[]; gcnt=0;
            for g=s:s+3
                gcnt=gcnt+1;
                good=find(normPeaks(:,g,2)>0);
                goodT{d,scnt,gcnt}=good;
                
                meanPeak=[meanPeak mean(normPeaks(good,g,2))];
                st=[st std(normPeaks(good,g,2))];
                ci=[ci 1.96*st(end)/sqrt(length(good))];
                NofCT(d,scnt,gcnt)=length(good);
            end
        else
            meanPeak=mean(normPeaks(:,s:s+3,2));
            st=std(normPeaks(:,s:s+3,2));
            ci=1.96*st/sqrt(size(normPeaks,1));
            for g=1:4
                NofCT(d,scnt,g)=size(normPeaks,1);
                goodT{d,scnt,g}=1:size(normPeaks,1);
            end
        end
        
        %         st=std(normPeaks(:,s:s+3,2));
        %         ci=1.96*st/sqrt(size(normPeaks,1));
        %             normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
        plot((freqs)+(d-1)*0.1,meanPeak,'color',cols(d),'linewidth',2)
        hold on
        for q=1:4
            plot([freqs(q)+(d-1)*0.1 freqs(q)+(d-1)*0.1],[meanPeak(q)-st(q) meanPeak(q)+st(q)],'color',cols(d))
            plot([freqs(q)+(d-1)*0.1 freqs(q)+(d-1)*0.1],[meanPeak(q)-ci(q) meanPeak(q)+ci(q)],'color',cols(d),'linewidth',5)
        end
        %         set(gca,'xtick',1:4)
        %         set(gca,'xticklabel',freqs)
        xlabel('frequency [Hz]')
        axis([0 9 -inf inf])
        title([int2str(spatials(scnt)),'um / f',num2str(harmlist(2))])
        
        
        meanPeakT(d,scnt,:)=meanPeak;
        stT(d,scnt,:)=st;
        ciT(d,scnt,:)=ci;
        
    end
    scnt=0;
    for s=1:4
        fcnt=fcnt+1;
        figure(fcnt)
        scnt=scnt+1;
        
        if version==2
            meanPeak=[]; st=[]; ci=[]; gcnt=0;
            for g=s:4:24
                gcnt=gcnt+1;
                good=find(normPeaks(:,g,2)>0);
                goodS{d,scnt,gcnt}=good;
                
                meanPeak=[meanPeak mean(normPeaks(good,g,2))];
                st=[st std(normPeaks(good,g,2))];
                ci=[ci 1.96*st(end)/sqrt(length(good))];
                NofCS(d,scnt,gcnt)=length(good);
            end
        else
            meanPeak=mean(normPeaks(:,s:4:24,2));
            st=std(normPeaks(:,s:4:24,2));
            ci=1.96*st/sqrt(size(normPeaks,1));
            for g=1:6
                NofCS(d,scnt,g)=size(normPeaks,1);
                goodS{d,scnt,g}=1:size(normPeaks,1);
            end
        end
        
        %         meanPeak=mean(normPeaks(:,s:4:24,2));%contains peaks normalized to each cell's f1
        %         st=std(normPeaks(:,s:4:24,2));
        %                 ci=1.96*st/sqrt(size(normPeaks,1));
        
        %             normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
        semilogx(([100 200 500 1000 2000 4000])+spatials*0.1,meanPeak,'color',cols(d),'linewidth',2)
        hold on
        %         patch([1:6 6:-1:1],[meanPeak-st fliplr(meanPeak+st)],cols(d),'facealpha',0.4)
        for q=1:6
            semilogx([spatials(q)+spatials(q)*0.1 spatials(q)+spatials(q)*0.1],[meanPeak(q)-st(q) meanPeak(q)+st(q)],'color',cols(d))
            semilogx([spatials(q)+spatials(q)*0.1 spatials(q)+spatials(q)*0.1],[meanPeak(q)-ci(q) meanPeak(q)+ci(q)],'color',cols(d),'linewidth',5)
            
        end
        
        %         set(gca,'xtick',1:6)
        %         set(gca,'xticklabel',spatials)
        xlabel('spatial period [um]')
        axis([50 5000 -inf inf])
        title([int2str(freqs(s)),'Hz / f',num2str(harmlist(2))])
        
        meanPeakS(d,scnt,:)=meanPeak;
        stS(d,scnt,:)=st;
        ciS(d,scnt,:)=ci;
        
    end
    
end



for d=1:3
    figure(100*d)
    for fcnt=1:6
        plot((freqs)+(fcnt-7)*0.1,reshape(meanPeakT(d,fcnt,:),4,1),'color',col2{d}(fcnt),'linewidth',2,'linestyle',style2{fcnt})
        hold on
    end
    legend(num2str(spatials'))
    for fcnt=1:6
        for q=1:4
            plot([freqs(q)+(fcnt-7)*0.1 freqs(q)+(fcnt-7)*0.1],[reshape(meanPeakT(d,fcnt,q),1,1)-reshape(stT(d,fcnt,q),1,1) reshape(meanPeakT(d,fcnt,q),1,1)+reshape(stT(d,fcnt,q),1,1)],'color',col2{d}(fcnt),'linestyle',style2{fcnt})
            plot([freqs(q)+(fcnt-7)*0.1 freqs(q)+(fcnt-7)*0.1],[reshape(meanPeakT(d,fcnt,q),1,1)-reshape(ciT(d,fcnt,q),1,1) reshape(meanPeakT(d,fcnt,q),1,1)+reshape(ciT(d,fcnt,q),1,1)],'color',col2{d}(fcnt),'linewidth',3,'linestyle',style2{fcnt})
            
        end
    end
    
    %     set(gca,'xtick',1:4)
    %     set(gca,'xticklabel',freqs)
    xlabel('frequency [Hz]')
    axis([0 9 -inf inf])
    title(['temporal tuning (',Dtype(d),')'])
    
    
    
    figure(1000*d)
    for fcnt=1:4
        semilogx((spatials)+spatials*0.05*(fcnt-3),reshape(meanPeakS(d,fcnt,:),6,1),'color',col1{d}(fcnt),'linewidth',2,'linestyle',style1{fcnt})
        hold on
    end
    legend(num2str(freqs'))
    for fcnt=1:4
        for q=1:6
            semilogx([spatials(q)+spatials(q)*0.05*(fcnt-3) spatials(q)+spatials(q)*0.05*(fcnt-3)],[reshape(meanPeakS(d,fcnt,q),1,1)-reshape(stS(d,fcnt,q),1,1) reshape(meanPeakS(d,fcnt,q),1,1)+reshape(stS(d,fcnt,q),1,1)],'color',col1{d}(fcnt),'linestyle',style1{fcnt})
            semilogx([spatials(q)+spatials(q)*0.05*(fcnt-3) spatials(q)+spatials(q)*0.05*(fcnt-3)],[reshape(meanPeakS(d,fcnt,q),1,1)-reshape(ciS(d,fcnt,q),1,1) reshape(meanPeakS(d,fcnt,q),1,1)+reshape(ciS(d,fcnt,q),1,1)],'color',col1{d}(fcnt),'linewidth',3,'linestyle',style1{fcnt})
        end
        
    end
    %     set(gca,'xtick',1:6)
    %     set(gca,'xticklabel',spatials)
    xlabel('spatial period [um]')
    axis([50 5000 -inf inf])
    title(['spatial tuning (',Dtype(d),')'])
end


for d=1:3
    tdata=collNorm{d};
    
    if version==2
        goods=reshape(goodT(d,:,:),24,1);
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            curr=curr(goods{m});
            meandata=[meandata mean(curr)];
        end
    else
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            meandata=[meandata mean(curr)];
        end
    end
    normthis=meandata/max(meandata);
    zi=reshape(normthis,length(unique(freqs)),length(unique(spatials)));
    zi=zi';
    zi=flipud(zi);
    figure(10000*d)
    set(gcf,'position',[1 31 1600 794])
    imagesc(1:length(unique(freqs)),1:length(unique(spatials)),zi)
    colormap(gray)
    set(gca,'YTick',1:length(unique(spatials)))
    set(gca,'XTick',1:length(unique(freqs)))
    set(gca,'XTickLabel',unique(freqs))
    xlabel('Hz')
    set(gca,'YTickLabel',fliplr(unique(spatials)))
    ylabel('period (um)')
end


statsS=zeros(3,4,6);%m-p,m-h, p-h
statsT=zeros(3,6,4); %6 spatial tunings, each 4 temporal freqs
cond=[1 2;1 3;2 3];
scnt=0;
for s=1:4:24 %go through temporal tuning
    scnt=scnt+1;
    sscnt=0;
    for ss=s:s+3
        sscnt=sscnt+1;
        for c=1:3
            if version==2
                data1=collNorm{cond(c,1)}(goodT{cond(c,1),scnt,sscnt},ss);
                data2=collNorm{cond(c,2)}(goodT{cond(c,2),scnt,sscnt},ss);
            else
                data1=collNorm{cond(c,1)}(:,ss);
                data2=collNorm{cond(c,2)}(:,ss);
            end
            p=ranksum(data1,data2);
            
            statsT(c,scnt,sscnt)=p; %comparison, figure,data point
        end
    end
end
scnt=0;
for s=1:4 %go through temporal tuning
    scnt=scnt+1;
    sscnt=0;
    for ss=s:4:24
        sscnt=sscnt+1;
        for c=1:3
            if version==2
                data1=collNorm{cond(c,1)}(goodS{cond(c,1),scnt,sscnt},ss);
                data2=collNorm{cond(c,2)}(goodS{cond(c,2),scnt,sscnt},ss);
            else
                data1=collNorm{cond(c,1)}(:,ss);
                data2=collNorm{cond(c,2)}(:,ss);
            end
            p=ranksum(data1,data2);
            
            statsS(c,scnt,sscnt)=p; %comparison, figure,data point
        end
    end
end



condlist={'w-p','w-h','p-h'};

fcnt=0;
for s=1:4:24
    
    fcnt=fcnt+1;
    figure(fcnt)
    yy=get(gca,'ytick');
    y=yy(end);
    icnt=0;
    for i=s:s+3
        icnt=icnt+1;
        for c=1:3
            curr=statsT(c,fcnt,icnt);
            curr=round(curr*1000)/1000;
            text(freqs(icnt),y+(y/8)*c,num2str(curr))
            text(freqs(icnt),y+(y/8)*(c+7),num2str(NofCT(c,fcnt,icnt)))
            if i==s
                text(0.5,y+(y/8)*c,condlist{c})
                text(0.5,y+(y/8)*(c+7),Dtype(c))
            end
        end
    end
    
    axis([0 9 yy(1) y+(y/8)*11])
end


% fcnt=0;
ssspatials=[100 300 600 1000 2000 4000];
for s=1:4
    
    fcnt=fcnt+1;
    figure(fcnt)
    yy=get(gca,'ytick');
    y=yy(end);
    icnt=0;
    for i=s:4:24
        icnt=icnt+1;
        for c=1:3
            curr=statsS(c,fcnt-6,icnt);
            curr=round(curr*1000)/1000;
            text(ssspatials(icnt),y+(y/8)*c,num2str(curr))
            text(ssspatials(icnt),y+(y/8)*(c+7),num2str(NofCS(c,fcnt-6,icnt)))
            if i==s
                text(10,y+(y/8)*c,condlist{c})
                text(10,y+(y/8)*(c+7),Dtype(c))
            end
        end
    end
    axis([50 5000 yy(1) y+(y/8)*11])
end





for f=1:6
    figure(f)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.emf']))
    
end
for f=1:4
    figure(f+6)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.emf']))
    
end

for d=1:3
    figure(100*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.emf']))
end
for d=1:3
    figure(10000*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_heatmap_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_heatmap_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
end
%% speed preference: tuning curves
% peaks are baseline corrected
clear
close all
superpath='f:\all_species_temp';
savepath=fullfile(superpath,'overview_pics');


names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[1; 2; 4; 8; 16];
types='wph';
cols='bgr';
titles={'mouse','pig','human'};

figure
clc
colldata=cell(4,3);
for t=1:length(types)
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    
    
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    load(fullfile(superpath,[Dtype,'_peaks6vel_black']))
    peaksB=peak6vel(:,2:6);
    load(fullfile(superpath,[Dtype,'_peaks6vel_white']))
    peaksW=peak6vel(:,2:6);
    
    peaksB=peaksB(goodB,:);
    su=sum(peaksB,2);
    genGoodB=find(su>0);
    peaksB=peaksB(genGoodB,:);
    goodB=goodB(genGoodB);
    peaksW=peaksW(goodW,:);
    su=sum(peaksW,2);
    genGoodW=find(su>0);
    peaksW=peaksW(genGoodW,:);
    goodW=goodW(genGoodW);
    
    display([types(t),': ',int2str(length(goodW)),' white, ',int2str(length(goodB)),' black, ',int2str(length(unique([goodB;goodW]))),' all'])
    
    %black
    currdata=peaksB./repmat(max(peaksB,[],2),1,5);
    me=mean(currdata);
    st=std(currdata);
    ci=1.96*st/sqrt(size(currdata,1));
    subplot(2,2,1)
    plot(speeds,me,'-o','color',cols(t))
    hold on
    for s=1:length(speeds)
        plot([speeds(s) speeds(s)],[me(s)-st(s) me(s)+st(s)],'-','color',cols(t),'linewidth',2)
        plot([speeds(s) speeds(s)],[me(s)-ci(s) me(s)+ci(s)],'-','color',cols(t),'linewidth',4)
    end
    axis([0 17 0 1.6])
    title('black')
    colldata{1,t}=currdata;
    
    %white
    currdata=peaksW./repmat(max(peaksW,[],2),1,5);
    me=mean(currdata);
    st=std(currdata);
    ci=1.96*st/sqrt(size(currdata,1));
    subplot(2,2,2)
    plot(speeds,me,'-o','color',cols(t))
    hold on
    for s=1:length(speeds)
        plot([speeds(s) speeds(s)],[me(s)-st(s) me(s)+st(s)],'-','color',cols(t),'linewidth',2)
        plot([speeds(s) speeds(s)],[me(s)-ci(s) me(s)+ci(s)],'-','color',cols(t),'linewidth',4)
    end
    axis([0 17 0 1.6])
    title('white')
    colldata{2,t}=currdata;
    
    %all
    allids=unique([goodB;goodW]);
    idsB=setdiff(goodB,goodW);
    idsW=setdiff(goodW,goodB);
    all=intersect(goodB,goodW);
    peaks=zeros(length(allids),5);
    for p=1:length(idsB)
        tmp=find(allids==idsB(p));
        tmp2=find(goodB==idsB(p));
        peaks(tmp,:)=peaksB(tmp2,:);
    end
    for p=1:length(idsW)
        tmp=find(allids==idsW(p));
        tmp2=find(goodW==idsW(p));
        peaks(tmp,:)=peaksW(tmp2,:);
    end
    for p=1:length(all)
        tmp=find(allids==all(p));
        tmp2=find(goodB==all(p));
        tmp3=find(goodW==all(p));
        Bdata=peaksB(tmp2,:);
        Wdata=peaksW(tmp3,:);
        maB=max(Bdata); maW=max(Wdata);
        [ma id]=max([maB maW]);
        if id==1
            peaks(tmp,:)=peaksB(tmp2,:);
        else
            peaks(tmp,:)=peaksW(tmp3,:);
        end
    end
    
    
    
    
    currdata=peaks./repmat(max(peaks,[],2),1,5);
    me=mean(currdata);
    st=std(currdata);
    ci=1.96*st/sqrt(size(currdata,1));
    subplot(2,2,3)
    plot(speeds,me,'-o','color',cols(t))
    hold on
    for s=1:length(speeds)
        plot([speeds(s) speeds(s)],[me(s)-st(s) me(s)+st(s)],'-','color',cols(t),'linewidth',2)
        plot([speeds(s) speeds(s)],[me(s)-ci(s) me(s)+ci(s)],'-','color',cols(t),'linewidth',4)
    end
    axis([0 17 0 1.6])
    title('bigger max')
    colldata{3,t}=currdata;
    
    
    
    %all2
    allval=zeros(size(info,1),2);
    black=cell2mat(info(goodB,24));
    white=cell2mat(info(goodW,25));
    allval(goodB,1)=black;
    allval(goodW,2)=white;
    [allval ids]=max(allval,[],2);
    tmp=find(allval>0);
    ids=ids(tmp);
    
    un=unique([goodB;goodW]);
    peaks=zeros(length(ids),5);
    for ii=1:length(ids)
        if ids(ii)==1
            tmp=find(goodB==un(ii));
            peaks(ii,:)=peaksB(tmp,:);
        else
            tmp=find(goodW==un(ii));
            peaks(ii,:)=peaksW(tmp,:);
        end
    end
    
    
    
    currdata=peaks./repmat(max(peaks,[],2),1,5);
    me=mean(currdata);
    st=std(currdata);
    ci=1.96*st/sqrt(size(currdata,1));
    subplot(2,2,4)
    plot(speeds,me,'-o','color',cols(t))
    hold on
    for s=1:length(speeds)
        plot([speeds(s) speeds(s)],[me(s)-st(s) me(s)+st(s)],'-','color',cols(t),'linewidth',2)
        plot([speeds(s) speeds(s)],[me(s)-ci(s) me(s)+ci(s)],'-','color',cols(t),'linewidth',4)
    end
    axis([0 17 0 1.6])
    title('bigger 50% cum speed (same data as histo)')
    colldata{4,t}=currdata;
    
end

pvals=zeros(4,3,5);

for c=1:4
    cnt=0;
    for cc=1:2
        for ccc=cc+1:3
            cnt=cnt+1;
            for s=1:5
                data1=colldata{c,cc}(:,s);
                data2=colldata{c,ccc}(:,s);
                p=ranksum(data1,data2);
                pvals(c,cnt,s)=p;
            end
        end
    end
end

for c=1:4
    subplot(2,2,c)
    for cc=1:3
        for s=1:5
            text(speeds(s),1.2+0.1*cc,num2str(round(pvals(c,cc,s)*100)/100))
        end
    end
end
%% plot mean heatmaps for cluster result
clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
% path1='e:\all_species_temp\newAttemptOnlyMouse\DG_PCanalysis';
path1='D:\_data\all_species\newAttemptOnlyHum\DG_PCanalysis';
file='clust_Kmeans_PC_F1';
load(fullfile(path1,file))
sel=1; clust=12; %actual clust
gau=1;
st=0;
fact=2;

if sel==1
    data=ID(:,clust-1);
end

% superpath='e:\all_species_temp';
superpath='D:\_data\all_species';
Dtype='wph';
normDG=[];
for ty=1:3
    load(fullfile(superpath,'newAttempt6','DG',[Dtype(ty),'_normPeaks']))
    normDG=[normDG;normPeaks];
end

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0.2 0.2 0.2; 0 1 1; 1 1 0;0.5 0 0];

clu=unique(data);
for clus=1:length(clu)
    currID=find(data==clu(clus));
    
    figure(2)
    subplot(4,5,clus)
    curr=normDG(currID,:,2);
    curr=mean(curr,1);
    curr=fliplr(reshape(curr,4,6));
    [y,x] = meshgrid(1:6,1:4);
    imagesc(x(:)',y(:)',curr');
    colormap('gray')
    title([int2str(subplot(4,5,clus)),' (mean; ',int2str(length(currID)),'cells)'])
    figure(3)
    subplot(4,5,clus)
    curr=normDG(currID,:,2);
    curr=std(curr,0,1);
    curr=fliplr(reshape(curr,4,6));
    [y,x] = meshgrid(1:6,1:4);
    imagesc(curr'*100,[0 100]);
    colormap('gray')
    title([int2str(clus),' (std)'])
    
    
    if gau==1
        
        figure(4)
        subplot(4,5,clus)
        curr=normDG(currID,:,2);
        curr=mean(curr,1);
        curr=fliplr(reshape(curr,4,6));
        [y,x] = meshgrid(1:6,1:4);
        imagesc(x(:)',y(:)',curr');
        hold on
        colormap('gray')
        title([int2str(clus),' (mean; ',int2str(length(currID)),'cells)'])
        
        
        figure(1)
        
        curr=normDG(currID,:,2);
        curr=mean(curr,1);
        dat=fliplr(reshape(curr,4,6));
        
        minData = mean(mean(dat))+st*std(reshape(dat,24,1));
        p = (dat - minData);
        tmp=find(p<0);
        p(tmp)=0;
        sumData = sum(sum(p));
        p=p/sumData;
        
        
        [y,x] = meshgrid(6:-1:1,1:4);
        xx=reshape(x,24,1);
        yy=reshape(y,24,1);
        pp=reshape(p,24,1);
        mx=dot(xx,pp); %centerpoint
        my=dot(yy,pp); %centerpoint
        
        a=dot((xx-mx).^2,pp);
        b=dot((xx-mx).*(yy-my),pp);
        c=dot((yy-my).^2,pp);
        co=[a,b;b,c];
        
        if fact>=1
            D=fact*sqrt(eigs(co));
        else
            D=sqrt(eigs(co));%Eigenvalues = two main axes
        end
        [V,dd]=eigs(co);
        ma=max(max(dd));
        [ma1 ma2]=find(dd==ma);
        %         ma1
        V=V(ma1,:); v=V(1)/V(2);
        
        Ang=atan2(V(2),V(1)); %used
        %                Ang=atan2(V(1),V(2));
        AngD=rad2deg(Ang); %angle in degrees
        
        ymin=get(gca,'ylim');
        ymin=ymin(1);
        
        centerX=mx;%x-axis center;1=left
        centerY=my;%y-axis center; 1=bottom
        
        
        
        centerX2=mx;%x-axis center;1=left
        centerY2=7-my;%y-axis center; 1=bottom
        Ang2=-Ang;
        [X,Y]=meshgrid(0.5:0.01:4.5,0.5:0.01:6.5);
        
        theta = pi/2-Ang;
        %
        aa = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
        bb = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
        cc = sin(theta)^2/2/(D(1)/1.5)^2 + cos(theta)^2/2/(D(2)/1.5)^2;
        allClust = exp( - (aa*(Y-centerY).^2 + 2*bb*(Y-centerY).*(X-centerX2) + cc*(X-centerX2).^2)) ;
        subplot(4,5,clus)
        contour(X,Y,allClust(:,:),1,'color','r')
        
        hold on
        plot(centerX,centerY,'or','markerfacecolor','r','markersize',2)
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',[1 2 4 8])
        set(gca,'ytick',1:6)
        set(gca,'yticklabel',[100 200 500 1000 2000 4000])
        axis([0.5 4.5 0.5 6.5])
        
        
        
        
        
        figure(4)
        curr=normDG(currID,:,2);
        curr=mean(curr,1);
        dat=fliplr(reshape(curr,4,6));
        dat=fliplr(dat);
        
        minData = mean(mean(dat))+st*std(reshape(dat,24,1));
        p = (dat - minData);
        tmp=find(p<0);
        p(tmp)=0;
        sumData = sum(sum(p));
        p=p/sumData;
        
        
        [y,x] = meshgrid(6:-1:1,1:4);
        xx=reshape(x,24,1);
        yy=reshape(y,24,1);
        pp=reshape(p,24,1);
        mx=dot(xx,pp); %centerpoint
        my=dot(yy,pp); %centerpoint
        
        a=dot((xx-mx).^2,pp);
        b=dot((xx-mx).*(yy-my),pp);
        c=dot((yy-my).^2,pp);
        co=[a,b;b,c];
        
        if fact>=1
            D=fact*sqrt(eigs(co));
        else
            D=sqrt(eigs(co));%Eigenvalues = two main axes
        end
        [V,dd]=eigs(co);
        ma=max(max(dd));
        [ma1 ma2]=find(dd==ma);
        %         ma1
        V=V(ma1,:); v=V(1)/V(2);
        
        Ang=atan2(V(2),V(1)); %used
        %                Ang=atan2(V(1),V(2));
        AngD=rad2deg(Ang); %angle in degrees
        
        ymin=get(gca,'ylim');
        ymin=ymin(1);
        
        centerX=mx;%x-axis center;1=left
        centerY=my;%y-axis center; 1=bottom
        
        
        
        centerX2=mx;%x-axis center;1=left
        centerY2=7-my;%y-axis center; 1=bottom
        Ang2=-Ang;
        [X,Y]=meshgrid(0.5:0.01:4.5,0.5:0.01:6.5);
        
        theta = pi/2-Ang;
        %
        aa = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
        bb = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
        cc = sin(theta)^2/2/(D(1)/1.5)^2 + cos(theta)^2/2/(D(2)/1.5)^2;
        allClust = exp( - (aa*(Y-centerY).^2 + 2*bb*(Y-centerY).*(X-centerX2) + cc*(X-centerX2).^2)) ;
        contour(X,Y,allClust(:,:),1,'color','r')
        hold on
        plot(centerX,centerY,'or','markerfacecolor','r','markersize',2)
        
    end
    
    
    
end

%% plot mean heatmaps for cluster result: separate for each species
clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
path1='e:\all_species_temp\newAttemptOnlyMouse\DG_PCanalysis';

file='clust_Kmeans_PC_F1';
load(fullfile(path1,file))
sel=1; clust=12; %actual clust
gau=1;
st=0;
fact=2;

if sel==1
    data=ID(:,clust-1);
end

superpath='e:\all_species_temp';
Dtype='wph';
cols='bgry';
normDG=[];
for ty=1:3
    load(fullfile(superpath,'newAttempt3','DG',[Dtype(ty),'_normPeaks']))
    normDG=[normDG;normPeaks];
end

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0.2 0.2 0.2; 0 1 1; 1 1 0;0.5 0 0];

clu=unique(data);
for clus=1:length(clu)
    currID=find(data==clu(clus));
    
    if gau==1
        
        for t=1:4
            figure(t)
            subplot(4,4,clus)
            if t==4
                tst=1:length(currID);
            else
                tst=find(tID(currID)==t);
            end
            
            if ~isempty(tst)
                curr=normDG(currID(tst),:,2);
                curr=mean(curr,1);
                curr=flipud(reshape(curr,4,6));
                [x,y] = meshgrid(1:6,1:4);
                imagesc(x(:)',y(:)',curr);
                hold on
                colormap('gray')
                title([int2str(clus),' (mean; ',int2str(length(tst)),'cells)'])
                
                
                
                curr=normDG(currID(tst),:,2);
                curr=mean(curr,1);
                dat=reshape(curr,4,6);
                dat=flipud(fliplr(dat));
                
                minData = mean(mean(dat))+st*std(reshape(dat,24,1));
                p = (dat - minData);
                tmp=find(p<0);
                p(tmp)=0;
                sumData = sum(sum(p));
                p=p/sumData;
                
                
                [x,y] = meshgrid(6:-1:1,1:4);
                xx=reshape(x,24,1);
                yy=reshape(y,24,1);
                pp=reshape(p,24,1);
                mx=dot(xx,pp); %centerpoint
                my=dot(yy,pp); %centerpoint
                
                a=dot((xx-mx).^2,pp);
                b=dot((xx-mx).*(yy-my),pp);
                c=dot((yy-my).^2,pp);
                co=[a,b;b,c];
                
                if fact>=1
                    D=fact*sqrt(eigs(co));
                else
                    D=sqrt(eigs(co));%Eigenvalues = two main axes
                end
                [V,dd]=eigs(co);
                ma=max(max(dd));
                [ma1 ma2]=find(dd==ma);
                %         ma1
                V=V(ma1,:); v=V(1)/V(2);
                
                Ang=atan2(V(2),V(1)); %used
                %                Ang=atan2(V(1),V(2));
                AngD=rad2deg(Ang); %angle in degrees
                
                ymin=get(gca,'ylim');
                ymin=ymin(1);
                
                centerX=mx;%x-axis center;1=left
                centerY=my;%y-axis center; 1=bottom
                
                
                
                centerX2=mx;%x-axis center;1=left
                centerY2=7-my;%y-axis center; 1=bottom
                Ang2=-Ang;
                [Y,X]=meshgrid(0.5:0.01:4.5,0.5:0.01:6.5);
                
                theta = pi/2-Ang;
                %
                if ~isempty(find(D==0))
                    tmp=find(D==0)
                    D(tmp)=0.1;
                end
                
                aa = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
                bb = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
                cc = sin(theta)^2/2/(D(1)/1.5)^2 + cos(theta)^2/2/(D(2)/1.5)^2;
                allClust = exp( - (aa*(Y-centerY).^2 + 2*bb*(Y-centerY).*(X-centerX2) + cc*(X-centerX2).^2)) ;
                
                contour(X,Y,allClust(:,:),1,'color','r')
                hold on
                plot(centerX,centerY,'or','markerfacecolor','r','markersize',2)
                
                
            end
            
            if t==1
                figure(100)
                subplot(4,4,clus)
                [x,y] = meshgrid(1:6,1:4);
                imagesc(x(:)',y(:)',zeros(6,4));
                colormap('gray')
                hold on
            end
            
            if ~isempty(tst)
                figure(100)
                subplot(4,4,clus)
                contour(X,Y,allClust(:,:),1,'color',cols(t))
                hold on
                plot(centerX,centerY,'o','color',cols(t),'markerfacecolor',cols(t),'markersize',2)
            end
        end
    end
    
    
end
%% pie plots for clusters
clear
close all

path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
load(fullfile(path1,'clust_Kmeans_PC_F1new'))
clust=14;
data=ID(:,clust-1);



figure
for t=1:3
    tot=find(tID==t);
    nums=[];
    for cc=1:clust
        ids=find(data==cc);
        tmp=find(tID(ids)==t);
        nums=[nums;length(tmp)];
    end
    subplot(2,2,t)
    labels={'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
    pie(nums,labels)
    colormap('gray')
    
end

%% re-number clusters and re-plot pies


clear
close all

path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
load(fullfile(path1,'clust_Kmeans_PC_F1new'))
clust=14;
data=ID(:,clust-1);

newID=[3 11 10 12 2 1 5 8 6 14 4 9 13 7];

data2=zeros(size(ID,1),1);
for c=1:clust
    tmp=find(data==c);
    data2(tmp)=newID(c);
end


figure
for t=1:3
    tot=find(tID==t);
    nums=[];
    for cc=1:clust
        ids=find(data2==cc);
        tmp=find(tID(ids)==t);
        nums=[nums;length(tmp)];
    end
    subplot(2,2,t)
    labels={'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
    pie(nums,labels)
    %     colormap('gray')
    
end

data=data2;
save(fullfile(path1,'clust_Kmeans_PC_F1new_corr'),'data','keepOrig','orig2clust','tID')
%% LF polarity
% calculated by finding first "real" peak and then 75% of its hight
% manual correction for random burst induced or very unstable peaks
% based on sigma=60 plots

close all
clear

badIDs=202:250;
path2='F:\all_species_temp\final_pics';

path1='F:\all_species_temp';
list='wph';
cols='bgr';

figure

for l=1:length(list);
    load(fullfile(path1,[list(l),'_info']))
    
    good=find(goodones==1);
    if l==2
        good=setdiff(good,badIDs);
    end
    LF=cell2mat(info(good,9));
    keep2=LF; id2=find(LF~=0);
    LF=length(find(LF~=0));
    lfOFF=find(keep2==-1);
    lfON=find(keep2==1);
    
    
    
    
    subplot(2,2,l)
    h=pie([length(lfOFF) length(lfON) ])
    hText = findobj(h,'Type','text'); % text handles
    percentValues = get(hText,'String'); % percent values
    str = {'off ';'on ';}; % text
    combinedstrings = strcat(str,percentValues);
    set(hText,{'String'},combinedstrings);
    title(['polartiy based on LF responses ',list(l),' (total: ',int2str(LF),' cells)'])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'polarity_LF.emf'))
close

%% 6vel per exp
clear
close all
superpath='f:\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

allmeans=[2.0015 3.1624 3.7842];

names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[0.2; 1; 2; 4; 8; 16];
types='wph';
cols='bgr';
col=colormap(hsv);
close all
col=col(1:4:end,:);
titles={'mouse','pig','human'};

collMean=zeros(20,3); NofE=zeros(3,1);
alldata=cell(3,20);
for t=1:length(types)
    
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    
    
    
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    allval=zeros(size(info,1),2);
    black=cell2mat(info(goodB,24));
    white=cell2mat(info(goodW,25));
    allval(goodB,1)=black;
    allval(goodW,2)=white;
    allvalorig=max(allval,[],2);
    
    dates=info(:,1);
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    start=[start;length(allvalorig)+1];
    udates=udates(ia);
    ends=start(2:end)-1;
    [udates ai]=sort(udates);
    start=start(ai);
    ends=ends(ai);
    NofE(t)=length(udates);
    
    su=length(udates)/4;
    if mod(su,1)==0
        su=su+1;
    else
        su=ceil(su);
    end
    collmeds=[];
    for u=1:length(udates)
        allval=allvalorig(start(u):ends(u));
        tmp=find(allval>0);
        allval=allval(tmp);
        alldata{t,u}=allval;
        allval=[allval;-1;17];
        
        
        
        
        
        if length(allval)>2
            figure(t)
            subplot(4,su,u)
            hist(allval,38)
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
            med=median(allval(1:end-2));
            collMean(u,t)=med;
            st=std(allval(1:end-2));
            conf = 1.96*st/sqrt(length(allval)-2);
            title([udates{u}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(allval)-2),')'])
            hold on
            plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
            plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
            yy=get(gca,'ylim');
            plot([med-conf med-conf],[0 yy(2)+1],'--k')
            plot([med+conf med+conf],[0 yy(2)+1],'--k')
            %     set(gca,'xtick',speeds)
            % set(gca,'xticklabel',speeds)
            xlabel('speed (mm/s)')
            axis([0 16 0 yy(2)+1])
            
            figure(t)
            subplot(4,su,su*4)
            c=cdfplot(allval(1:end-2))
            set(c,'color',cols(t),'linestyle','-','linewidth',2)
            hold on
            xlabel('speed (mm/s)')
            ylabel('probability')
            title('better bar')
            figure(4)
            subplot(1,2,1)
            c=cdfplot(allval(1:end-2))
            set(c,'color',cols(t),'linestyle','-','linewidth',2)
            hold on
            xlabel('speed (mm/s)')
            ylabel('probability')
            title('better bar')
            figure(4)
            subplot(1,2,2)
            plot(t,med,'o','markerfacecolor',cols(t),'markeredgecolor',cols(t))
            hold on
            %             plot([t t],[med-st med+st],'-','color',cols(t))
            %             plot([t t],[med-conf med+conf],'-','color',cols(t),'linewidth',2)
            axis([0 4 0 inf])
            collmeds=[collmeds;med];
        end
    end
    plot([t-0.2 t+0.2],[mean(collmeds) mean(collmeds)],'-','color',cols(t),'linewidth',2) %mean of means
    plot([t t],[mean(collmeds)-std(collmeds) mean(collmeds)+std(collmeds)],'-','color',cols(t)) %mean of means
    plot([t-0.2 t+0.2],[allmeans(t) allmeans(t)],'-','color',cols(t)) %previously calculated population mean
    subplot(1,2,1)
    plot(mean(collmeds), 0.5,'o','color',cols(t),'markerfacecolor',cols(t),'markersize',10)
    plot(allmeans(t),0.5,'o','color',cols(t),'markersize',10)
end
for t=1:length(types)
    uucnt=0;
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    allval=zeros(size(info,1),2);
    black=cell2mat(info(goodB,24));
    white=cell2mat(info(goodW,25));
    allval(goodB,1)=black;
    allval(goodW,2)=white;
    allvalorig=max(allval,[],2);
    
    dates=info(:,1);
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    start=[start;length(allvalorig)+1];
    udates=udates(ia);
    ends=start(2:end)-1;
    [udates ai]=sort(udates);
    start=start(ai);
    ends=ends(ai);
    NofE(t)=length(udates);
    
    leg=cell(NofE(t),1);
    
    su=length(udates)/4;
    if mod(su,1)==0
        su=su+1;
    else
        su=ceil(su);
    end
    collmeds=[];
    for u=1:length(udates)
        allval=allvalorig(start(u):ends(u));
        tmp=find(allval>0);
        allval=allval(tmp);
        allval=[allval;-1;17];
        
        
        
        
        if length(allval)>2
            uucnt=uucnt+1;
            leg{uucnt}=[udates{u},' / ',num2str(median(allval(1:end-2)))];
            figure(100+t)
            c=cdfplot(allval(1:end-2))
            set(c,'color',col(u,:),'linestyle','-','linewidth',2)
            hold on
            xlabel('speed (mm/s)')
            ylabel('probability')
            title('better bar')
            
            axis([0 16 0 1])
            collmeds=[collmeds;med];
        end
    end
    leg=leg(1:uucnt);
    legend(leg)
end

p1=ranksum(collMean(1:NofE(1),1),collMean(1:NofE(2),2)); %means for M vs P
p2=ranksum(collMean(1:NofE(1),1),collMean(1:NofE(3),3));  %means for M vs H
p3=ranksum(collMean(1:NofE(2),2),collMean(1:NofE(3),3));  %means for P vs H
% set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(savepath,'speed_all_noSlow.bmp'))
% saveas(gcf,fullfile(savepath,'speed_all_noSlow.fig'))
% close

% within species
goodExps=cell(3,1); pall=cell(3,1);

for t=1:3
    load(fullfile(superpath,[types(t),'_info']))
    dates=info(:,1);
    udates=unique(dates);
    
    
    if t==2
        gooddates=[];
        for dd=1:length(udates)
            if isempty(regexp(udates{dd},'20140723d'))
                gooddates=[gooddates;dd];
            end
        end
        udates=udates(gooddates);
    end
    alldates=udates;
    
    
    goodExp=[];
    
    for e=1:20
        
        if ~isempty(alldata{t,e})
            goodExp=[goodExp;e];
        end
    end
    udates=udates(goodExp);
    goodExps{t}=goodExp;
    pvals=[];
    for g=1:length(alldates)-1
        for g2=g+1:length(alldates)
            if sum(ismember(udates,alldates{g})==1) && sum(ismember(udates,alldates{g2})==1)
                p=ranksum(alldata{t,g},alldata{t,g2});
            else
                p=-1 ;
            end
            pvals=[pvals;p];
        end
    end
    pall{t}=pvals;
    matP=zeros(length(alldates),length(alldates));
    numOfComp=[0 length(alldates)-1:-1:1];
    for g=1:length(alldates)-1
        if sum(ismember(udates,alldates{g})==1)
            matP(g,g)=1;
            currID1=sum(numOfComp(1:g))+1;
            currID2=sum(numOfComp(1:g+1));
            matP(g,g+1:end)=pvals(currID1:currID2);
            matP(g+1:end,g)=pvals(currID1:currID2);
        else
            matP(g,g)=-1;
            matP(g,g+1:end)=repmat(-1,length(alldates)-g,1);
            matP(g+1:end,g)=repmat(-1,length(alldates)-g,1);
        end
    end
    if sum(ismember(udates,alldates{g+1})==1)
        matP(g+1,g+1)=1;
    else
        matP(g+1,g+1)=-1;
    end
    matP=1-matP;
    tmp0=find(matP==2);
    tmp1=find(matP>0.95);
    tmp2=find(matP>0.99);
    tmp3=find(matP<0.95); matP(tmp3)=0;
    matP(tmp1)=0.8; matP(tmp2)=1; matP(tmp0)=0.4;
    figure
    imagesc(matP,[0 1])
    set(gca,'xticklabel',udates)
    set(gca,'yticklabel',udates)
    colormap('gray')
    title(['p-values across ',types(t),' experiments'])
end
figure
subplot(1,2,1)
for t=1:3
    plot(t,pall{t},'o','color',cols(t))
    hold on
    plot([t-0.2 t+0.2],[median(pall{t}) median(pall{t})],'-','color',cols(t))
end

plot([0 4],[0.05 0.05],'-k')
plot([0 4],[0.01 0.01],'-k')
axis([0 4 0 1])
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'m','p','h'})
title('ranksum p-values between all experiments within each species')



% across species
pcross=cell(3,1); cnt=0;
for t=1:2
    for tt=t+1:3
        pvals=[]; cnt=cnt+1;
        for g=1:length(goodExps{t})
            for g2=1:length(goodExps{tt})
                p=ranksum(alldata{t,goodExps{t}(g)},alldata{tt,goodExps{tt}(g2)});
                pvals=[pvals;p];
            end
        end
        pcross{cnt}=pvals;
        malength=max(length(goodExps{t}),length(goodExps{tt}));
        matP=zeros(length(goodExps{t}),length(goodExps{tt}));
        numOfComp=[0 length(goodExp)-1:-1:1];
        
    end
end


subplot(1,2,2)
for t=1:3
    plot(t,pcross{t},'o','color',cols(t))
    hold on
    plot([t-0.2 t+0.2],[median(pcross{t}) median(pcross{t})],'-','color',cols(t))
end

plot([0 4],[0.05 0.05],'-k')
plot([0 4],[0.01 0.01],'-k')
axis([0 4 0 1])
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'m-p','m-h','p-h'})
title('ranksum p-values between all experiments across species')
%% chirp per exp

% for pig, exclude 20140723d

clear
close all
nfft=4;
nlist=[2000 4000 6000 8000 10000 15000];
superpath='f:\all_species_temp';
load(fullfile(superpath,'newAttempt2','chirp',['chirpRatios_',int2str(nlist(nfft))]))
Dtype='wph';

step=5;
data=colldata2;
fr=frs2{1}(1,:);
dataset='smooth'; %data2
% dataset='noSmooth_noNorm'; %data1a
% dataset='noSmooth_Norm'; %data1b
step=5;
cols=colormap(hsv);
close all
cols=cols(1:4:end,:);
col='bgr';
alldata=cell(3,15);
expNum=zeros(3,1);

for t=1:3
    figure(t)
    load(fullfile(superpath,[Dtype(t),'_info']))
    datanow=data{t};
    tog=togos{t};
    dates=info(tog,1);
    udates=unique(dates);
    if t==2
        for u=1:length(udates)
            if ~isempty(regexp(udates{u},'20140723d'))
                excl=u;
            end
        end
        ulist=1:length(udates);
        ulist=setdiff(ulist,excl);
        udates=udates(ulist);
    end
    expNum(t)=length(udates);
    
    means=[]; leg=cell(length(udates),1);
    for d=1:length(udates)
        start=0;m=0;
        while start==0
            m=m+1;
            if ~isempty(regexp(dates{m},udates{d}))
                start=m;
            end
        end
        stop=0;
        while stop==0 && m<length(dates)
            m=m+1;
            if isempty(regexp(dates{m},udates{d}))
                stop=m-1;
            end
        end
        if stop==0
            stop=length(dates);
        end
        
        currdata=datanow(start:stop,:);
        plot(fr,mean(currdata),'-','color',cols(d,:),'linewidth',2)
        hold on
        leg{d}=[udates{d},' / ',int2str(size(currdata,1)),' cells'];
        means=[means;mean(currdata)];
        alldata{t,d}=currdata;
    end
    plot(fr,mean(means),'--k','linewidth',3)
    plot(fr,mean(datanow),'-k','linewidth',4)
    leg{d+1}='mean of means';
    leg{d+2}='mean of all cells';
    legend(leg)
    title(Dtype(t))
    axis([1 8 0 1])
end

figure
for t=1:3
    load(fullfile(superpath,[Dtype(t),'_info']))
    datanow=data{t};
    tog=togos{t};
    dates=info(tog,1);
    udates=unique(dates);
    if t==2
        for u=1:length(udates)
            if ~isempty(regexp(udates{u},'20140723d'))
                excl=u;
            end
        end
        ulist=1:length(udates);
        ulist=setdiff(ulist,excl);
        udates=udates(ulist);
    end
    expNum(t)=length(udates);
    
    means=[]; leg=cell(length(udates),1);
    for d=1:length(udates)
        start=0;m=0;
        while start==0
            m=m+1;
            if ~isempty(regexp(dates{m},udates{d}))
                start=m;
            end
        end
        stop=0;
        while stop==0 && m<length(dates)
            m=m+1;
            if isempty(regexp(dates{m},udates{d}))
                stop=m-1;
            end
        end
        if stop==0
            stop=length(dates);
        end
        
        currdata=datanow(start:stop,:);
        plot(fr,mean(currdata),'-','color',col(t),'linewidth',1)
        hold on
        means=[means;mean(currdata)];
        alldata{t,d}=currdata;
    end
    plot(fr,mean(means),'--','color',col(t),'linewidth',3)
    plot(fr,mean(datanow),'-','color',col(t),'linewidth',4)
    axis([1 8 0 1])
end

%pvals within
within=zeros(3,53);
for t=1:3
    cnt=0;
    pvals=zeros(1,53);
    for e=1:expNum(t)-1
        for ee=e:expNum(t)
            cnt=cnt+1;
            dat1=alldata{t,e};
            dat2=alldata{t,ee};
            for i=1:53
                now=i:i+step-1;
                if i+step-1 > size(dat1,2)
                    now=i:size(dat1,2);
                end
                p=ranksum(mean(dat1(:,now),2),mean(dat2(:,now),2));
                pvals(cnt,i)=p;
            end
        end
    end
    within(t,:)=median(pvals);
end
figure
plot(fr,within')
hold on
plot([fr(1) fr(end)],[0.05 0.05])
legend('w','p','h')
title('median p-values after comparison of all exp. within a species')
axis([1 8 0 1])

%pvals across
across=zeros(3,53);
tcnt=0;
for t=1:2
    for tt=t+1:3
        tcnt=tcnt+1;
        cnt=0;
        pvals=zeros(1,53);
        for e=1:expNum(t)-1
            for ee=1:expNum(tt)
                cnt=cnt+1;
                dat1=alldata{t,e};
                dat2=alldata{tt,ee};
                for i=1:53
                    p=ranksum(dat1(:,i),dat2(:,i));
                    pvals(cnt,i)=p;
                end
            end
        end
        across(tcnt,:)=median(pvals);
    end
end
figure
plot(fr,across')
hold on
plot([fr(1) fr(end)],[0.05 0.05])
legend('w-p','w-h','p-h')
title('median p-values after comparison of all exp. across species')
axis([1 8 0 1])
%% missing number for overview table
clear
close all

path1='F:\all_species_temp';
types='wph';

FLASH=zeros(3,1);
LF=zeros(3,1);
DS=zeros(3,1);
DG=zeros(3,1);
totCells=zeros(3,1);


for t=1:length(types)
    load(fullfile(path1,[types(t),'_info']))
    goodones=cell2mat(info(:,3));
    goodones=find(goodones>0);
    new=info(goodones,:);
    dates=new(:,1);
    udates=unique(dates);
    
    starts=[];  stops=[];
    for u=1:length(udates)
        m=0; found=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(dates{m},udates{u}))
                found=1;
            end
        end
        starts=[starts;m];
        
        found=0;
        while found==0 && m<length(dates)
            m=m+1;
            if isempty(regexp(dates{m},udates{u}))
                found=1;
            end
        end
        if found==1
            stops=[stops;m-1];
        else
            stops=[stops;m];
        end
    end
    
    
    %flash
    flash=cell2mat(new(:,26));
    ID=find(flash==1);
    for s=1:length(starts)
        tmp=find(ID>=starts(s));
        tmp2=find(ID(tmp)<=stops(s));
        FLASH(t,s)=length(tmp2);
    end
    
    %LF
    ID=cell2mat(new(:,9));
    ID=find(ID~=0);
    for s=1:length(starts)
        tmp=find(ID>=starts(s));
        tmp2=find(ID(tmp)<=stops(s));
        LF(t,s)=length(tmp2);
    end
    
    
    %DS
    ID=cell2mat(new(:,12));
    ID=find(ID==1);
    for s=1:length(starts)
        tmp=find(ID>=starts(s));
        tmp2=find(ID(tmp)<=stops(s));
        DS(t,s)=length(tmp2);
    end
    
    %DG
    load(fullfile(path1,'newAttempt3','DG',[types(t),'_normPeaks']))
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    for s=1:length(starts)
        tmp=find(ID>=starts(s));
        tmp2=find(ID(tmp)<=stops(s));
        DG(t,s)=length(tmp2);
    end
    
    %totcells
    for s=1:length(starts)
        tmp=starts(s):stops(s);
        totCells(t,s)=length(tmp);
    end
end

FLASH=FLASH';
LF=LF';
DG=DG';
DS=DS';
totCells=totCells';

flash=sum(FLASH);
lf=sum(LF);
dg=sum(DG);
ds=sum(DS);
tot=sum(totCells);
%% total responding cells

clear
close all

path1='F:\all_species_temp';
types='wph';

FLASH=zeros(3,1);
LF=zeros(3,1);
DS=zeros(3,1);
DG=zeros(3,1);
totCells=zeros(3,1);

clc
for t=1:length(types)
    load(fullfile(path1,[types(t),'_info']))
    goodones=cell2mat(info(:,3));
    goodones=find(goodones>0);
    
    if t==2
        badIDs=202:250;
        bads=[];
        for b=1:length(badIDs)
            tmp=find(goodones==badIDs(b));
            bads=[bads;tmp];
        end
        goodones=setdiff(goodones,bads);
    end
    
    new=info(goodones,:);
    dates=new(:,1);
    udates=unique(dates);
    
    
    
    %flash
    flash=cell2mat(new(:,26));
    ID=find(flash==1);
    goodIDs=ID;
    
    
    %LF
    ID=cell2mat(new(:,9));
    ID=find(ID~=0);
    goodIDs=[goodIDs;ID];
    
    
    
    %DG
    load(fullfile(path1,'newAttempt3','DG',[types(t),'_normPeaks']))
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    % chirp
    load(fullfile(path1,'newAttempt2','chirp',['chirpRatios_',int2str(8000)]))
    togo=togos{t};
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    
    %6vel
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    ID=[];
    for to=1:length(goodB)
        tmp=find(goodones==goodB(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    ID=[];
    for to=1:length(goodW)
        tmp=find(goodones==goodW(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    
    goodIDs=unique(goodIDs);
    
    
    
    display([types(t),': ',int2str(length(goodIDs)),' out of ',int2str(length(goodones))])
end
%% total responding cells per exp

clear
close all

path1='F:\all_species_temp';
types='wph';

FLASH=zeros(3,1);
LF=zeros(3,1);
DS=zeros(3,1);
DG=zeros(3,1);
totCells=zeros(3,1);
numPerExp=zeros(1,3);


clc
for t=1:length(types)
    load(fullfile(path1,[types(t),'_info']))
    goodones=cell2mat(info(:,3));
    goodones=find(goodones>0);
    
    if t==2
        badIDs=202:250;
        bads=[];
        for b=1:length(badIDs)
            tmp=find(goodones==badIDs(b));
            bads=[bads;tmp];
        end
        goodones=setdiff(goodones,bads);
    end
    
    new=info(goodones,:);
    dates=new(:,1);
    
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    start=[start;length(dates)+1];
    udates=udates(ia);
    ends=start(2:end)-1;
    
    
    
    %flash
    flash=cell2mat(new(:,26));
    ID=find(flash==1);
    goodIDs=ID;
    
    
    %LF
    ID=cell2mat(new(:,9));
    ID=find(ID~=0);
    goodIDs=[goodIDs;ID];
    
    
    
    %DG
    load(fullfile(path1,'newAttempt3','DG',[types(t),'_normPeaks']))
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    % chirp
    load(fullfile(path1,'newAttempt2','chirp',['chirpRatios_',int2str(8000)]))
    togo=togos{t};
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    
    %6vel
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    ID=[];
    for to=1:length(goodB)
        tmp=find(goodones==goodB(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    ID=[];
    for to=1:length(goodW)
        tmp=find(goodones==goodW(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    
    goodIDs=unique(goodIDs);
    
    for u=1:length(udates)
        tmp=find(goodIDs>=start(u));
        tmp2=find(goodIDs(tmp)<=ends(u));
        numPerExp(u,t)=length(tmp2);
    end
    
    display([types(t),': ',int2str(length(goodIDs)),' out of ',int2str(length(goodones))])
end

%% Spatio-temporal per exp
clear
close all

path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
load(fullfile(path1,'clust_Kmeans_PC_F1new_corr'))
types='wph';



for t=1:3
    figure
    load(fullfile('F:\all_species_temp',[types(t),'_info']))
    load(fullfile('F:\all_species_temp\newAttempt3\DG',[types(t),'_normPeaks']))
    dates=info(togo,1);
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    udates=udates(ia);
    start=[start;length(dates)+1];
    
    tot=find(tID==t);
    dataT=data(tot);
    
    for d=1:length(udates)
        nowData=dataT(start(d):start(d+1)-1);
        
        nums=[];
        for cc=1:14
            ids=find(nowData==cc);
            nums=[nums;length(ids)];
        end
        figure(t)
        subplot(4,4,d)
        labels={'1','2','3','4','5','6','7','8','9','10','11','12','13','14'};
        pie(nums,labels)
        
        figure(10+t)
        subplot(4,4,d)
        numlab=cell(14,1);
        for n=1:length(nums)
            numlab{n}=int2str(nums(n));
        end
        pie(nums,numlab)
        title(udates(d))
    end
    figure(t)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('F:\all_species_temp\final_pics',['_clustPercentage_',types(t),'.emf']))
    figure(t+10)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('F:\all_species_temp\final_pics',['_clustNumbers_',types(t),'.emf']))
end

%% find upsilon cells (Petrusca 2007)

close all
clear

path1='F:\all_species_temp\newAttempt3';

load(fullfile(path1,'DG_PCanalysis','clust_Kmeans_PC_F1new_corr'))
load(fullfile(path1,'DG','h_normPeaks'))
load(fullfile('F:\all_species_temp','polarity'))
load(fullfile('F:\all_species_temp','h_vel6_new_Black'))
load(fullfile('F:\all_species_temp','h_chirp'))
load(fullfile('F:\all_species_temp','h_flash60sigma'))
load(fullfile('F:\all_species_temp\newAttempt\DG','h_DGrates'))

dataH=find(tID==3);
dataH=data(dataH);

tmp=find(type==3);
IDh=flashID(tmp);
pol=onoff(tmp);


% all ratios

rat=normPeaks(:,[10 11 12 14 15 16 18 19 20],3)./normPeaks(:,[10 11 12 14 15 16 18 19 20],2);

int=[];
for r=1:9
    tmp=find(rat(:,r)>=1);
    int=[int;tmp];
end
int=unique(int);
clust=dataH(int);
[clust id]=sort(clust);
int=int(id);

%flash
Upol=[];
for c=1:length(clust)
    realID=togo(int(c));
    flash=find(IDh==realID);
    if ~isempty(flash)
        Upol=[Upol;pol(flash)];
        figure
        plot(flash_ratesSmoo2{realID})
        title(int2str(clust(c)))
    else
        Upol=[Upol;0];
    end
end

Upol=Upol(id);


%bar
for c=1:length(clust)
    realID=togo(int(c));
    if ~isempty(vel6_newB{realID})
        figure
        plot(vel6_newB{realID})
        title(int2str(clust(c)))
    end
end



%chirp
for c=1:length(clust)
    realID=togo(int(c))
    if ~isempty(chirpContr{realID})
        figure
        plot(chirpContr{realID})
        title(int2str(clust(c)))
    end
end


%DG
ii=17
tmp=find(clust==1);
keep=int(tmp);
clu=clust(tmp);
tmp=find(clust==6);
keep=[keep;int(tmp)];
clu=[clu;clust(tmp)];
for c=1:length(keep)
    %     realID=togo(int(c));
    if ~isempty(DGRates{keep(c),4,ii})
        figure
        plot(DGRates{keep(c),4,ii})
        title(int2str(clu(c)))
    end
end

ran=randi(293,10,1);
ran=setdiff(ran,clust);

for c=1:length(ran)
    %     realID=togo(int(c));
    if ~isempty(DGRates{ran(c),4,ii})
        figure
        plot(DGRates{ran(c),4,ii})
        title('random')
    end
end


%% Upsilon: plot some stuff
% close all
clear

path1='e:\all_species_temp\newAttempt3';

load(fullfile(path1,'DG_PCanalysis','clust_Kmeans_PC_F1new_corr'))
load(fullfile(path1,'DG','h_normPeaks'))
load(fullfile('e:\all_species_temp','polarity'))
load(fullfile('e:\all_species_temp','h_vel6_new_Black'))
load(fullfile('e:\all_species_temp','h_chirp'))
load(fullfile('e:\all_species_temp','h_flash60sigma'))
% load(fullfile('e:\all_species_temp\newAttempt\DG','h_DGrates'))

dataH=find(tID==3);
dataH=data(dataH);

tmp=find(type==3);
IDh=flashID(tmp);
pol=onoff(tmp);


% all ratios

rat=normPeaks(:,[10 11 12 14 15 16 18 19 20],3)./normPeaks(:,[10 11 12 14 15 16 18 19 20],2);

int=[];
for r=1:9
    tmp=find(rat(:,r)>=1);
    int=[int;tmp];
end
int=unique(int);
clust=dataH(int);
[clust id]=sort(clust);
int=int(id);


% F1=normPeaks(int,[3 7 11 15 19 23],2);
% F2=normPeaks(int,[3 7 11 15 19 23],3);
% F1=normPeaks(int,[3 7 11 15 19 23]-1,2);
% F2=normPeaks(int,[3 7 11 15 19 23]-1,3);

% list=[10 11 12 14 15 16 18 19 20];
list=9:21;

equi2=[1 2 3 4 1 2 3 4 1 2 3 4];
rat=normPeaks(int,list,3)./normPeaks(int,list,2);


F1=[]; F2=[];
for f=1:length(int)
    for r=1:size(rat,2)
        tmp=find(rat(f,r)>=1);
        if ~isempty(tmp)
            
            F1=[F1; normPeaks(int(f),[1:4:21]+equi2(r)-1,2)];
            F2=[F2; normPeaks(int(f),[1:4:21]+equi2(r)-1,3)];
        end
    end
end


F1=[]; F2=[];
for f=1:length(int)
    
    tmp=find(rat(f,:)>=1);
    if ~isempty(tmp)
        curr1=[]; curr2=[];
        for s=1:6
            now=[1 2 3 4]+(4*(s-1));
            
            
            [t1 t2]=max(normPeaks(int(f),now,3));
            if t1==0
                [t1 t2]=max(normPeaks(int(f),now,2));
            end
            
            curr1=[curr1 normPeaks(int(f),now(t2),2)];
            curr2=[curr2 normPeaks(int(f),now(t2),3)];
        end
        F1=[F1;curr1];
        F2=[F2; curr2];
    end
    
end



for f=1:29
    figure
    semilogx([100 200 500 1000 2000 4000],F1(f,:),'-k')
    hold on
    plot([100 200 500 1000 2000 4000],F2(f,:),'-r')
    title(int2str(togo(int(f))))
    axis([50 5000 0 1])
end

figure
curr=normPeaks(int,:,2);
curr=mean(curr,1);
curr=fliplr(reshape(curr,4,6));
[y,x] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr');
colormap('gray')

figure
curr=normPeaks(int,:,3);
curr=mean(curr,1);
curr=fliplr(reshape(curr,4,6));
[y,x] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr');
colormap('gray')



for i=1:length(int)
    figure
    curr=normPeaks(int(i),:,2);
    curr=fliplr(reshape(curr,4,6));
    [y,x] = meshgrid(1:6,1:4);
    imagesc(x(:)',y(:)',curr');
    colormap('gray')
end
%% characterize distance invariant cells (cluster 13 and 14)
close all
clear

path1='F:\all_species_temp\newAttempt3';

load(fullfile('F:\all_species_temp','temporalDG_maxPerSpeed'))

load(fullfile(path1,'DG_PCanalysis','clust_Kmeans_PC_F1new_corr'))
load(fullfile(path1,'DG','h_normPeaks'))
load(fullfile('F:\all_species_temp','polarity'))
load(fullfile('F:\all_species_temp','h_vel6_new_Black'))
load(fullfile('F:\all_species_temp','h_chirpFFT'))
chirpTogo=togo;
load(fullfile('F:\all_species_temp','h_info'))

load(fullfile('F:\all_species_temp','h_flash60sigma'))
load(fullfile('F:\all_species_temp\newAttempt\DG','h_DGrates'))

% human
dataH=find(tID==3);
dataH=data(dataH);

tmp=find(type==3);


IDh=flashID(tmp);
pol=onoff(tmp);

% find distance-invariant cells
tmp=find(dataH==13);
tmp2=find(dataH==14);
DistInv=[tmp;tmp2];
clu=[repmat(13,length(tmp),1);repmat(14,length(tmp2),1)];

%flash
DIpol=[];
for c=1:length(DistInv)
    realID=togo(DistInv(c));
    flash=find(IDh==realID);
    if ~isempty(flash)
        DIpol=[DIpol;pol(flash)];
        figure
        plot(flash_ratesSmoo2{realID})
        title(int2str(clu(c)))
    else
        DIpol=[DIpol;0];
    end
end



%bar
for c=1:length(DistInv)
    realID=togo(DistInv(c));
    if ~isempty(vel6_newB{realID})
        figure
        plot(vel6_newB{realID})
        title(int2str(clu(c)))
    end
end



%chirp
figure
for c=1:length(DistInv)
    realID=togo(DistInv(c));
    tmp=find(chirpTogo==realID);
    if ~isempty(tmp)
        plot(chirpFFTsmooth{tmp,6,4},'color',cols(clu(c)-12));
        hold on
    end
end


% speed tuning
speed=colldata(:,3);
speeds=zeros(length(speed{1,1}),5);

for s=1:5
    speeds(:,s)=speed{s};
end

figure
currS=speeds(DistInv,:);
for s=1:size(currS,1)
    plot([1 2 4 8 16],currS(s,:))
    hold on
end
plot([1 2 4 8 16],mean(currS),'-k','linewidth',2)

rest=setdiff(1:length(speeds),DistInv);
currS=speeds(rest,:);

plot([1 2 4 8 16],mean(currS),'-c','linewidth',2)


%% characterize distance invariant cells (cluster 13 and 14), PIG
close all
clear

path1='F:\all_species_temp\newAttempt3';

load(fullfile(path1,'DG_PCanalysis','clust_Kmeans_PC_F1new_corr'))
load(fullfile(path1,'DG','p_normPeaks'))
load(fullfile('F:\all_species_temp','polarity'))
load(fullfile('F:\all_species_temp','p_vel6_new_Black'))
load(fullfile('F:\all_species_temp','p_chirpFFT'))
chirpTogo=togo;
load(fullfile('F:\all_species_temp','p_info'))
load(fullfile('F:\all_species_temp','p_flash60sigma'))
load(fullfile('F:\all_species_temp\newAttempt\DG','p_DGrates'))

% pig
dataH=find(tID==2);
dataH=data(dataH);
cols='br';

tmp=find(type==2);


IDh=flashID(tmp);
pol=onoff(tmp);

% find distance-invariant cells
tmp=find(dataH==13);
tmp2=find(dataH==14);
DistInv=[tmp;tmp2];
clu=[repmat(13,length(tmp),1);repmat(14,length(tmp2),1)];

%flash
DIpol=[];
for c=1:length(DistInv)
    realID=togo(DistInv(c));
    flash=find(IDh==realID);
    if ~isempty(flash)
        DIpol=[DIpol;pol(flash)];
        figure
        plot(flash_ratesSmoo2{realID})
        title(int2str(clu(c)))
    else
        DIpol=[DIpol;0];
    end
end



%bar
for c=1:length(DistInv)
    realID=togo(DistInv(c));
    if ~isempty(vel6_newB{realID})
        figure
        plot(vel6_newB{realID})
        title(int2str(clu(c)))
    end
end



%chirp
figure
for c=1:length(DistInv)
    realID=togo(DistInv(c));
    tmp=find(chirpTogo==realID)
    if ~isempty(tmp)
        plot(chirpFFTsmooth{tmp,6,4},'color',cols(clu(c)-12));
        hold on
    end
end
%% full-field in different clusters
close all
clear

path1='F:\all_species_temp\newAttempt3';

load(fullfile(path1,'DG_PCanalysis','clust_Kmeans_PC_F1new_corr'))
load(fullfile(path1,'DG','h_normPeaks'))

load(fullfile('F:\all_species_temp','polarity'))


% human
dataH=find(tID==3);
dataH=data(dataH);

tmp=find(type==3);


IDh=flashID(tmp);
pol=onoff(tmp);

clu=unique(dataH);
FlashPerClustH=zeros(14,2);
for c=1:length(clu)
    cc=clu(c);
    tst=find(dataH==cc);
    FlashPerClustH(cc,1)=length(tst);
    pol1=[];
    for t=1:length(tst)
        realID=togo(tst(t));
        tmp=find(IDh==realID);
        if ~isempty(tmp)
            pol1=[pol1;pol(tmp)];
        end
    end
    FlashPerClustH(cc,2)=length(pol1);
end

FlashPerClustH(:,3)=round(100./FlashPerClustH(:,1).*FlashPerClustH(:,2));


% pig
load(fullfile(path1,'DG','p_normPeaks'))

dataH=find(tID==2);
dataH=data(dataH);

tmp=find(type==2);


IDh=flashID(tmp);
pol=onoff(tmp);

clu=unique(dataH);
FlashPerClustP=zeros(14,2);
for c=1:length(clu)
    cc=clu(c);
    tst=find(dataH==cc);
    FlashPerClustP(cc,1)=length(tst);
    pol1=[];
    for t=1:length(tst)
        realID=togo(tst(t));
        tmp=find(IDh==realID);
        if ~isempty(tmp)
            pol1=[pol1;pol(tmp)];
        end
    end
    FlashPerClustP(cc,2)=length(pol1);
end

FlashPerClustP(:,3)=round(100./FlashPerClustP(:,1).*FlashPerClustP(:,2));




% mouse
load(fullfile(path1,'DG','w_normPeaks'))

dataH=find(tID==1);
dataH=data(dataH);

tmp=find(type==1);


IDh=flashID(tmp);
pol=onoff(tmp);

clu=unique(dataH);
FlashPerClustW=zeros(14,2);
for c=1:length(clu)
    cc=clu(c);
    tst=find(dataH==cc);
    FlashPerClustW(cc,1)=length(tst);
    pol1=[];
    for t=1:length(tst)
        realID=togo(tst(t));
        tmp=find(IDh==realID);
        if ~isempty(tmp)
            pol1=[pol1;pol(tmp)];
        end
    end
    FlashPerClustW(cc,2)=length(pol1);
end

FlashPerClustW(:,3)=round(100./FlashPerClustW(:,1).*FlashPerClustW(:,2));

%% speed tuning vs latency: bar
close all
clear

load(fullfile('F:\all_species_temp','latency'))
load(fullfile('F:\all_species_temp','p_info'))
pInfo=info;
load(fullfile('F:\all_species_temp','w_info'))
wInfo=info;
load(fullfile('F:\all_species_temp','h_info'))
hInfo=info;

%pig
pblack=[]; pblackID=[];
pwhite=[]; pwhiteID=[];
for p=1:length(pInfo)
    if ~isempty(pInfo{p,16}) && pInfo{p,16}>0
        pblack=[pblack;pInfo{p,16}];
        pblackID=[pblackID;p];
    end
    if ~isempty(pInfo{p,17}) && pInfo{p,17}>0
        pwhite=[pwhite;pInfo{p,17}];
        pwhiteID=[pwhiteID;p];
    end
end

tmp=find(type==2);
pLat=allLat(tmp,:);
platID=flashID(tmp);

pbID=[]; pblID=[];
for b=1:length(pblackID)
    if ~isempty(find(platID==pblackID(b)))
        pbID=[pbID;b];
        pblID=[pblID;find(platID==pblackID(b))];
    end
end
pwID=[]; pwlID=[];
for b=1:length(pwhiteID)
    if ~isempty(find(platID==pwhiteID(b)))
        pwID=[pwID;b];
        pwlID=[pwlID;find(platID==pwhiteID(b))];
    end
end

figure(1)
subplot(2,2,1)
scatter(pblack(pbID),pLat(pblID,1),'ok','markerfacecolor','k')
hold on
scatter(pblack(pbID),pLat(pblID,2),'ok')
scatter(pblack(pbID),pLat(pblID,3),'dk')
scatter(pblack(pbID),pLat(pblID,4),'dk','markerfacecolor','k')
scatter(pwhite(pwID),pLat(pwlID,1),'or','markerfacecolor','r')
scatter(pwhite(pwID),pLat(pwlID,2),'or')
scatter(pwhite(pwID),pLat(pwlID,3),'dr')
scatter(pwhite(pwID),pLat(pwlID,4),'dr','markerfacecolor','r')


tmp=unique([pwlID;pblID]);
figure(2)
subplot(2,2,1)
scatter(pLat(tmp,1),pLat(tmp,2),'.r')
hold on
scatter(pLat(tmp,4),pLat(tmp,3),'.r')
plot([0 1800],[0 1800],'-k')
xlabel('neg. contrast')
ylabel('pos. contrast')


figure(3)
subplot(2,2,1)
rat1b=pLat(pblID,1)./pLat(pblID,2);
rat2b=pLat(pblID,4)./pLat(pblID,3);
rat1w=pLat(pwlID,1)./pLat(pwlID,2);
rat2w=pLat(pwlID,4)./pLat(pwlID,3);
plot(pblack(pbID),rat1b,'o')
hold on
plot(pblack(pbID),rat2b,'or')
plot(pwhite(pwID),rat1w,'o')
plot(pwhite(pwID),rat2w,'or')
plot([0 12],[1 1],'-k')
xlabel('speed pref')
ylabel('latBlack/latWhite')
medB1=[]; medW1=[];
medB2=[]; medW2=[];
for s=0:0.5:11.5
    tmp=find(pblack(pbID)>=s);
    
    tmp2=find(pblack(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1b(ii)>0);
    tst2=ii(tst);
    tst=find(rat1b(tst2)~=inf);
    tst2=tst2(tst);
    medB1=[medB1;median(rat1b(tst2))];
    tst=find(rat2b(ii)>0);
    tst2=ii(tst);
    tst=find(rat2b(tst2)~=inf);
    tst2=tst2(tst);
    medB2=[medB2;median(rat2b(tst2))];
    
    tmp=find(pwhite(pwID)>=s);
    
    tmp2=find(pwhite(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1w(ii)>0);
    tst2=ii(tst);
    tst=find(rat1w(tst2)~=inf);
    tst2=tst2(tst);
    medW1=[medW1;median(rat1w(tst2))];
    tst=find(rat2w(ii)>0);
    tst2=ii(tst);
    tst=find(rat2w(tst2)~=inf);
    tst2=tst2(tst);
    medW2=[medW2;median(rat2w(tst2))];
end
for s=1:length(medB1)
    if ~isnan(medB1(s))
        if medB1(s)<1
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medB2(s))
        if medB2(s)<1
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'--m','linewidth',2)
        end
    end
    if ~isnan(medW1(s))
        if medW1(s)<1
            plot([s/2-0.5 s/2],[medW1(s) medW1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medW1(s) medW1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medW2(s))
        if medW2(s)<1
            plot([s/2-0.5 s/2],[medW2(s) medW2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medW2(s) medW2(s)],'--m','linewidth',2)
        end
    end
end
allrat=[rat1b;rat2b;rat1w;rat2w]; tmp=find(allrat>0); tmp2=find(allrat(tmp)~=inf);
plot([0 12],[median(allrat(tmp(tmp2))) median(allrat(tmp(tmp2)))],'-g','linewidth',2)


% mouse
wblack=[]; wblackID=[];
wwhite=[]; wwhiteID=[];
for p=1:length(wInfo)
    if ~isempty(wInfo{p,16}) && wInfo{p,16}>0
        wblack=[wblack;wInfo{p,16}];
        wblackID=[wblackID;p];
    end
    if ~isempty(wInfo{p,17}) && wInfo{p,17}>0
        wwhite=[wwhite;wInfo{p,17}];
        wwhiteID=[wwhiteID;p];
    end
end

tmp=find(type==1);
wLat=allLat(tmp,:);
wlatID=flashID(tmp);


wbID=[]; wblID=[];
for b=1:length(wblackID)
    if ~isempty(find(wlatID==wblackID(b)))
        wbID=[wbID;b];
        wblID=[wblID;find(wlatID==wblackID(b))];
    end
end
wwID=[]; wwlID=[];
for b=1:length(wwhiteID)
    if ~isempty(find(wlatID==wwhiteID(b)))
        wwID=[wwID;b];
        wwlID=[wwlID;find(wlatID==wwhiteID(b))];
    end
end

figure(1)
subplot(2,2,2)
scatter(wblack(wbID),wLat(wblID,1),'ok','markerfacecolor','k')
hold on
scatter(wblack(wbID),wLat(wblID,2),'ok')
scatter(wblack(wbID),wLat(wblID,3),'dk')
scatter(wblack(wbID),wLat(wblID,4),'dk','markerfacecolor','k')
scatter(wwhite(wwID),wLat(wwlID,1),'or','markerfacecolor','r')
scatter(wwhite(wwID),wLat(wwlID,2),'or')
scatter(wwhite(wwID),wLat(wwlID,3),'dr')
scatter(wwhite(wwID),wLat(wwlID,4),'dr','markerfacecolor','r')



tmp=unique([wwlID;wblID]);
figure(2)
subplot(2,2,2)
scatter(wLat(tmp,1),wLat(tmp,2),'.r')
hold on
scatter(wLat(tmp,4),wLat(tmp,3),'.r')
plot([0 1800],[0 1800],'-k')
xlabel('neg. contrast')
ylabel('pos. contrast')



figure(3)
subplot(2,2,2)
rat1b=wLat(wblID,1)./wLat(wblID,2);
rat2b=wLat(wblID,4)./wLat(wblID,3);
rat1w=wLat(wwlID,1)./wLat(wwlID,2);
rat2w=wLat(wwlID,4)./wLat(wwlID,3);
plot(wblack(wbID),rat1b,'o')
hold on
plot(wblack(wbID),rat2b,'or')
plot(wwhite(wwID),rat1w,'o')
plot(wwhite(wwID),rat2w,'or')
plot([0 12],[1 1],'-k')
xlabel('speed pref')
ylabel('latBlack/latWhite')
medB1=[]; medW1=[];
medB2=[]; medW2=[];
for s=0:0.5:11.5
    tmp=find(wblack(wbID)>=s);
    
    tmp2=find(wblack(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1b(ii)>0);
    tst2=ii(tst);
    tst=find(rat1b(tst2)~=inf);
    tst2=tst2(tst);
    medB1=[medB1;median(rat1b(tst2))];
    tst=find(rat2b(ii)>0);
    tst2=ii(tst);
    tst=find(rat2b(tst2)~=inf);
    tst2=tst2(tst);
    medB2=[medB2;median(rat2b(tst2))];
    
    tmp=find(wwhite(wwID)>=s);
    
    tmp2=find(wwhite(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1w(ii)>0);
    tst2=ii(tst);
    tst=find(rat1w(tst2)~=inf);
    tst2=tst2(tst);
    medW1=[medW1;median(rat1w(tst2))];
    tst=find(rat2w(ii)>0);
    tst2=ii(tst);
    tst=find(rat2w(tst2)~=inf);
    tst2=tst2(tst);
    medW2=[medW2;median(rat2w(tst2))];
end
for s=1:length(medB1)
    if ~isnan(medB1(s))
        if medB1(s)<1
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medB2(s))
        if medB2(s)<1
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'--m','linewidth',2)
        end
    end
    if ~isnan(medW1(s))
        if medW1(s)<1
            plot([s/2-0.5 s/2],[medW1(s) medW1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medW1(s) medW1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medW2(s))
        if medW2(s)<1
            plot([s/2-0.5 s/2],[medW2(s) medW2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medW2(s) medW2(s)],'--m','linewidth',2)
        end
    end
end
allrat=[rat1b;rat2b;rat1w;rat2w]; tmp=find(allrat>0); tmp2=find(allrat(tmp)~=inf);
plot([0 12],[median(allrat(tmp(tmp2))) median(allrat(tmp(tmp2)))],'-g','linewidth',2)




% human
hblack=[]; hblackID=[];
hwhite=[]; hwhiteID=[];
for p=1:length(hInfo)
    if ~isempty(hInfo{p,16}) && hInfo{p,16}>0
        hblack=[hblack;hInfo{p,16}];
        hblackID=[hblackID;p];
    end
    if ~isempty(hInfo{p,17}) && hInfo{p,17}>0
        hwhite=[hwhite;hInfo{p,17}];
        hwhiteID=[hwhiteID;p];
    end
end

tmp=find(type==3);
hLat=allLat(tmp,:);
hlatID=flashID(tmp);


hbID=[]; hblID=[];
for b=1:length(hblackID)
    if ~isempty(find(hlatID==hblackID(b)))
        hbID=[hbID;b];
        hblID=[hblID;find(hlatID==hblackID(b))];
    end
end
hwID=[]; hwlID=[];
for b=1:length(hwhiteID)
    if ~isempty(find(hlatID==hwhiteID(b)))
        hwID=[hwID;b];
        hwlID=[hwlID;find(hlatID==hwhiteID(b))];
    end
end

figure(1)
subplot(2,2,3)
scatter(hblack(hbID),hLat(hblID,1),'ok','markerfacecolor','k')
hold on
scatter(hblack(hbID),hLat(hblID,2),'ok')
scatter(hblack(hbID),hLat(hblID,3),'dk')
scatter(hblack(hbID),hLat(hblID,4),'dk','markerfacecolor','k')
scatter(hwhite(hwID),hLat(hwlID,1),'or','markerfacecolor','r')
scatter(hwhite(hwID),hLat(hwlID,2),'or')
scatter(hwhite(hwID),hLat(hwlID,3),'dr')
scatter(hwhite(hwID),hLat(hwlID,4),'dr','markerfacecolor','r')



tmp=unique([hwlID;hblID]);
figure(2)
subplot(2,2,3)
scatter(hLat(tmp,1),hLat(tmp,2),'.r')
hold on
scatter(hLat(tmp,4),hLat(tmp,3),'.r')
plot([0 1800],[0 1800],'-k')
xlabel('neg. contrast')
ylabel('pos. contrast')



figure(3)
subplot(2,2,3)
rat1b=hLat(hblID,1)./hLat(hblID,2); %black bar neg contrast
rat2b=hLat(hblID,4)./hLat(hblID,3);%black bar pos contrast
rat1w=hLat(hwlID,1)./hLat(hwlID,2);%white bar neg contrast
rat2w=hLat(hwlID,4)./hLat(hwlID,3);%white bar pos contrast
plot(hblack(hbID),rat1b,'o')
hold on
plot(hblack(hbID),rat2b,'or')
plot(hwhite(hwID),rat1w,'o')
plot(hwhite(hwID),rat2w,'or')
plot([0 12],[1 1],'-k')
xlabel('speed pref')
ylabel('latBlack/latWhite')
medB1=[]; medW1=[];
medB2=[]; medW2=[];
for s=0:0.5:11.5
    tmp=find(hblack(hbID)>=s);
    
    tmp2=find(hblack(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1b(ii)>0);
    tst2=ii(tst);
    tst=find(rat1b(tst2)~=inf);
    tst2=tst2(tst);
    medB1=[medB1;median(rat1b(tst2))];
    tst=find(rat2b(ii)>0);
    tst2=ii(tst);
    tst=find(rat2b(tst2)~=inf);
    tst2=tst2(tst);
    medB2=[medB2;median(rat2b(tst2))];
    
    tmp=find(hwhite(hwID)>=s);
    
    tmp2=find(hwhite(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1w(ii)>0);
    tst2=ii(tst);
    tst=find(rat1w(tst2)~=inf);
    tst2=tst2(tst);
    medW1=[medW1;median(rat1w(tst2))];
    tst=find(rat2w(ii)>0);
    tst2=ii(tst);
    tst=find(rat2w(tst2)~=inf);
    tst2=tst2(tst);
    medW2=[medW2;median(rat2w(tst2))];
end
for s=1:length(medB1)
    if ~isnan(medB1(s))
        if medB1(s)<1
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medB2(s))
        if medB2(s)<1
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'--m','linewidth',2)
        end
    end
    if ~isnan(medW1(s))
        if medW1(s)<1
            plot([s/2-0.5 s/2],[medW1(s) medW1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medW1(s) medW1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medW2(s))
        if medW2(s)<1
            plot([s/2-0.5 s/2],[medW2(s) medW2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medW2(s) medW2(s)],'--m','linewidth',2)
        end
    end
end
allrat=[rat1b;rat2b;rat1w;rat2w]; tmp=find(allrat>0); tmp2=find(allrat(tmp)~=inf);
plot([0 12],[median(allrat(tmp(tmp2))) median(allrat(tmp(tmp2)))],'-g','linewidth',2)

for s=1:3
    subplot(2,2,s)
    axis([0 12 0.4 8])
end
%% speed tuning vs latency: DG
close all
clear

load(fullfile('F:\all_species_temp','latency'))
load(fullfile('F:\all_species_temp','p_info'))
load(fullfile('F:\all_species_temp','p_DGspeed'))
pInfo=info; pSpeed=DGspeed; ptogo=togo;
load(fullfile('F:\all_species_temp','p_info'))
load(fullfile('F:\all_species_temp','p_DGspeed'))
wInfo=info; wSpeed=DGspeed; wtogo=togo;
load(fullfile('F:\all_species_temp','h_info'))
load(fullfile('F:\all_species_temp','p_DGspeed'))
hInfo=info; hSpeed=DGspeed; htogo=togo;


%pig

tmp=find(type==2);
pLat=allLat(tmp,:);
platID=flashID(tmp);

pbID=[]; pblID=[];
for b=1:length(pSpeed)
    if ~isempty(find(platID==ptogo(b)))
        pbID=[pbID;b];
        pblID=[pblID;find(platID==ptogo(b))];
    end
end


figure(1)
subplot(2,2,1)
scatter(pSpeed(pbID),pLat(pblID,1),'ok','markerfacecolor','k')
hold on
scatter(pSpeed(pbID),pLat(pblID,2),'ok')
scatter(pSpeed(pbID),pLat(pblID,3),'dk')
scatter(pSpeed(pbID),pLat(pblID,4),'dk','markerfacecolor','k')



figure(2)
subplot(2,2,1)
scatter(pLat(pblID,1),pLat(pblID,2),'or')
hold on
scatter(pLat(pblID,4),pLat(pblID,3),'dr')
plot([0 1800],[0 1800],'-k')
xlabel('neg. contrast')
ylabel('pos. contrast')


figure(3)
subplot(2,2,1)
rat1b=pLat(pblID,1)./pLat(pblID,2);
rat2b=pLat(pblID,4)./pLat(pblID,3);

plot(pSpeed(pbID),rat1b,'o')
hold on
plot(pSpeed(pbID),rat2b,'or')
plot([0 12],[1 1],'-k')
xlabel('speed pref')
ylabel('latBlack/latWhite')
medB1=[];
medB2=[];
for s=0:0.5:11.5
    tmp=find(pSpeed(pbID)>=s);
    
    tmp2=find(pSpeed(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1b(ii)>0);
    tst2=ii(tst);
    tst=find(rat1b(tst2)~=inf);
    tst2=tst2(tst);
    medB1=[medB1;median(rat1b(tst2))];
    tst=find(rat2b(ii)>0);
    tst2=ii(tst);
    tst=find(rat2b(tst2)~=inf);
    tst2=tst2(tst);
    medB2=[medB2;median(rat2b(tst2))];
    
end
for s=1:length(medB1)
    if ~isnan(medB1(s))
        if medB1(s)<1
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medB2(s))
        if medB2(s)<1
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'--m','linewidth',2)
        end
    end
    
end
allrat=[rat1b;rat2b]; tmp=find(allrat>0); tmp2=find(allrat(tmp)~=inf);
plot([0 12],[median(allrat(tmp(tmp2))) median(allrat(tmp(tmp2)))],'-g','linewidth',2)


% mouse

tmp=find(type==1);
wLat=allLat(tmp,:);
wlatID=flashID(tmp);


wbID=[]; wblID=[];
for b=1:length(wSpeed)
    if ~isempty(find(wlatID==wtogo(b)))
        wbID=[wbID;b];
        wblID=[wblID;find(wlatID==wtogo(b))];
    end
end


figure(1)
subplot(2,2,2)
scatter(wSpeed(wbID),wLat(wblID,1),'ok','markerfacecolor','k')
hold on
scatter(wSpeed(wbID),wLat(wblID,2),'ok')
scatter(wSpeed(wbID),wLat(wblID,3),'dk')
scatter(wSpeed(wbID),wLat(wblID,4),'dk','markerfacecolor','k')



figure(2)
subplot(2,2,2)
scatter(wLat(wblID,1),wLat(wblID,2),'or')
hold on
scatter(wLat(wblID,4),wLat(wblID,3),'dr')
plot([0 1800],[0 1800],'-k')
xlabel('neg. contrast')
ylabel('pos. contrast')


figure(3)
subplot(2,2,2)
rat1b=wLat(wblID,1)./wLat(wblID,2);
rat2b=wLat(wblID,4)./wLat(wblID,3);

plot(wSpeed(wbID),rat1b,'o')
hold on
plot(wSpeed(wbID),rat2b,'or')
plot([0 12],[1 1],'-k')
xlabel('speed pref')
ylabel('latBlack/latWhite')
medB1=[];
medB2=[];
for s=0:0.5:11.5
    tmp=find(wSpeed(wbID)>=s);
    
    tmp2=find(wSpeed(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1b(ii)>0);
    tst2=ii(tst);
    tst=find(rat1b(tst2)~=inf);
    tst2=tst2(tst);
    medB1=[medB1;median(rat1b(tst2))];
    tst=find(rat2b(ii)>0);
    tst2=ii(tst);
    tst=find(rat2b(tst2)~=inf);
    tst2=tst2(tst);
    medB2=[medB2;median(rat2b(tst2))];
    
end
for s=1:length(medB1)
    if ~isnan(medB1(s))
        if medB1(s)<1
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medB2(s))
        if medB2(s)<1
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'--m','linewidth',2)
        end
    end
    
end
allrat=[rat1b;rat2b]; tmp=find(allrat>0); tmp2=find(allrat(tmp)~=inf);
plot([0 12],[median(allrat(tmp(tmp2))) median(allrat(tmp(tmp2)))],'-g','linewidth',2)




% human


tmp=find(type==3);
hLat=allLat(tmp,:);
hlatID=flashID(tmp);


hbID=[]; hblID=[];
for b=1:length(hSpeed)
    if ~isempty(find(hlatID==htogo(b)))
        hbID=[hbID;b];
        hblID=[hblID;find(hlatID==htogo(b))];
    end
end


figure(1)
subplot(2,2,3)
scatter(hSpeed(hbID),hLat(hblID,1),'ok','markerfacecolor','k')
hold on
scatter(hSpeed(hbID),hLat(hblID,2),'ok')
scatter(hSpeed(hbID),hLat(hblID,3),'dk')
scatter(hSpeed(hbID),hLat(hblID,4),'dk','markerfacecolor','k')



figure(2)
subplot(2,2,3)
scatter(hLat(hblID,1),hLat(hblID,2),'or')
hold on
scatter(hLat(hblID,4),hLat(hblID,3),'dr')
plot([0 1800],[0 1800],'-k')
xlabel('neg. contrast')
ylabel('pos. contrast')


figure(3)
subplot(2,2,3)
rat1b=hLat(hblID,1)./hLat(hblID,2);
rat2b=hLat(hblID,4)./hLat(hblID,3);

plot(hSpeed(hbID),rat1b,'o')
hold on
plot(hSpeed(hbID),rat2b,'or')
plot([0 12],[1 1],'-k')
xlabel('speed pref')
ylabel('latBlack/latWhite')
medB1=[];
medB2=[];
for s=0:0.5:11.5
    tmp=find(hSpeed(hbID)>=s);
    
    tmp2=find(hSpeed(tmp)<s+0.5);
    ii=tmp(tmp2);
    tst=find(rat1b(ii)>0);
    tst2=ii(tst);
    tst=find(rat1b(tst2)~=inf);
    tst2=tst2(tst);
    medB1=[medB1;median(rat1b(tst2))];
    tst=find(rat2b(ii)>0);
    tst2=ii(tst);
    tst=find(rat2b(tst2)~=inf);
    tst2=tst2(tst);
    medB2=[medB2;median(rat2b(tst2))];
    
end
for s=1:length(medB1)
    if ~isnan(medB1(s))
        if medB1(s)<1
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'-k','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB1(s) medB1(s)],'--k','linewidth',2)
        end
    end
    if ~isnan(medB2(s))
        if medB2(s)<1
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'-m','linewidth',2)
        else
            plot([s/2-0.5 s/2],[medB2(s) medB2(s)],'--m','linewidth',2)
        end
    end
    
end
allrat=[rat1b;rat2b]; tmp=find(allrat>0); tmp2=find(allrat(tmp)~=inf);
plot([0 12],[median(allrat(tmp(tmp2))) median(allrat(tmp(tmp2)))],'-g','linewidth',2)


%% polarity of pig DS-cells
clear
close all

path1='F:\all_species_temp';
load(fullfile(path1,'polarity'))
load(fullfile(path1,'p_info'))


DS=[];
for i=1:length(info)
    if info{i,12}==1
        if goodones(i)==1
            DS=[DS;i];
        end
    end
end


pig=find(type==2);
flash=flashID(pig);
pol=onoff(pig);

polDS=[];
for d=1:length(DS)
    tmp=find(flash==DS(d));
    
    if ~isempty(tmp)
        polDS=[polDS;pol(tmp)];
    else
        polDS=[polDS;NaN];
    end
end

onoff=100/length(polDS)*length(find(polDS==2))
on=100/length(polDS)*length(find(polDS==1))
off=100/length(polDS)*length(find(polDS==-1))
none=100/length(polDS)*length(find(isnan(polDS)))
%% LF of ON-OFF cells
clear
close all


path1='F:\all_species_temp';

ty='wph';
load(fullfile(path1,'polarity'))

LFonoffW=[]; LFonoffP=[]; LFonoffH=[];
for t=1:3
    load(fullfile(path1,[ty(t),'_info']))
    
    id=find(type==t);
    fl=flashID(id);
    pol=onoff(id);
    
    ONOFF=find(pol==2);
    onoffID=fl(ONOFF);
    
    for o=1:length(onoffID)
        lf=info{onoffID(o),9};
        if goodones(onoffID(o))==1
            if t==1
                LFonoffW=[LFonoffW;lf];
            elseif t==2
                LFonoffP=[LFonoffP;lf];
            else
                LFonoffH=[LFonoffH;lf];
            end
        end
    end
end



dominant=zeros(3,2);
pLF=[-1 1];

for t=1:3
    if t==1
        data=LFonoffW;
    elseif t==2
        data=LFonoffP;
    else
        data=LFonoffH;
    end
    cnt=0;
    for i=pLF
        cnt=cnt+1;
        tmp=find(data==i);
        dominant(t,cnt)=round(100/length(data)*length(tmp));
    end
    
    
end
%% PER EXP Latency
clear
close all

path1='F:\all_species_temp';

dtype='wph';

col=colormap(hsv);
close all
col=col(1:4:end,:);

load(fullfile(path1,'latency'))
for d=1:3
    figure
    load(fullfile(path1,[dtype(d),'_info']))
    
    currID=find(type==d);
    currFlash=flashID(currID);
    currLat=allLat(currID,:);
    
    dates=info(currFlash,1);
    udates=unique(dates);
    
    leg=cell(1,1)
    for u=1:length(udates)
        list=[];
        for t=1:length(dates)
            if ~isempty(regexp(dates{t},udates{u}))
                list=[list;t];
            end
        end
        
        
        now=reshape(currLat(list,:),length(list)*4,1);
        tmp=find(now>0);
        now=now(tmp);
        med=median(now);
        st=std(now);
        hold on
        plot(1, med,'o','color',col(u,:),'markersize',5)
        leg{u}=[udates{u},' / ',num2str(med)];
    end
    
    allmed=reshape(currLat,length(currLat)*4,1);
    tmp=find(allmed>0);
    allmed=median(allmed(tmp));
    leg{u+1}=num2str(allmed);
    
    plot(1,allmed,'ok','markerfacecolor','k')
    legend(leg)
    title(dtype(d))
    axis([0.5 1.5 100 500])
    
    figure
    
    
    
    for u=1:length(udates)
        list=[];
        for t=1:length(dates)
            if ~isempty(regexp(dates{t},udates{u}))
                list=[list;t];
            end
        end
        
        subplot(4,4,u)
        now=reshape(currLat(list,:),length(list)*4,1);
        tmp=find(now>0);
        now=now(tmp);
        med=median(now);
        st=std(now);
        now=[now;0;2000];
        [h,x]=hist(now,40);
        hist(now,40)
        hold on
        plot(med,max(h)+1,'or','markerfacecolor','r')
        plot([med-st med+st],[max(h)+1 max(h)+1],'-r')
        allmed=reshape(currLat,length(currLat)*4,1);
        tmp=find(allmed>0);
        allmed=median(allmed(tmp));
        plot(allmed,max(h)+2,'og','markerfacecolor','g')
        title(udates{u})
        axis([0 2000 0 max(h)+4])
        
        
    end
    %     set(gcf,'position',[1 31 1600 794])
    %     saveas(gcf,fullfile('F:\all_species_temp\final_pics',['_latency_',dtype(d),'.emf']))
end


%% PER EXP Polarity


clear
close all

path1='F:\all_species_temp';

dtype='wph';

load(fullfile(path1,'polarity'))
for d=1:3
    figure
    load(fullfile(path1,[dtype(d),'_info']))
    
    currID=find(type==d);
    currFlash=flashID(currID);
    currPol=onoff(currID);
    
    dates=info(currFlash,1);
    udates=unique(dates);
    
    for u=1:length(udates)
        list=[];
        for t=1:length(dates)
            if ~isempty(regexp(dates{t},udates{u}))
                list=[list;t];
            end
        end
        
        subplot(4,4,u)
        re=[-1 1 2 0];
        now=currPol(list);
        num=[];
        for r=1:length(re)
            tmp=find(now==re(r));
            num=[num;length(tmp)];
        end
        
        
        h=pie(num)
        hText = findobj(h,'Type','text'); % text handles
        percentValues = get(hText,'String'); % percent values
        str = {'off ';'on ';'onoff ';'none '}; % text
        tmp=find(num>0);
        leg=cell(length(tmp),1);
        for tt=1:length(tmp)
            leg{tt}=strcat(str(tmp(tt)),int2str(num(tmp(tt))));
        end
        combinedstrings = strcat(str(tmp),percentValues);
        set(hText,{'String'},combinedstrings);
        title(udates{u})
        
        
    end
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('F:\all_species_temp\final_pics',['_polarity_',dtype(d),'.emf']))
end
%% PER EXP LF


clear
close all

path1='F:\all_species_temp';
badIDs=202:250;
dtype='wph';

for d=1:3
    figure
    load(fullfile(path1,[dtype(d),'_info']))
    
    good=find(goodones==1);
    if d==2
        good=setdiff(good,badIDs);
    end
    LF=cell2mat(info(good,9));
    
    dates=info(good,1);
    udates=unique(dates);
    
    for u=1:length(udates)
        list=[];
        for t=1:length(dates)
            if ~isempty(regexp(dates{t},udates{u}))
                list=[list;t];
            end
        end
        
        subplot(4,4,u)
        re=[-1 1];
        now=LF(list);
        num=[];
        for r=1:length(re)
            tmp=find(now==re(r));
            num=[num;length(tmp)];
        end
        
        if ~isempty(find(num>0))
            
            h=pie(num)
            hText = findobj(h,'Type','text'); % text handles
            percentValues = get(hText,'String'); % percent values
            str = {'off ';'on '}; % text
            tmp=find(num>0);
            leg=cell(length(tmp),1);
            for tt=1:length(tmp)
                leg{tt}=strcat(str(tmp(tt)),int2str(num(tmp(tt))));
            end
            combinedstrings = strcat(str(tmp),percentValues);
            set(hText,{'String'},combinedstrings);
            title(udates{u})
        end
        
    end
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('F:\all_species_temp\final_pics',['_polarityLF_',dtype(d),'.emf']))
end
%% PER EXP speed from DG


clear
close all
superpath='f:\all_species_temp';
cols='bgr';
Dtype='wph';
spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
speed=spat.*freq/1000; speed=speed';
% uspeed=unique(speed);
% uspeed=[0.5 1 2 4 8 16 32];
uspeed=[1 2 4 8 16];
% uspeed=[1 2 4 8];

col=colormap(hsv);
close all
col=col(1:4:end,:);


allmeans=[2.98 3.57 3.99];

ids=[];
for u=1:length(uspeed)
    tmp=find(speed==uspeed(u));
    ids=[ids;tmp];
end

speed=speed(ids);

NofE=zeros(3,1);
alldata=cell(3,1);
for d=1:3
    load(fullfile(superpath,[Dtype(d),'_DGspeed']))
    load(fullfile(superpath,[Dtype(d),'_info']))
    if d==1
        mall=DGspeed;
    elseif d==2
        pall=DGspeed;
    else
        hall=DGspeed;
    end
    
    allvalorig=DGspeed;
    dates=info(togo,1);
    udates2=unique(info(:,1));
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    start=[start;length(allvalorig)];
    udates=udates(ia);
    ends=start(2:end)-1;
    [udates ai]=sort(udates);
    start=start(ai);
    ends=ends(ai);
    NofE(d)=length(udates);
    
    leg=cell(1,1);
    collmeds=[]; ucnt=0;
    for u=1:length(udates2)
        if sum(ismember(udates,udates2{u}))==1
            ucnt=ucnt+1;
            allval=allvalorig(start(ucnt):ends(ucnt));
            tmp=find(allval>0);
            allval=allval(tmp);
            alldata{d,u}=allval;
            med=median(allval);
            
            
            if length(allval)>2
                figure(d)
                c=cdfplot(allval)
                set(c,'color',col(ucnt,:),'linestyle','-','linewidth',2)
                hold on
                leg{ucnt}=[udates{ucnt}, ' / ',num2str(med),' (',int2str(length(allval)),' cells)'];
            end
            axis([0 12 0 1])
            xlabel('pref. speed (mm/s)')
            ylabel('% cells')
            title(Dtype(d))
            
            
            if length(allval)>2
                
                
                figure(4)
                subplot(1,2,1)
                c=cdfplot(allval)
                set(c,'color',cols(d),'linestyle','-','linewidth',2)
                hold on
                subplot(1,2,2)
                plot(d,med,'o','markerfacecolor',cols(d),'markeredgecolor',cols(d))
                hold on
                %             plot([t t],[med-st med+st],'-','color',cols(t))
                %             plot([t t],[med-conf med+conf],'-','color',cols(t),'linewidth',2)
                axis([0 4 0 inf])
                collmeds=[collmeds;med];
            end
        else
            alldata{d,u}=[];
            
        end
    end
    
    
    plot([d-0.2 d+0.2],[mean(collmeds) mean(collmeds)],'-','color',cols(d),'linewidth',2) %mean of means
    plot([d d],[mean(collmeds)-std(collmeds) mean(collmeds)+std(collmeds)],'-','color',cols(d)) %mean of means
    plot([d-0.2 d+0.2],[allmeans(d) allmeans(d)],'-','color',cols(d)) %previously calculated population mean
    subplot(1,2,1)
    plot(mean(collmeds), 0.5,'o','color',cols(d),'markerfacecolor',cols(d),'markersize',10)
    plot(allmeans(d),0.5,'o','color',cols(d),'markersize',10)
    
    figure(d)
    legend(leg)
end



% within species
goodExps=cell(3,1); pall=cell(3,1);

for t=1:3
    load(fullfile(superpath,[Dtype(t),'_info']))
    load(fullfile(superpath,[Dtype(t),'_DGspeed']))
    
    dates=info(:,1);
    udates=unique(dates);
    
    
    if t==2
        gooddates=[];
        for dd=1:length(udates)
            if isempty(regexp(udates{dd},'20140723d'))
                gooddates=[gooddates;dd];
            end
        end
        udates=udates(gooddates);
    end
    alldates=udates;
    
    goodExp=[];
    for e=1:15
        if ~isempty(alldata{t,e})
            goodExp=[goodExp;e];
        end
    end
    udates=udates(goodExp);
    goodExps{t}=goodExp;
    pvals=[];
    for g=1:length(alldates)-1
        for g2=g+1:length(alldates)
            if sum(ismember(udates,alldates{g})==1) && sum(ismember(udates,alldates{g2})==1)
                p=ranksum(alldata{t,g},alldata{t,g2});
            else
                p=-1 ;
            end
            pvals=[pvals;p];
        end
    end
    pall{t}=pvals;
    matP=zeros(length(alldates),length(alldates));
    numOfComp=[0 length(alldates)-1:-1:1];
    for g=1:length(alldates)-1
        if sum(ismember(udates,alldates{g})==1)
            matP(g,g)=1;
            currID1=sum(numOfComp(1:g))+1;
            currID2=sum(numOfComp(1:g+1));
            matP(g,g+1:end)=pvals(currID1:currID2);
            matP(g+1:end,g)=pvals(currID1:currID2);
        else
            matP(g,g)=-1;
            matP(g,g+1:end)=repmat(-1,length(alldates)-g,1);
            matP(g+1:end,g)=repmat(-1,length(alldates)-g,1);
        end
    end
    if sum(ismember(udates,alldates{g+1})==1)
        matP(g+1,g+1)=1;
    else
        matP(g+1,g+1)=-1;
    end
    matP=1-matP;
    tmp0=find(matP==2);
    tmp1=find(matP>0.95);
    tmp2=find(matP>0.99);
    tmp3=find(matP<0.95); matP(tmp3)=0;
    matP(tmp1)=0.8; matP(tmp2)=1; matP(tmp0)=0.4;
    figure
    imagesc(matP,[0 1])
    set(gca,'xticklabel',alldates)
    set(gca,'yticklabel',alldates)
    colormap('gray')
    title(['p-values across ',Dtype(t),' experiments'])
end
figure
subplot(1,2,1)
for t=1:3
    plot(t,pall{t},'o','color',cols(t))
    hold on
    plot([t-0.2 t+0.2],[median(pall{t}) median(pall{t})],'-','color',cols(t))
end

plot([0 4],[0.05 0.05],'-k')
plot([0 4],[0.01 0.01],'-k')
axis([0 4 0 1])
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'m','p','h'})
title('ranksum p-values between all experiments within each species')



% across species
pcross=cell(3,1); cnt=0;
for t=1:2
    for tt=t+1:3
        pvals=[]; cnt=cnt+1;
        for g=1:length(goodExps{t})
            for g2=1:length(goodExps{tt})
                p=ranksum(alldata{t,goodExps{t}(g)},alldata{tt,goodExps{tt}(g2)});
                pvals=[pvals;p];
            end
        end
        pcross{cnt}=pvals;
        malength=max(length(goodExps{t}),length(goodExps{tt}));
        matP=zeros(length(goodExps{t}),length(goodExps{tt}));
        numOfComp=[0 length(goodExp)-1:-1:1];
        
    end
end


subplot(1,2,2)
for t=1:3
    plot(t,pcross{t},'o','color',cols(t))
    hold on
    plot([t-0.2 t+0.2],[median(pcross{t}) median(pcross{t})],'-','color',cols(t))
end

plot([0 4],[0.05 0.05],'-k')
plot([0 4],[0.01 0.01],'-k')
axis([0 4 0 1])
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'m-p','m-h','p-h'})
title('ranksum p-values between all experiments across species')
%% bootstrap DG speed
% result: m-p (0.038), m-h (0.003), p-h (0.381)
% result: m (0.507), p (0.502), h (0.511)
clear
close all
superpath='f:\all_species_temp';

types='wph';
cols='bgr';
titles={'mouse','pig','human'};

speedpref=[]; realIDs=[]; tID=[]; dates=[];
for t=1:length(types)
    
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_DGspeed']))
    
    
    date=info(togo,1);
    
    allvalorig=DGspeed;
    tmp=find(allvalorig>0);
    date=date(tmp);
    dates=[dates;date];
    allvalorig=allvalorig(tmp);
    realIDs=[realIDs;tmp];
    tID=[tID;repmat(t,length(tmp),1)];
    speedpref=[speedpref;allvalorig];
end


dlength=[];
for t=1:3
    tpm=find(tID==t);
    
    dlength=[dlength;length(tpm)];
end

samplesize=round(length(tID)/length(unique(dates)));
% samplesize=15;

% within species
within=zeros(3,1000);
for t=1:3
    datanow=find(tID==t);
    datanow=speedpref(datanow);
    for it=1:1000
        topick1=randi(size(datanow,1),samplesize,1);
        d1=datanow(topick1,:);
        topick2=randi(size(datanow,1),samplesize,1);
        d2=datanow(topick2,:);
        
        p=ranksum(mean(d1,2),mean(d2,2));
        within(t,it)=p;
    end
end


% across species
across=zeros(3,1000);
cnt=0;
for t=1:2
    for tt=t+1:3
        cnt=cnt+1;
        datanow1=find(tID==t);
        datanow1=speedpref(datanow1);
        datanow2=find(tID==tt);
        datanow2=speedpref(datanow2);
        for it=1:1000
            topick1=randi(size(datanow1,1),samplesize,1);
            d1=datanow1(topick1,:);
            topick2=randi(size(datanow2,1),samplesize,1);
            d2=datanow2(topick2,:);
            
            p=ranksum(mean(d1,2),mean(d2,2));
            across(cnt,it)=p;
        end
    end
end

clc
median(across,2)
median(within,2)
%% heatmap per cluster and species


clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
file='clust_Kmeans_PC_F1new_corr';
load(fullfile(path1,file))
sel=0; clust=14; %actual clust
gau=1;
st=0;
fact=2;

if sel==1
    data=ID(:,clust-1);
end

superpath='F:\all_species_temp';
Dtype='wph';
normDG=[];
for ty=1:3
    load(fullfile(superpath,'newAttempt6','DG',[Dtype(ty),'_normPeaks']))
    normDG=[normDG;normPeaks];
end

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0.2 0.2 0.2; 0 1 1; 1 1 0;0.5 0 0];

clu=unique(data);
for clus=1:length(clu)
    currID=find(data==clu(clus));
    
    for t=1:3
        figure(t)
        tmp=find(tID(currID)==t);
        currnow=currID(tmp);
        
        if ~isempty(currnow)
            subplot(4,5,clus)
            curr=normDG(currnow,:,2);
            curr=mean(curr,1);
            curr=flipud(reshape(curr,4,6));
            [y,x] = meshgrid(1:6,1:4);
            imagesc(y(:)',x(:)',curr);
            colormap('gray')
            
            
            if gau==1
                
                figure(10+t)
                subplot(4,5,clus)
                curr=normDG(currnow,:,2);
                curr=mean(curr,1);
                curr=flipud(reshape(curr,4,6));
                [y,x] = meshgrid(1:6,1:4);
                imagesc(y(:)',x(:)',curr);
                hold on
                colormap('gray')
                title([int2str(clus),' (mean; ',int2str(length(currID)),'cells)'])
                
                
                
                
                curr=normDG(currnow,:,2);
                curr=mean(curr,1);
                dat=flipud(reshape(curr,4,6));
                dat=flipud(dat);
                
                minData = mean(mean(dat))+st*std(reshape(dat,24,1));
                p = (dat - minData);
                tmp=find(p<0);
                p(tmp)=0;
                sumData = sum(sum(p));
                p=p/sumData;
                
                
                [y,x] = meshgrid(6:-1:1,1:4);
                xx=reshape(x,24,1);
                yy=reshape(y,24,1);
                pp=reshape(p,24,1);
                mx=dot(xx,pp); %centerpoint
                my=dot(yy,pp); %centerpoint
                
                a=dot((xx-mx).^2,pp);
                b=dot((xx-mx).*(yy-my),pp);
                c=dot((yy-my).^2,pp);
                co=[a,b;b,c];
                
                if fact>=1
                    D=fact*sqrt(eigs(co));
                else
                    D=sqrt(eigs(co));%Eigenvalues = two main axes
                end
                [V,dd]=eigs(co);
                ma=max(max(dd));
                [ma1 ma2]=find(dd==ma);
                %         ma1
                V=V(ma1,:); v=V(1)/V(2);
                
                Ang=atan2(V(2),V(1)); %used
                %                Ang=atan2(V(1),V(2));
                AngD=rad2deg(Ang); %angle in degrees
                
                ymin=get(gca,'ylim');
                ymin=ymin(1);
                
                centerX=mx;%x-axis center;1=left
                centerY=my;%y-axis center; 1=bottom
                
                
                
                centerX2=5-mx;%x-axis center;1=left
                centerY2=7-my;%y-axis center; 1=bottom
                Ang2=-Ang;
                [X,Y]=meshgrid(0.5:0.01:4.5,0.5:0.01:6.5);
                
                theta = pi/2-Ang;
                %
                aa = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
                bb = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
                cc = sin(theta)^2/2/(D(1)/1.5)^2 + cos(theta)^2/2/(D(2)/1.5)^2;
                allClust = exp( - (aa*(Y-centerY2).^2 + 2*bb*(Y-centerY2).*(X-centerX2) + cc*(X-centerX2).^2)) ;
                contour(Y,X,allClust(:,:),1,'color','r')
                hold on
                plot(7-centerY,5-centerX,'or','markerfacecolor','r','markersize',2)
                
            end
        end
        
    end
end

%% PER EXP distribution in spatio-temporal clusters



clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
file='clust_Kmeans_PC_F1new_corr';
load(fullfile(path1,file))
sel=0; clust=14; %actual clust
gau=1;
st=0;
fact=2;

if sel==1
    data=ID(:,clust-1);
end

superpath='F:\all_species_temp';
Dtype='wph';

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0.2 0.2 0.2; 0 1 1; 1 1 0;0.5 0 0];

clu=unique(data);

for t=1:3
    load(fullfile(superpath,[Dtype(t),'_info']))
    load(fullfile(superpath,'newAttempt6','DG',[Dtype(t),'_normPeaks']))
    dates=info(togo,1);
    
    currnow=find(tID==t);
    currdata=data(currnow);
    
    figure(t)
    
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    start=[start;length(currnow)+1];
    udates=udates(ia);
    
    for u=1:length(udates)
        subplot(4,4,u)
        
        nums=[];
        for clus=1:length(clu)
            ID=find(currdata==clus);
            tmp=find(ID>=start(u));
            tmp2=find(ID(tmp)<start(u+1));
            nums=[nums;length(tmp2)];
        end
        
        bar(nums)
        title(udates{u})
    end
end
%% electrode maps

clear
close all

path1='F:\all_species_temp';
types='wph';

FLASH=zeros(3,1);
LF=zeros(3,1);
DS=zeros(3,1);
DG=zeros(3,1);
totCells=zeros(3,1);

clc
for t=1:length(types)
    
    if t==1
        currpath='S:\data\Katja\WT';
    elseif t==2
        currpath='S:\data\Katja\Pig';
        
    else
        currpath='S:\data\Katja\Human';
    end
    
    
    load(fullfile(path1,[types(t),'_info']))
    goodones=cell2mat(info(:,3));
    goodones=find(goodones>0);
    
    if t==2
        badIDs=202:250;
        bads=[];
        for b=1:length(badIDs)
            tmp=find(goodones==badIDs(b));
            bads=[bads;tmp];
        end
        goodones=setdiff(goodones,bads);
    end
    
    new=info(goodones,:);
    dates=new(:,1);
    udates=unique(dates);
    
    
    
    %flash
    flash=cell2mat(new(:,26));
    ID=find(flash==1);
    goodIDs=ID;
    
    
    %LF
    ID=cell2mat(new(:,9));
    ID=find(ID~=0);
    goodIDs=[goodIDs;ID];
    
    
    
    %DG
    load(fullfile(path1,'newAttempt3','DG',[types(t),'_normPeaks']))
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    % chirp
    load(fullfile(path1,'newAttempt2','chirp',['chirpRatios_',int2str(8000)]))
    togo=togos{t};
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    
    %6vel
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    ID=[];
    for to=1:length(goodB)
        tmp=find(goodones==goodB(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    ID=[];
    for to=1:length(goodW)
        tmp=find(goodones==goodW(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    
    goodIDs=unique(goodIDs);
    
    
    
    % now go through all dates
    dates=info(goodones,1);
    units=info(goodones,2);
    udates=unique(dates);
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    start=[start;length(dates)+1];
    udates=udates(ia);
    
    CHmap=[0 24 26 29 32 35 37 0; 21 22 25 30 31 36 39 40; 19 20 23 28 33 38 41 42; 16 17 18 27 34 43 44 45; 15 14 13 4 57 48 47 46;12 11 8 3 58 53 50 49;10 9 6 1 60 55 52 51;0 7 5 2 59 56 54 0];
    CHmap=flipud(CHmap);
    for u=1:length(udates)
        figure
        un=units(start(u):start(u+1)-1);
        upath=fullfile(currpath,udates{u},'units');
        ulist=dir(fullfile(upath,'*mat'));
        CH=[];
        for uu=1:length(un)
            m=0; found=0;
            while found==0
                m=m+1;
                if ~isempty(regexp(ulist(m).name,un{uu}))
                    found=1;
                end
            end
            name=ulist(m).name;
            tmp=regexp(name,'CH');
            tmp2=regexp(name(tmp:end),'_');
            num=str2num(name(tmp+2:tmp2(1)+tmp-2));
            CH=[CH;num];
        end
        
        
        goodies=[];
        for s=start(u):start(u+1)-1
            if ~isempty(find(goodIDs==s))
                goodies=[goodies;1];
            else
                goodies=[goodies;0]; %responding cells
            end
        end
        
        for xx=1:8
            for yy=1:8
                plot(xx,yy,'o','color',[0.6 0.6 0.6],'markerfacecolor',[0.6 0.6 0.6],'markersize',10)
                hold on
            end
        end
        plot(1,1,'ow','markerfacecolor','w','markersize',10)
        plot(8,1,'ow','markerfacecolor','w','markersize',10)
        plot(1,8,'ow','markerfacecolor','w','markersize',10)
        plot(8,8,'ow','markerfacecolor','w','markersize',10)
        
        bad=find(goodies==0);
        for b=1:length(bad)
            [xx,yy]=find(CHmap==CH(bad(b)));
            plot(xx,yy,'o','color','r','markerfacecolor','r','markersize',45)
        end
        nice=find(goodies==1);
        for n=1:length(nice)
            [yy,xx]=find(CHmap==CH(nice(n)));
            plot(xx,yy,'o','color','g','markerfacecolor','g','markersize',28)
        end
        axis([0 9 0 9])
        title([types(t),' / ', udates{u}])
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile('F:\all_species_temp\CHmaps',[types(t),'_',udates{u},'.emf']))
        
    end
    
end
%% bootstrap speed pref.bar
% result: m-p (0.038), m-h (0.003), p-h (0.381)
% result: m (0.507), p (0.502), h (0.511)
clear
close all
superpath='f:\all_species_temp';

types='wph';
cols='bgr';
titles={'mouse','pig','human'};

speedpref=[]; realIDs=[]; tID=[]; dates=[];
for t=1:length(types)
    
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    date=info(:,1);
    allval=zeros(size(info,1),2);
    black=cell2mat(info(goodB,24));
    white=cell2mat(info(goodW,25));
    allval(goodB,1)=black;
    allval(goodW,2)=white;
    allvalorig=max(allval,[],2);
    tmp=find(allvalorig>0);
    date=date(tmp);
    dates=[dates;date];
    allvalorig=allvalorig(tmp);
    realIDs=[realIDs;tmp];
    tID=[tID;repmat(t,length(tmp),1)];
    speedpref=[speedpref;allvalorig];
end

save(fullfile(superpath,'newAttempt','speed_bar','speedBar'),'speedpref','tID','realIDs')

dlength=[];
for t=1:3
    tpm=find(tID==t);
    
    dlength=[dlength;length(tpm)];
end

samplesize=round(length(tID)/length(unique(dates)));
% samplesize=round(37/2);

% within species
within=zeros(3,1000);
for t=1:3
    datanow=find(tID==t);
    datanow=speedpref(datanow);
    for it=1:1000
        topick1=randi(size(datanow,1),samplesize,1);
        d1=datanow(topick1,:);
        topick2=randi(size(datanow,1),samplesize,1);
        d2=datanow(topick2,:);
        
        p=ranksum(mean(d1,2),mean(d2,2));
        within(t,it)=p;
    end
end


% across species
across=zeros(3,1000);
cnt=0;
for t=1:2
    for tt=t+1:3
        cnt=cnt+1;
        datanow1=find(tID==t);
        datanow1=speedpref(datanow1);
        datanow2=find(tID==tt);
        datanow2=speedpref(datanow2);
        for it=1:1000
            topick1=randi(size(datanow1,1),samplesize,1);
            d1=datanow1(topick1,:);
            topick2=randi(size(datanow2,1),samplesize,1);
            d2=datanow2(topick2,:);
            
            p=ranksum(mean(d1,2),mean(d2,2));
            across(cnt,it)=p;
        end
    end
end

clc
median(across,2)
median(within,2)
%% plot example cells: flash
clear
close all
load(fullfile('F:\all_species_temp','h_flashSpikes'))
load(fullfile('F:\all_species_temp','h_flash40sigma'))

ids=[23 56 103 119 481 850];
stim=[2000 4000 6000 8000 1000];



for i=1:length(ids)
    data=flash_raster(ids(i),:);
    figure
    for d=1:length(data)
        spikes=data{d};
        for s=1:length(spikes)
            plot([spikes(s) spikes(s)],[d-0.4 d+0.4],'-')
            hold on
        end
    end
    for st=1:length(stim)
        plot([stim(st) stim(st)],[0 d+1],'-k')
    end
    set(gcf,'position',[1601 1 1600 824])
    saveas(gcf,fullfile('F:\all_species_temp\final_pics',['_exFlash_',int2str(ids(i)),'_raster.emf']))
    close
    
    figure
    plot(flash_ratesSmoo{ids(i)})
    hold on
    for st=1:length(stim)
        plot([stim(st) stim(st)],[0 max(flash_ratesSmoo{ids(i)})+1],'-k')
    end
    set(gcf,'position',[1601 1 1600 824])
    saveas(gcf,fullfile('F:\all_species_temp\final_pics',['_exFlash_',int2str(ids(i)),'_rate40sigma.emf']))
    close
    
end
%% plot example cells: Y-like 694
% plot 500 um 2 Hz and 4000 um 2 Hz
close all
clear

path1='e:\all_species_temp\newAttempt3';

% load(fullfile(path1,'DG_PCanalysis','clust_Kmeans_PC_F1new_corr'))
load(fullfile(path1,'DG','h_normPeaks'))
% load(fullfile('F:\all_species_temp','polarity'))
% load(fullfile('F:\all_species_temp','h_vel6_new_Black'))
% load(fullfile('F:\all_species_temp','h_chirp'))
% load(fullfile('F:\all_species_temp','h_flash60sigma'))
load(fullfile('e:\all_species_temp\newAttempt\DG','h_DGrates'))
load(fullfile('e:\all_species_temp\newAttempt\DG','h_DGFFT_binary'))
fft=FFTabs;
load(fullfile('e:\all_species_temp\newAttempt\DG','h_DGFFT_freq'))

id=694;
id2=find(togo==id);

figure
plot(DGRates{id2,4,10})
hold on
fs2=500;
T=12;
phase=pi/2;
frq1=2;
t=1/fs2:1/fs2:T;
resp1=sin(2*pi*frq1*t+phase);
plot(t*1000,resp1*max(DGRates{id2,4,10})/2+max(DGRates{id2,4,10})/2,'-k');

figure
plot(DGRates{id2,4,22})
hold on
fs2=4000;
T=12;
phase=pi/2;
frq1=2;
t=1/fs2:1/fs2:T;
resp1=sin(2*pi*frq1*t+phase);
plot(t*1000,resp1*max(DGRates{id2,4,22})/2+max(DGRates{id2,4,22})/2,'-k');
%% stimuli fig

path1='S:\data\Katja\Human\20131210a\HEKA';
path2='F:\all_species_temp\final_pics';
hlist=dir(fullfile(path1,'*phys'));

chirp=82;
prot=read_header_field_heka(path1,hlist(chirp).name,'Stimulus Protocol');

int=prot(3:end-1,4);

time=prot(2:end-2,1);
timediff=diff(time);
firstpart=find(timediff>50);

firstpart=firstpart(1);

su=sum(timediff);
timediff=round(timediff);
roundsum=sum(timediff);
changing=ceil(roundsum-su);
stepping=floor(length(timediff)/changing);
for ti=1:stepping:length(timediff)
    timediff(ti)=16;
end

int2=[]; idds=[];
for ii=1:firstpart
    tmp=repmat(int(ii),round(timediff(ii)),1);
    int2=[int2;tmp];
    idds=[idds;round(timediff(ii))];
end

figure
plot(int2)
saveas(gcf,fullfile(path2,'_stim_chirp.emf'))


figure
subplot(2,1,1)
fs2=100;
T=12;
phase=pi/2;
frq1=1;
t=1/fs2:1/fs2:T;
resp1=sin(2*pi*frq1*t+phase);
imagesc(repmat(resp1,100,1))
colormap('gray')
subplot(2,1,2)
fs2=100;
T=12;
phase=pi/2;
frq1=0.25;
t=1/fs2:1/fs2:T;
resp1=sin(2*pi*frq1*t+phase);
imagesc(repmat(resp1,100,1))
colormap('gray')
saveas(gcf,fullfile(path2,'_stim_DG.emf'))

wnf=22;
prot=read_header_field_heka(path1,hlist(wnf).name,'Stimulus Protocol');

figure
plot(prot(2:end-1,1),prot(2:end-1,4))
saveas(gcf,fullfile(path2,'_stim_WNF.emf'))
%% example speed

%% plot for each cell whole flash
clear
close all

path1='e:\all_species_temp';
epath='S:\data\Katja\Human';

load(fullfile(path1,['h_info']))
goodones=cell2mat(info(:,3));
goodones=find(goodones>0);


new=info(goodones,:);
dates=new(:,1);
units=new(:,2);
udates=unique(dates);



%flash
flash=cell2mat(new(:,26));
ID=find(flash==1);
goodIDs=ID;

for g=1:length(goodIDs)
    upath=fullfile(epath,dates{goodIDs(g)},'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,units{goodIDs(g)}))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlash'))
            flashfiles=[flashfiles;h];
        end
    end
    
    figure
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(spikes,40,15000);
        allrate=[allrate;rate];
        subplot(1,2,1)
        plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
        hold on
    end
    axis([0 15000 0 f+1])
    title(ulist(m).name,'interpreter','none')
    subplot(1,2,2)
    meanrate=mean(allrate,1);
    plot(meanrate)
    
    set(gcf,'position',[1601 1 1600 824])
    realID=goodones(goodIDs(g));
    saveas(gcf,fullfile('E:\all_species_temp\overviewFlashResponses',[int2str(realID),'_',ulist(m).name(1:end-4),'.bmp']))
    close
end

%% plot for each cell whole drifting-grating
clear
close all

path1='e:\all_species_temp';
epath='S:\data\Katja\Human';

load(fullfile(path1,['h_info']))
goodones=cell2mat(info(:,3));
goodones=find(goodones>0);


new=info(goodones,:);
dates=new(:,1);
units=new(:,2);
udates=unique(dates);



%flash
flash=cell2mat(new(:,26));
ID=find(flash==1);
goodIDs=ID;

for g=1:length(goodIDs)
    upath=fullfile(epath,dates{goodIDs(g)},'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,units{goodIDs(g)}))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'DriftGrat_Contr1_2000um%dur_12000_speed_4'))
            flashfiles=[flashfiles;h];
        end
    end
    
    figure
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(spikes,40,13000);
        subplot(2,2,[2 4])
        if ~isempty(find(rate>0))
            plot(rate/max(rate)+f)
        end
        hold on
        allrate=[allrate;rate];
        subplot(2,2,1)
        plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
        hold on
    end
    plot([0 100],[f+1 f+1],'-r')
    axis([0 13000 0 f*3])
    title(ulist(m).name,'interpreter','none')
    subplot(2,2,3)
    meanrate=mean(allrate,1);
    plot(meanrate)
    
    set(gcf,'position',[1601 1 1600 824])
    realID=goodones(goodIDs(g));
    saveas(gcf,fullfile('E:\all_species_temp\overviewDGResponses',[int2str(realID),'_',ulist(m).name(1:end-4),'.bmp']))
    close
end
%% plot for each cell whole drifting-grating2
clear
close all

path1='e:\all_species_temp';
epath='S:\data\Katja\Human';

load(fullfile(path1,['h_info']))
goodones=cell2mat(info(:,3));
goodones=find(goodones>0);


new=info(goodones,:);
dates=new(:,1);
units=new(:,2);
udates=unique(dates);



%flash
flash=cell2mat(new(:,26));
ID=find(flash==1);
goodIDs=ID;

for g=1:length(goodIDs)
    upath=fullfile(epath,dates{goodIDs(g)},'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,units{goodIDs(g)}))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'DriftGrat_Contr1_4000um%dur_12000_speed_4'))
            flashfiles=[flashfiles;h];
        end
    end
    
    figure
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(spikes,40,13000);
        subplot(2,2,[2 4])
        if ~isempty(find(rate>0))
            plot(rate/max(rate)+f)
        end
        hold on
        allrate=[allrate;rate];
        subplot(2,2,1)
        plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
        hold on
    end
    plot([0 100],[f+1 f+1],'-r')
    axis([0 13000 0 f*3])
    title(ulist(m).name,'interpreter','none')
    subplot(2,2,3)
    meanrate=mean(allrate,1);
    plot(meanrate)
    
    set(gcf,'position',[1601 1 1600 824])
    realID=goodones(goodIDs(g));
    saveas(gcf,fullfile('E:\all_species_temp\overviewDGResponses2',[int2str(realID),'_',ulist(m).name(1:end-4),'.bmp']))
    close
end

%% plot for each cell whole chirp
clear
close all

path1='e:\all_species_temp';
epath='S:\data\Katja\Human';

load(fullfile(path1,['h_info']))
goodones=cell2mat(info(:,3));
goodones=find(goodones>0);


new=info(goodones,:);
dates=new(:,1);
units=new(:,2);
udates=unique(dates);



%flash
load(fullfile(path1,'newAttempt2','chirp',['chirpRatios_',int2str(8000)]))
togo=togos{3};
ID=[];
for to=1:length(togo)
    tmp=find(goodones==togo(to));
    ID=[ID;tmp];
end
goodIDs=ID;

for g=1:length(goodIDs)
    upath=fullfile(epath,dates{goodIDs(g)},'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,units{goodIDs(g)}))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'chirp'))
            flashfiles=[flashfiles;h];
        end
    end
    
    figure
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(spikes,40,13000);
        subplot(2,2,4)
        if ~isempty(find(rate>0))
            plot(rate/max(rate)+f)
        end
        hold on
        allrate=[allrate;rate];
        subplot(2,2,[1 2])
        plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
        hold on
    end
    plot([0 100],[f+1 f+1],'-r')
    axis([0 13000 0 f*3])
    title(ulist(m).name,'interpreter','none')
    subplot(2,2,3)
    meanrate=mean(allrate,1);
    plot(meanrate)
    
    set(gcf,'position',[1601 1 1600 824])
    realID=goodones(goodIDs(g));
    saveas(gcf,fullfile('E:\all_species_temp\overviewChirpResponses',[int2str(realID),'_',ulist(m).name(1:end-4),'.bmp']))
    close
end
%% NEW: plot examples spatio-temporal

clear
close all
path1='e:\all_species_temp\newAttempt3\DG_PCanalysis';
file='clust_Kmeans_PC_F1new_corr';
load(fullfile(path1,file))
sel=0; clust=14; %actual clust
gau=1;
st=0;
fact=2;


myexamples=[388 411 473];
mystimuli=[3 11 23];
myfreq=[4 4 4];



if sel==1
    data=ID(:,clust-1);
end

superpath='e:\all_species_temp';
Dtype='wph';
cols='bgr';
normDG=[];
for ty=1:3
    load(fullfile(superpath,'newAttempt3','DG',[Dtype(ty),'_normPeaks']))
    normDG=[normDG;normPeaks];
end

path2='S:\Paper factory\Katja\5_HUMPIG\newFigs';
path3='E:\all_species_temp';
load(fullfile(path3,'h_DGrates'));

load(fullfile('e:\all_species_temp\newAttempt\DG','h_DGFFT_binary'))
fft=FFTabs;
load(fullfile('e:\all_species_temp\newAttempt\DG','h_DGFFT_freq'))

ttID=find(tID==3);
frs=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];


for ex=1:length(myexamples)
    tmp=find(ttID==myexamples(ex));
    tmp3=togo(tmp);
    tmp2=myexamples(ex);
    
    
    
    figure
    
    curr=normDG(tmp2,:,2);
    curr=flipud(reshape(curr,4,6));
    [x,y] = meshgrid(1:6,1:4);
    imagesc(x(:)',y(:)',curr);
    hold on
    colormap('gray')
    set(gcf,'position',[1601 1 1600 824])
    saveas(gcf,fullfile(path2,['ex_spatiotemp_',int2str(myexamples(ex)),'_heatmap.emf']))
    close
    
    
    figure
    for ss=1:length(mystimuli)
        subplot(2,2,ss)
        rate=DGrates{tmp3,mystimuli(ss)};
        plot(rate)
    end
    set(gcf,'position',[1601 1 1600 824])
    saveas(gcf,fullfile(path2,['ex_spatiotemp_',int2str(myexamples(ex)),'_rates.emf']))
    close
    
    figure
    for ss=1:length(mystimuli)
        subplot(2,2,ss)
        freq=FFTabs{tmp,6,mystimuli(ss)};
        FFTresp=fft{tmp,6,mystimuli(ss)};
        
        
        
        tp1=find(freq>=0.35);
        tp2=find(freq<=30);
        data=FFTresp(tp1(1):tp2(end));
        freq2=freq(tp1(1):tp2(end));
        
        tempData=data;
        excluded=[];
        for f=[0.5 1 2 3];
            tp=find(freq2==frs(mystimuli(ss))*f);
            excluded=[excluded tp-2:tp+2];
        end
        included=1:length(tempData);
        start=find(freq2<=frs(mystimuli(ss))*0.5);
        start=start(end)-5;
        if start<0
            start=1;
        end
        endpoint=find(freq2==frs(mystimuli(ss))*3);
        endpoint=endpoint+5;
        included=included(start:endpoint);
        % included=setdiff(included,excluded-start);
        included=setdiff(included,excluded);
        tempData=tempData(included);
        
        bckr=mean(tempData);
        bckrSTD=std(tempData);
        thr=bckr+3*bckrSTD;
        
        
        
        
        
        
        plot(freq,FFTresp)
        hold on
        tst=find(freq==myfreq(ss));
        val=FFTresp(tst);
        plot(myfreq(ss),val,'or')
        plot([0 10],[thr thr],'-g')
        axis([0 10 0 0.01])
    end
    set(gcf,'position',[1601 1 1600 824])
    saveas(gcf,fullfile(path2,['ex_spatiotemp_',int2str(myexamples(ex)),'_ffts.emf']))
    close
    
    
end
%% NEW: plot examples for flash
clear
close all

path1='e:\all_species_temp';
epath='S:\data\Katja\Human';

load(fullfile(path1,['h_info']))
goodones=cell2mat(info(:,3));
goodones=find(goodones>0);


new=info(goodones,:);
dates=new(:,1);
units=new(:,2);
udates=unique(dates);



%flash
flash=cell2mat(new(:,26));
ID=find(flash==1);
goodIDs=ID;


toplot=[149 390 481 855 103 867 357];

for g=1:length(toplot)
    otherID=find(goodones==toplot(g));
    currID=find(goodones(goodIDs)==toplot(g));
    upath=fullfile(epath,dates{otherID},'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    
    
    
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,units{otherID}))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlash'))
            flashfiles=[flashfiles;h];
        end
    end
    
    figure
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(spikes,40,15000);
        rate=rate(200:10500);
        allrate=[allrate;rate];
        spikes=spikes(find(spikes>=200));
        spikes=spikes-200+1;
        spikes=spikes(find(spikes<=10300));
        subplot(2,1,1)
        for s=1:length(spikes)
            plot([spikes(s) spikes(s)],[f-0.4 f+0.4],'-b')
            hold on
            
        end
        %         plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
        %         hold on
    end
    plot([1800 1800],[0 f+1],'-k')
    plot([3800 3800],[0 f+1],'-k')
    plot([5800 5800],[0 f+1],'-k')
    plot([7800 7800],[0 f+1],'-k')
    plot([9800 9800],[0 f+1],'-k')
    axis([0 10500 0 f+1])
    title(ulist(m).name,'interpreter','none')
    subplot(2,1,2)
    meanrate=mean(allrate,1);
    plot(meanrate)
    axis([0 10500 0 inf])
    
    set(gcf,'position',[1601 1 1600 824])
    saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_flash_',int2str(toplot(g)),'_',ulist(m).name(1:end-4),'.emf']))
    close
end
%% NEW: plot examples for chirp
clear
close all

path1='e:\all_species_temp';
epath='S:\data\Katja\Human';

load(fullfile(path1,['h_info']))
goodones=cell2mat(info(:,3));
goodones=find(goodones>0);


new=info(goodones,:);
dates=new(:,1);
units=new(:,2);
udates=unique(dates);

toplot=[342 609 721 757 763 867];


load(fullfile(path1,'newAttempt2','chirp',['chirpRatios_',int2str(8000)]))
togo=togos{3};
ID=[];
for to=1:length(togo)
    tmp=find(goodones==togo(to));
    ID=[ID;tmp];
end
goodIDs=ID;

for g=1:length(toplot)
    otherID=find(goodones==toplot(g));
    currID=find(goodones(goodIDs)==toplot(g));
    upath=fullfile(epath,dates{otherID},'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,units{otherID}))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'chirp'))
            flashfiles=[flashfiles;h];
        end
    end
    prot=read_header_field_heka(fullfile(epath,dates{otherID},'HEKA'),unit{1,2}{flashfiles(1),1},'Stimulus Protocol');
    start=prot(2,1);
    start=start-200;
    
    
    if g==1
        
        
        int=prot(2:end-1,4);
        int=[128; int];
        
        time=prot(2:end-2,1);
        time=time-start;
        tmp=find(time<=10000);
        time=time(tmp);
        timediff=diff(time);
        timediff=timediff(1:482);
        
        
        su=sum(timediff);
        timediff=round(timediff);
        roundsum=sum(timediff);
        changing=ceil(roundsum-su);
        stepping=floor(length(timediff)/changing);
        for ti=1:stepping:length(timediff)+2
            timediff(ti)=16;
        end
        tmp=find(timediff>40);
        timediff(tmp)=16;
        
        int2=[]; idds=[];
        for ii=1:length(timediff)-1
            tmp=repmat(int(ii),round(timediff(ii)),1);
            int2=[int2;tmp];
            idds=[idds;round(timediff(ii))];
        end
        int2=[repmat(128,round(time(1)),1); int2];
        int2=[int2;repmat(128,10000-length(int2),1)];
        
        figure
        plot(int2)
        set(gcf,'position',[1601 1 1600 824])
        saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_chirp_stimulus.emf']))
        close
    end
    
    figure
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(spikes,10,13000);
        rate=rate(start:10000+start);
        allrate=[allrate;rate];
        subplot(2,1,1)
        spikes=spikes(find(spikes>=start));
        spikes=spikes-start+1;
        spikes=spikes(find(spikes<=10000));
        for s=1:length(spikes)
            plot([spikes(s) spikes(s)],[f-0.4 f+0.4],'-b')
            hold on
            
        end
        %         plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
        %         hold on
    end
    %     plot([0 100],[f+1 f+1],'-r')
    axis([0 10000 0 f+1])
    title(ulist(m).name,'interpreter','none')
    subplot(2,1,2)
    meanrate=mean(allrate,1);
    plot(meanrate)
    axis([0 10000 0 inf])
    
    set(gcf,'position',[1601 1 1600 824])
    realID=goodones(goodIDs(g));
    saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_chirp_',int2str(toplot(g)),'_',ulist(m).name(1:end-4),'.emf']))
    close
end
%% NEW: plot Y-cells: work on it!

close all
clear

path1='e:\all_species_temp\newAttempt3';

load(fullfile(path1,'DG_PCanalysis','clust_Kmeans_PC_F1new_corr'))
load(fullfile(path1,'DG','h_normPeaks'))
load(fullfile('e:\all_species_temp','polarity'))
load(fullfile('e:\all_species_temp','h_vel6_new_Black'))
load(fullfile('e:\all_species_temp','h_chirp'))
load(fullfile('e:\all_species_temp\newAttempt\DG','h_DGFFT_binary'))
fft=FFTabs;
load(fullfile('e:\all_species_temp\newAttempt\DG','h_DGFFT_freq'))
load(fullfile('e:\all_species_temp','h_DGrates'));

dataH=find(tID==3);
dataH=data(dataH);

tmp=find(type==3);
IDh=flashID(tmp);
pol=onoff(tmp);


% all ratios

rat=normPeaks(:,[10 11 12 14 15 16 18 19 20],3)./normPeaks(:,[10 11 12 14 15 16 18 19 20],2);

int=[];
for r=1:9
    tmp=find(rat(:,r)>=1);
    int=[int;tmp];
end
int=unique(int);
clust=dataH(int);
[clust id]=sort(clust);
int=int(id);

list=9:21;
list2=1:24;


equi2=[1 2 3 4 1 2 3 4 1 2 3 4];
rat=normPeaks(int,:,3)./normPeaks(int,:,2);


% F1=[]; F2=[];
%
% for f=1:length(int)
%     for r=1:size(rat,2)
%         tmp=find(rat(f,r)>=1);
%         if ~isempty(tmp)
%
%             F1=[F1; normPeaks(int(f),[1:4:21]+equi2(r)-1,2)];
%             F2=[F2; normPeaks(int(f),[1:4:21]+equi2(r)-1,3)];
%         end
%     end
% end

totake=[29 28 25 23 19 13 11 5 4 3 1];
int=int(totake);
rat=rat(totake,:);
F500=[]; F4000=[];
allFs=zeros(length(int),6);
for f=1:length(int)
    
    
    
    curr1=[]; curr2=[];
    for s=1:6 %go through spatials
        %         if ~isempty(find(rat(f,(s-1)*4+1:(s-1)*4+4)>1))
        %             now=find(rat(f,(s-1)*4+1:(s-1)*4+4)>1);
        %             %             tmp=find(rat(f,:)>=1);
        %             %             tmp2=list(tmp);
        %             %             now2=intersect(now,tmp2);
        
        now=1:4;
        
        currlist=(s-1)*4+1:(s-1)*4+4;
        [t1 t2]=max(normPeaks(int(f),currlist(now),3));
        t2=currlist(now(t2));
        if t1==0
            [t1 t2]=max(normPeaks(int(f),now,2));
            t2=currlist(now(t2));
        end
        
        idx=t2;
        allFs(f,s)=idx;
        if s==3
            F500=[F500;idx];
        elseif s==6
            
            F4000=[F4000;idx];
        end
    end
end

F500=[19 14 11 14 14 19 13 15 16 16 19];
tosave=[5 8 11];
temp=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
for f=1:length(int)
    f500=fft{int(f),6,F500(f)};
    f4000=fft{int(f),6,F4000(f)};
    ma=max(max(f500(10:301)), max(f4000(10:301)));
    rate500=DGrates{togo(int(f)),F500(f)}(3000:14000);
    rate4000=DGrates{togo(int(f)),F4000(f)}(3000:14000);
    ma2=max(max(rate500),max(rate4000));
    figure
    subplot(2,2,1)
    plot(rate500)
    axis([0 11000 0 ma2+ma2/10])
    subplot(2,2,3)
    
    plot(FFTabs{1,6,1},f500)
    hold on
    tmp=temp(F500(f));
    tst=find(FFTabs{1,6,1}==tmp);
    tst2=find(FFTabs{1,6,1}==tmp*2);
    plot(tmp,f500(tst),'or')
    plot(tmp*2,f500(tst2),'or')
    
    axis([0 20 0 ma+ma/10])
    title([int2str(togo(int(f))),' / ',int2str(F500(f))])
    
    
    subplot(2,2,2)
    plot(rate4000)
    axis([0 11000 0 ma2+ma2/10])
    subplot(2,2,4)
    
    plot(FFTabs{1,6,1},f4000)
    hold on
    tmp=temp(F4000(f));
    tst=find(FFTabs{1,6,1}==tmp);
    tst2=find(FFTabs{1,6,1}==tmp*2);
    plot(tmp,f4000(tst),'or')
    plot(tmp*2,f4000(tst2),'or')
    axis([0 20 0 ma+ma/10])
    title([int2str(togo(int(f))),' / ',int2str(F4000(f))])
    
    if ~isempty(find(tosave==f))
        set(gcf,'position',[1601 1 1600 824])
        saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_Ycell_',int2str(togo(int(f))),'.emf']))
        
    end
end



for i=1:length(int(tosave))
    f1=normPeaks(int(tosave(i)),:,2);
    f2=normPeaks(int(tosave(i)),:,3);
    f2f1=f2./f1;
    figure
    curr=flipud(reshape(f2f1,4,6));
    [x,y] = meshgrid(1:6,1:4);
    imagesc(x(:)',y(:)',curr,[0 1]);
    hold on
    colormap('gray')
    saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_Ycell_',int2str(togo(int(tosave(i)))),'_ratiosF2F1.emf']))
    
    
    figure
    curr=flipud(reshape(f1,4,6));
    [x,y] = meshgrid(1:6,1:4);
    imagesc(x(:)',y(:)',curr);
    hold on
    colormap('gray')
    saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_Ycell_',int2str(togo(int(tosave(i)))),'_F1heatmap.emf']))
    
    
    
    
end

% four=fft{int(tosave(1)),6,1};
% figure
% plot(FFTabs{1,6,1},four)


% for f=1:length(int)
%     figure
%     subplot(1,2,1)
%     f500=fft{int(f),6,F500(f)};
%     plot(FFTabs{1,6,1},f500)
%     axis([0 20 0 inf])
%     title([int2str(togo(int(f))),' / ',int2str(F500(f))])
%     subplot(1,2,2)
%     f4000=fft{int(f),6,F4000(f)};
%     plot(FFTabs{1,6,1},f4000)
%     axis([0 20 0 inf])
%     title([int2str(togo(int(f))),' / ',int2str(F4000(f))])
% end


% for f=1:length(int)
%     figure
%     for s=1:4
%         subplot(2,2,s)
%         f500=fft{int(f),6,8+s};
%         plot(FFTabs{1,6,1},f500)
%         axis([0 20 0 inf])
%         title(int2str(togo(int(f))))
%     end
% end
%
%
% for f=1:length(int)
%     figure
%     for s=1:4
%         subplot(2,2,s)
%         f4000=fft{int(f),6,20+s};
%         plot(FFTabs{1,6,1},f4000)
%
%         title(int2str(togo(int(f))))
%     end
% end

% for f=1:length(int)
%     figure
%     cnt=0;
%     for s=1:6
%         if allFs(f,s)>0
%             cnt=cnt+1;
%         subplot(2,2,cnt)
%         curr=fft{int(f),6,allFs(f,s)};
%         plot(FFTabs{1,6,1},curr)
%         title([int2str(togo(int(f))),' / ',int2str(allFs(f,s))])
%         end
%     end
% end
%% scatter plot speed pref bar vs DG

clear
close all
superpath='e:\all_species_temp';


speeds=[0.2; 1; 2; 4; 8; 16];


load(fullfile(superpath,['h','_info']))

goodB=[]; goodW=[];
for i=1:length(info)
    if ~isempty(info{i,24})
        if info{i,24}>0
            goodB=[goodB;i];
        end
    end
    if ~isempty(info{i,25})
        if info{i,25}>0
            goodW=[goodW;i];
        end
    end
end

allval=zeros(size(info,1),2);
black=cell2mat(info(goodB,24));
white=cell2mat(info(goodW,25));
allval(goodB,1)=black;
allval(goodW,2)=white;
allval=max(allval,[],2);
barID=find(allval>0);
allval=allval(barID);





spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
speed=spat.*freq/1000; speed=speed';
uspeed=[1 2 4 8 16];

ids=[];
for u=1:length(uspeed)
    tmp=find(speed==uspeed(u));
    ids=[ids;tmp];
end



load(fullfile(superpath,['h','_DGspeed']))



dgID=[]; doubleBar=[];
for b=1:length(barID)
    tmp=find(togo==barID(b));
    if ~isempty(tmp)
        dgID=[dgID;tmp];
        doubleBar=[doubleBar b];
    end
end



figure
scatter(allval(doubleBar),DGspeed(dgID))
axis([0 10 0 10])
fitobject = fit(allval(doubleBar),DGspeed(dgID),'poly1');
aFit=fitobject.p1;
bFit=fitobject.p2;

hold on
plot([0 10],[0 10],'-k')
% plot([0 10],[bFit aFit*10+bFit],'-k')

xlabel('bar speed [mm/s]')
ylabel('DG speed [mm/s]')
%% more examples for speed tuning


clear
close all
superpath='e:\all_species_temp';


speeds=[0.2; 1; 2; 4; 8; 16];


load(fullfile(superpath,['h','_info']))
load(fullfile(superpath,['h','_peaks6vel_black']))
blackAll=peak6vel;
load(fullfile(superpath,['h','_vel6_new_Black']))
% load(fullfile(superpath,['h','_DGrates']))


goodB=[]; goodW=[];
for i=1:length(info)
    if ~isempty(info{i,24})
        if info{i,24}>0
            goodB=[goodB;i];
        end
    end
    if ~isempty(info{i,25})
        if info{i,25}>0
            goodW=[goodW;i];
        end
    end
end

allval=zeros(size(info,1),2);
black=cell2mat(info(goodB,24));
white=cell2mat(info(goodW,25));
allval(goodB,1)=black;
allval(goodW,2)=white;
[allval idbar]=max(allval,[],2);
barID=find(allval>0);
allval=allval(barID);
idbar=idbar(barID);




load(fullfile(superpath,['h','_DGspeed']))

dgID=[]; doubleBar=[];
for b=1:length(barID)
    tmp=find(togo==barID(b));
    if ~isempty(tmp)
        dgID=[dgID;tmp]; %ids for DG
        doubleBar=[doubleBar b]; %ids for bar
    end
end





spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
speed=spat.*freq/1000; speed=speed';
uspeed=[1 2 4 8 16];

ids=[];
for u=1:length(uspeed)
    tmp=find(speed==uspeed(u));
    ids=[ids;tmp];
end

speed=speed(ids);

load(fullfile(superpath,'newAttempt2','DG',['h','_FFTpeaks']))
load(fullfile(superpath,'newAttempt3','DG',['h','_normPeaks']))
data=reshape(FFTpeaks(:,:,2),size(FFTpeaks,1),24);
data=data./repmat(max(data,[],2),1,24);
data=data(:,ids);


prot=read_header_field_heka('S:\data\Katja\Human\20120224\HEKA','20120224_C1#0014_OPA1_bar_6veloc_downup_black%speed_200_gap_3000%_FW0ND4.phys','Stimulus Protocol');
cons=prot(1:end,2);
changes=find(cons==-1);
startpoints=round(prot(changes(2:7),1));

newdata=[];
for s=1:length(uspeed)
    tmp=find(speed==uspeed(s));
    newdata=[newdata max(data(:,tmp),[],2)];
end

[~,sortID]=sort(allval(doubleBar));
% [~,minID]=min(allval(doubleBar));
% [~,maxID]=max(allval(doubleBar));
minID=31;
maxID=20;

minIDbar=barID(doubleBar(minID));
minIDdg=dgID(minID);
maxIDbar=barID(doubleBar(maxID));
maxIDdg=dgID(maxID);

figure
subplot(3,2,1)
plot([1 2 4 8 16],blackAll(minIDbar,2:end)/max(blackAll(minIDbar,2:end)),'-ok');
hold on
plot([1 2 4 8 16],newdata(minIDdg,:)/max(newdata(minIDdg,:)),'-ob');
axis([0 17 0 1])
title(['cell ',int2str(minIDbar), ' (',int2str(minIDdg),') / ',num2str(allval(doubleBar(minID))),', ',num2str(DGspeed(minIDdg))])
subplot(3,2,3)
plot(vel6_newB{minIDbar})
hold on
for s=1:length(startpoints)
    plot([startpoints(s) startpoints(s)],[0 max(vel6_newB{minIDbar})],'-r')
end
axis([15000 45000 0 Inf])
subplot(3,2,5)
curr=flipud(reshape(normPeaks(minIDdg,:,2),4,6));
[x,y] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr,[0 1]);
colormap('gray')
subplot(3,2,2)
plot([1 2 4 8 16],blackAll(maxIDbar,2:end)/max(blackAll(maxIDbar,2:end)),'-ok');
hold on
plot([1 2 4 8 16],newdata(maxIDdg,:)/max(newdata(maxIDdg,:)),'-ob');
axis([0 17 0 1])
title(['cell ',int2str(maxIDbar), ' (',int2str(maxIDdg),') / ',num2str(allval(doubleBar(maxID))),', ',num2str(DGspeed(maxIDdg))])
subplot(3,2,4)
plot(vel6_newB{maxIDbar})
hold on
for s=1:length(startpoints)
    plot([startpoints(s) startpoints(s)],[0 max(vel6_newB{maxIDbar})],'-r')
end
axis([15000 45000 0 Inf])
subplot(3,2,6)
curr=flipud(reshape(normPeaks(maxIDdg,:,2),4,6));
[x,y] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr,[0 1]);
colormap('gray')
set(gcf,'position',[1601 1 1600 824])
saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_speed_',int2str(minIDbar),'_',int2str(maxIDbar),'.emf']))


minID=4;
maxID=23;

minIDbar=barID(doubleBar(minID));
minIDdg=dgID(minID);
maxIDbar=barID(doubleBar(maxID));
maxIDdg=dgID(maxID);

figure
subplot(3,2,1)
plot([1 2 4 8 16],blackAll(minIDbar,2:end)/max(blackAll(minIDbar,2:end)),'-ok');
hold on
plot([1 2 4 8 16],newdata(minIDdg,:)/max(newdata(minIDdg,:)),'-ob');
axis([0 17 0 1])
title(['cell ',int2str(minIDbar), ' (',int2str(minIDdg),') / ',num2str(allval(doubleBar(minID))),', ',num2str(DGspeed(minIDdg))])
subplot(3,2,3)
plot(vel6_newB{minIDbar})
hold on
for s=1:length(startpoints)
    plot([startpoints(s) startpoints(s)],[0 max(vel6_newB{maxIDbar})],'-r')
end
axis([15000 45000 0 Inf])
subplot(3,2,5)
curr=flipud(reshape(normPeaks(minIDdg,:,2),4,6));
[x,y] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr,[0 1]);
colormap('gray')
subplot(3,2,2)
plot([1 2 4 8 16],blackAll(maxIDbar,2:end)/max(blackAll(maxIDbar,2:end)),'-ok');
hold on
plot([1 2 4 8 16],newdata(maxIDdg,:)/max(newdata(maxIDdg,:)),'-ob');
axis([0 17 0 1])
title(['cell ',int2str(maxIDbar), ' (',int2str(maxIDdg),') / ',num2str(allval(doubleBar(maxID))),', ',num2str(DGspeed(maxIDdg))])
subplot(3,2,4)
plot(vel6_newB{maxIDbar})
hold on
for s=1:length(startpoints)
    plot([startpoints(s) startpoints(s)],[0 max(vel6_newB{maxIDbar})],'-r')
end
axis([15000 45000 0 Inf])
subplot(3,2,6)
curr=flipud(reshape(normPeaks(maxIDdg,:,2),4,6));
[x,y] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr,[0 1]);
colormap('gray')
set(gcf,'position',[1601 1 1600 824])
saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_speed_',int2str(minIDbar),'_',int2str(maxIDbar),'.emf']))



[~,sortID2]=sort(DGspeed);
minID=132;
maxID=87;

figure
subplot(2,2,1)
plot([1 2 4 8 16],newdata(minID,:)/max(newdata(minID,:)),'-ob');
axis([0 17 0 1])
title(['cell ',int2str(minID),' / ',num2str(DGspeed(minID))])
subplot(2,2,3)
curr=flipud(reshape(normPeaks(minID,:,2),4,6));
[x,y] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr,[0 1]);
colormap('gray')
subplot(2,2,2)
plot([1 2 4 8 16],newdata(maxID,:)/max(newdata(maxID,:)),'-ob');
axis([0 17 0 1])
title(num2str(DGspeed(maxID)))
subplot(2,2,4)
curr=flipud(reshape(normPeaks(maxID,:,2),4,6));
[x,y] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr,[0 1]);
colormap('gray')
set(gcf,'position',[1601 1 1600 824])
saveas(gcf,fullfile('S:\Paper factory\Katja\5_HUMPIG\newFigs',['ex_speed_',int2str(minID),'_',int2str(maxID),'.emf']))