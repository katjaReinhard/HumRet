clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
path1='C:\Users\KatjaR\OneDrive - imec\HumRet';
load(fullfile(path1,'h_responding'))

%% flash over experiment, step10

load(fullfile(path1,['h_info']))
good=find(goodones==1);
good2=[];
for i=1:length(good)
    if info{good(i),8}~=0
        good2=[good2;good(i)];
    end
end

rDates=info(good2,1);
rUnits=info(good2,2);

cd('R:\Users\Katja\human_retina\_scripts')

peaks=zeros(length(good2),60);
for g=1:length(good2)
    
    currUnit=rUnits{g};
    currDate=rDates{g};
    
    
    
    dpath=fullfile(path1,'Human',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    file=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlashBW%dur_2000_gap_2000%_'))
            file=[file;h];
        end
    end
    
    di=diff(file);
    tmp=find(di>1);
    if tmp(1)<5
        file=file(tmp(1)+1:end);
    end
    
    
    
    allrate=[];
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        rate=MEA_spikerates(spikes,10,11000);
        rate=rate(1:10000);
        allrate=[allrate;rate];
        
    end
    
    
    if mod(f,10)>5
        endR=f-mod(f,10);
        add=1;
    elseif mod(f,10)>0
        endR=f-mod(f,10);
        add=0;
    else
        endR=f;
        add=0;
    end
    
    cnt=0;
    for j=1:10:endR
        cnt=cnt+1;
        if j==endR
            if add==0
                curr=allrate(j:end,:);
            else
                curr=allrate(j:j+9,:);
            end
        else
            curr=allrate(j:j+9,:);
        end
        curr=mean(curr);
        bckgr=mean(curr(200:1900));
        curr=curr-bckgr;
        ma=max(curr);
        peaks(g,cnt)=ma;
    end
    if add==1
        cnt=cnt+1;
        curr=allrate(endR+10:end,:);
        curr=mean(curr);
        bckgr=mean(curr(200:1900));
        curr=curr-bckgr;
        ma=max(curr);
        peaks(g,cnt)=ma;
    end
    peaks(g,cnt+1:end)=-1;
    
    
end


peaks(find(isnan(peaks)))=0;


%%

normPeaks=zeros(length(peaks),17);
for i=1:length(peaks)
    curr=peaks(i,:);
    tmp=find(curr>0);
    tmp2=find(curr==-1);
    curr=curr(tmp);
    curr=curr/max(curr);
    normPeaks(i,tmp)=curr;
    normPeaks(i,tmp2)=-1;
end


figure
for i=1:16
    curr=peaks(:,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
    plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
    text(i,200,[int2str(length(tmp)),'/',int2str(length(tmp2))])
end
axis([0 17 0 210])


figure
for i=1:16
    curr=normPeaks(:,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
    plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
    text(i,1.3,[int2str(length(tmp)),'/',int2str(length(tmp2))])
end
axis([0 17 0 1.5])


%% flash over experiment, step5

load(fullfile(path1,['h_info']))
good=find(goodones==1);
good2=[];
for i=1:length(good)
    if info{good(i),8}~=0
        good2=[good2;good(i)];
    end
end

rDates=info(good2,1);
rUnits=info(good2,2);

cd('R:\Users\Katja\human_retina\_scripts')

peaks=zeros(length(good2),60);
for g=1:length(good2)
    
    currUnit=rUnits{g};
    currDate=rDates{g};
    
    
    
    dpath=fullfile(path1,'Human',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    file=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlashBW%dur_2000_gap_2000%_'))
            file=[file;h];
        end
    end
    
    di=diff(file);
    tmp=find(di>1);
    if tmp(1)<5
        file=file(tmp(1)+1:end);
    end
    
    
    
    allrate=[];
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        rate=MEA_spikerates(spikes,10,11000);
        rate=rate(1:10000);
        allrate=[allrate;rate];
        
    end
    
    
    if mod(f,5)>3
        endR=f-mod(f,5);
        add=1;
    elseif mod(f,5)>0
        endR=f-mod(f,5);
        add=0;
    else
        endR=f;
        add=0;
    end
    
    cnt=0;
    for j=1:5:endR
        cnt=cnt+1;
        if j==endR
            if add==0
                curr=allrate(j:end,:);
            else
                curr=allrate(j:j+4,:);
            end
        else
            curr=allrate(j:j+4,:);
        end
        curr=mean(curr);
        bckgr=mean(curr(200:1900));
        stback=std(curr(200:1900));
        curr=curr-bckgr;
        ma=max(curr);
        if ma>1.5*bckgr+stback
            peaks(g,cnt)=ma;
        end
    end
    if add==1
        cnt=cnt+1;
        curr=allrate(endR+5:end,:);
        curr=mean(curr);
        bckgr=mean(curr(200:1900));
        stback=std(curr(200:1900));
        curr=curr-bckgr;
        ma=max(curr);
        if ma>1.5*bckgr+stback
            peaks(g,cnt)=ma;
        end
    end
    peaks(g,cnt+1:end)=-1;
    
    
end


peaks(find(isnan(peaks)))=0;


%%

normPeaks=zeros(length(peaks),32);
for i=1:length(peaks)
    curr=peaks(i,:);
    tmp=find(curr>0);
    tmp2=find(curr==-1);
    curr=curr(tmp);
    curr=curr/max(curr);
    normPeaks(i,tmp)=curr;
    normPeaks(i,tmp2)=-1;
end

forhist=zeros(32,2);
figure
for i=1:32
    curr=peaks(:,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
    plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
%     text(i,260,[int2str(length(tmp2))])
%     text(i,265,[int2str(length(tmp))])
forhist(i,1)=length(tmp);
forhist(i,2)=length(tmp2);
end
axis([0 33 -10 270])

figure
bar(forhist)


figure
for i=1:32
    curr=normPeaks(:,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
    plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
%     text(i,1.3,[int2str(length(tmp2))])
%     text(i,1.5,[int2str(length(tmp))])
end
axis([0 33 0 1.7])

%% flash over experiment, step5, new version

load(fullfile(path1,['h_info']))
good=find(goodones==1);
good2=[];
for i=1:length(good)
    if info{good(i),8}~=0
        good2=[good2;good(i)];
    end
end

rDates=info(good2,1);
rUnits=info(good2,2);

cd('R:\Users\Katja\human_retina\_scripts')

peaks=zeros(length(good2),60);
for g=1:length(good2)
    
    currUnit=rUnits{g};
    currDate=rDates{g};
    
    
    
    dpath=fullfile(path1,'Human',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    file=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlashBW%dur_2000_gap_2000%_'))
            file=[file;h];
        end
    end
    
    di=diff(file);
    tmp=find(di>1);
    if tmp(1)<5
        file=file(tmp(1)+1:end);
    end
    
    
    
    allrate=[];
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        rate=MEA_spikerates(spikes,10,11000);
        rate=rate(1:10000);
        allrate=[allrate;rate];
        
    end
    
    
    if mod(f,5)>3
        endR=f-mod(f,5);
        add=1;
    elseif mod(f,5)>0
        endR=f-mod(f,5);
        add=0;
    else
        endR=f;
        add=0;
    end
    
    cnt=0;
    for j=1:5:endR
        cnt=cnt+1;
        if j==endR
            if add==0
                curr=allrate(j:end,:);
            else
                curr=allrate(j:j+4,:);
            end
        else
            curr=allrate(j:j+4,:);
        end
        curr=mean(curr);
        bckgr=mean(curr(200:1900));
        stback=std(curr(200:1900));
        curr=curr-bckgr;
        ma=max(curr);
        if ma>1.5*bckgr+stback
            peaks(g,cnt)=ma;
        end
    end
    if add==1
        cnt=cnt+1;
        curr=allrate(endR+5:end,:);
        curr=mean(curr);
        bckgr=mean(curr(200:1900));
        stback=std(curr(200:1900));
        curr=curr-bckgr;
        ma=max(curr);
        if ma>1.5*bckgr+stback
            peaks(g,cnt)=ma;
        end
    end
    peaks(g,cnt+1:end)=-1;
    tst=peaks(g,:);
    tmp=find(tst>0);
    if tmp(1)>1
        peaks(g,1:tmp(1)-1)=-1;
    end
    
end


peaks(find(isnan(peaks)))=0;


%%

normPeaks=zeros(length(peaks),32);
for i=1:length(peaks)
    curr=peaks(i,:);
    tmp=find(curr>0);
    tmp2=find(curr==-1);
    curr=curr(tmp);
    curr=curr/max(curr);
    normPeaks(i,tmp)=curr;
    normPeaks(i,tmp2)=-1;
end

forhist=zeros(32,2);
figure
for i=1:32
    curr=peaks(:,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
    plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
%     text(i,260,[int2str(length(tmp2))])
%     text(i,265,[int2str(length(tmp))])
forhist(i,1)=length(tmp);
forhist(i,2)=length(tmp2);
end
axis([0 33 -10 270])

figure
bar(forhist)


figure
for i=1:32
    curr=normPeaks(:,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
    plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
%     text(i,1.3,[int2str(length(tmp2))])
%     text(i,1.5,[int2str(length(tmp))])
end
axis([0 33 0 1.7])


%% flash over experiment, step5, EJ firing rate

load(fullfile(path1,['h_info']))
good=find(goodones==1);
good2=[];
for i=1:length(good)
    if info{good(i),8}~=0
        good2=[good2;good(i)];
    end
end

rDates=info(good2,1);
rUnits=info(good2,2);

% cd('R:\Users\Katja\human_retina\_scripts')
cd('C:\Users\KatjaR\OneDrive - imec\HumRet\_scripts')

peaks=zeros(length(good2),60);
backgrounds=zeros(length(good2),60);
for g=6
% for g=1:length(good2)
%     
    currUnit=rUnits{g};
    currDate=rDates{g};
    
    
    
    dpath=fullfile(path1,'data',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    file=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlashBW%dur_2000_gap_2000%_'))
            file=[file;h];
        end
    end
    
    di=diff(file);
    tmp=find(di>1);
    if tmp(1)<5
        file=file(tmp(1)+1:end);
    end
    
    
    
    sigma = 20;
    flength = 5*sigma;
    x = linspace(-flength / 2, flength / 2, flength);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    
    %     figure
    allrate=[];allrate2=[];
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        spikes=spikes(find(spikes<10000));
        spikes=round(spikes*10);
        spikes=spikes(spikes>0);
        rate=zeros(1,100000);
        rate(spikes)=1;
        rate(spikes(diff(spikes)==0))=2;
                rate2=MEA_spikerates(spikes/10,60,6000);
        allrate=[allrate;rate];
                allrate2=[allrate2;rate2];
    end
    
    
    
    if mod(f,5)>3
        endR=f-mod(f,5);
        add=1;
    elseif mod(f,5)>0
        endR=f-mod(f,5);
        add=0;
    else
        endR=f;
        add=0;
    end
    
    newrate = []; newrate2 = []; newrateB = zeros(32,size(allrate,2)/10);
    cnt=0;
    for j=1:5:endR
        cnt=cnt+1;
        if j==endR
            if add==0
                curr=allrate(j:end,:);
                curr2=allrate2(j:end,:);
            else
                curr=allrate(j:j+4,:);
                 curr2=allrate2(j:j+4,:);
            end
        else
            curr=allrate(j:j+4,:);
               curr2=allrate2(j:j+4,:);
        end
        now = mean(curr); cntr = 0;
        for w = 1:10:size(curr,2)
            cntr = cntr+1;
           if ~isempty(find(curr(w:w+9))>0)
              newrateB(cnt,cntr)=1;
           end
        end
        newrate2=[newrate2;mean(curr2)];
        curr = filter (gaussFilter,1, mean(curr));
        newrate = [newrate;curr];
        bckgr=mean(curr(2000:19000));
        stback=std(curr(2000:19000));
        curr=curr-bckgr;
        ma=max(curr);
        if ma>3.5*bckgr+stback
            peaks(g,cnt)=ma;
            backgrounds(g,cnt)=mean(bckgr);
        end
    end
%     figure
cnt = 0;
for c = 13:25
    cnt = cnt+1;
    figure
%     subplot(6,6,c)
    plot(newrate2(c,:))
    hold on
    plot([2000 2000],[0 20],'-k')
    plot([4000 4000],[0 20],'-k')
    plot([1800 6000],[0 0],'-k')
    axis([1800 6000 0 21])
    axis off
    saveas(gcf,fullfile('C:\Users\KatjaR\OneDrive - imec\HumRet\Figs',['final_exampleFFFlash_',int2str(cnt),'.emf']))
    close
end

% figure
% for x=1:32
%     for y=1:size(newrateB,2)
%         if newrateB(x,y)==1
%         plot([y y],[x-0.3 x+0.3],'-k')
%         hold on
%         end
%     end
% end
    
    if add==1
        cnt=cnt+1;
        curr=allrate(endR+5:end,:);
        curr = filter (gaussFilter,1, mean(curr));
        bckgr=mean(curr(2000:19000));
        stback=std(curr(2000:19000));
        curr=curr-bckgr;
        ma=max(curr);
        if ma>3.5*bckgr+stback
            peaks(g,cnt)=ma;
        end
    end
    peaks(g,cnt+1:end)=-1;
    tst=peaks(g,:);
    tmp=find(tst>0);
    if tmp(1)>1
        peaks(g,1:tmp(1)-1)=-1;
    end
    
end


peaks(find(isnan(peaks)))=0;
peaks=peaks*10000;
backgrounds=backgrounds*10000;
%%
figure
for c = 1:size(newrate,1)
    subplot(6,6,c)
    plot(newrate2(c,:))
    axis([0 6000 0 12])
end
%%

close all
normPeaks=zeros(length(peaks),32);
for i=1:length(peaks)
    curr=peaks(i,:);
    tmp=find(curr>0);
    tmp2=find(curr==-1);
    curr=curr(tmp);
    curr=curr/max(curr);
    normPeaks(i,tmp)=curr;
    normPeaks(i,tmp2)=-1;
end

forhist=zeros(13,2);
figure
for i=1:13
    curr=peaks(take,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
%     plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
%     text(i,260,[int2str(length(tmp2))])
%     text(i,265,[int2str(length(tmp))])
forhist(i,1)=length(tmp);
forhist(i,2)=length(tmp2);
end
axis([0 14 0 140])

figure
bar(forhist)


figure
for i=1:14
    curr=normPeaks(take,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
%     plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
%     text(i,1.3,[int2str(length(tmp2))])
%     text(i,1.5,[int2str(length(tmp))])
end
axis([0 33 0 1.7])


%% flash
close all



load(fullfile(path1,['h_info']))
good=find(goodones==1);
good2=[];
for i=1:length(good)
    if info{good(i),8}~=0
        good2=[good2;good(i)];
    end
end

rDates=info(good2,1);
rUnits=info(good2,2);

cd('R:\Users\Katja\human_retina\_scripts')

peaks=zeros(length(good2),1);




for g=1:length(rDates)
    g
     currUnit=rUnits{g};
    currDate=rDates{g};
    
    ID2=good2(g);
    startF=info{ID2,4};
    stopF=info{ID2,5};
    
    dpath=fullfile(path1,'data',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    file=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlash'))
            file=[file;h];
        end
    end
     tmp=find(file>startF);
    file=file(tmp);
    tmp=find(file<=stopF);
    file=file(tmp);
    
    sigma = 20;
    flength = 5*sigma;
    x = linspace(-flength / 2, flength / 2, flength);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    
    %     figure
    allrate=[];
   for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        spikes=spikes(find(spikes<10000));
        spikes=round(spikes*10);
        spikes=spikes(spikes>0);
        rate=zeros(1,100000);
        rate(spikes)=1;
        rate(spikes(diff(spikes)==0))=2;
        %         rate2=MEA_spikerates(spikes/10,40,6000);
        allrate=[allrate;rate];
        %         allrate2=[allrate2;rate2];
    end
    yfilt = filter (gaussFilter,1, mean(allrate));
    peaks(g)=max(yfilt)*10000;
end

figure
plot([1 2],[mean(peaks) mean(peaks)],'-k')
hold on
plot([1 2],[median(peaks) median(peaks)],'-c')
plot([1.5 1.5],[mean(peaks)-std(peaks) mean(peaks)+std(peaks)],'-k','linewidth',3)
plot([1.5 1.5],[min(peaks) max(peaks)],'-k')
    

figure
hist(peaks)
figure
hist(peaks,20)
%% development only for retinas that have been recorded for a long time




normPeaks=zeros(length(peaks),32);
take=[]; takeDates=cell(2,1);
cnt=0;
for i=1:length(peaks)
    curr=peaks(i,:);
    if curr(14)~=-10000
        cnt=cnt+1;
        take=[take;i];
        takeDates{cnt}=rDates{i};
    end
    tmp=find(curr>0);
    tmp2=find(curr==-1);
    curr=curr(tmp);
    curr=curr/max(curr);
    normPeaks(i,tmp)=curr;
    normPeaks(i,tmp2)=-1;
end

close all


forhist=zeros(13,2);
figure
for i=1:13
    curr=peaks(take,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'ok')
    hold on
    plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
%     plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
%     text(i,260,[int2str(length(tmp2))])
%     text(i,265,[int2str(length(tmp))])
forhist(i,1)=length(tmp);
forhist(i,2)=length(tmp2);
end
axis([0 14 0 140])

figure
bar(forhist)



figure
for i=1:13
    curr=backgrounds(take,i);
    tmp=find(curr>0);
    tmp2=find(curr>-1);
    dat=curr(tmp);
    plot(i,mean(dat),'og')
    hold on
%     plot(i,median(dat),'oc')
    plot([i i],[mean(dat)-std(dat) mean(dat)+std(dat)],'-k','linewidth',3)
%     plot([i i],[min(dat) max(dat)],'-k','linewidth',1)
%     text(i,260,[int2str(length(tmp2))])
%     text(i,265,[int2str(length(tmp))])

end
axis([0 14 -5 20])


