clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
load(fullfile(path1,'w_responding'))

%% flash over experiment, step10

load(fullfile(path1,['w_info']))
good=find(goodones==1);

rDates=info(good,1);
rUnits=info(good,2);

cd('R:\Users\Katja\human_retina\_scripts')

backgrRate=zeros(length(good),60);
for g=1:length(good)
    
    currUnit=rUnits{g};
    currDate=rDates{g};
    
    
    
    dpath=fullfile(path1,'WT',currDate);
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
    
    file=file(find(file>412));
    file=file(11:end);
    
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
        backgrRate(g,cnt)=bckgr;
    end
    
    if add==1
        cnt=cnt+1;
        curr=allrate(endR+5:end,:);
        curr=mean(curr);
        bckgr=mean(curr(200:1900));
        backgrRate(g,cnt)=bckgr;
    end
    backgrRate(g,cnt+1:end)=-1;
end
keep=backgrRate;
backgrRate(find(isnan(backgrRate)))=-1;

%%
close all
normPeaks=zeros(length(backgrRate),15);
for i=1:length(backgrRate)
    curr=backgrRate(i,:);
    tmp=find(curr>-1);
    tmp2=find(curr==-1);
    curr=curr(tmp);
    curr=curr/max(curr);
    normPeaks(i,tmp)=curr;
    normPeaks(i,tmp2)=-1;
end

forhist=zeros(15,2);
figure
for i=1:15
    curr=backgrRate(:,i);
    tmp=find(curr>-1);
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
% forhist(i,2)=length(tmp2);
end
axis([0 16 -5 50])

figure
bar(forhist)


figure
for i=1:15
    curr=normPeaks(:,i);
    tmp=find(curr>-1);
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
axis([0 16 -0.1 1.1])

