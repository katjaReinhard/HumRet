%%
clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
load(fullfile(path1,'h_responding'))

%%
load(fullfile(path1,['h_info']))
good=find(goodones==1);
% rDates=info(resp,1);
% rUnits=info(resp,2);
rDates=info(good,1);
rUnits=info(good,2);

cd('R:\Users\Katja\human_retina\_scripts')

cnt=0;
exnone=0;
exwith=0;
% figure
spontRate=zeros(length(good),20);
for g=1:length(good)
    
    cnt=cnt+1;
    
    currUnit=rUnits{g};
    currDate=rDates{g};
    
    %     startF=info{resp(g),4};
    %     stopF=info{resp(g),5};
    startF=info{good(g),4};
    stopF=info{good(g),5};
    
    
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
        if ~isempty(regexp(heka{h},'spont'))
            file=[file;h];
        end
    end
    %     tmp=find(file>=startF);
    %     file=file(tmp);
    %     tmp=find(file<=stopF);
    %     file=file(tmp);
    
    
    allrate=[];
    for r=1:length(file)
        spikes=unit{1,2}{file(r),2};
        spikes=round(spikes);
        rate=zeros(1,65000);
        rate(spikes)=1;
        rate(spikes(diff(spikes)==0))=2;
        allrate=[allrate rate(1:60000)];
        rate2=MEA_spikerates(spikes,10,65000);
        spontRate(g,r)=mean(rate2(1:60000));
    end
    spontRate(g,r+1:end)=-1;
    
    
    meanrate=allrate;
    
    
    
    
    Y=fft(meanrate,15000)/numel(meanrate);
    f=1000/2*linspace(0,1,15000/2+1);
    c=2*abs(Y(1:15000/2+1));
    c=c';
    
    tmp1=find(f>0.8);
    tmp1=tmp1(1);
    tmp3=find(f<=15);
    tmp3=tmp3(end);
    f=f(tmp1:tmp3);
    c=c(tmp1:tmp3);
    
    
    bckr=mean(c);
    bckrSTD=std(c);
    thr=bckr+5*bckrSTD;
    tmp=find(c>thr);
    if ~isempty(tmp)
        oscillations(g,1)=1;
%         if exwith<3
%             exwith=exwith+1;
           figure
           plot(f,c)
           title('osc')
            axis([0 15 0 0.0007])
%         end
    else
        oscillations(g,1)=0;
        if exnone<3 && std(c)>0
            exnone=exnone+1;
           figure
           plot(f,c)
           title('none')
           axis([0 15 0 0.0007])
        end
    end
    oscillations(g,2)=good(g);
    if ismember(resp,good(g))
        oscillations(g,3)=1;
    else
        oscillations(g,3)=0;
    end
    goodDates{g}=currDate;
    
    %      if cnt>16
    %          figure
    %          cnt=1;
    %      end
    %      subplot(4,4,cnt)
    %      plot(f,c)
    %      hold on
    %      plot([0 32],[thr thr],'-r')
    
    
end

%% plot development of firing rate
% figure
spont=zeros(length(spontRate),20);
for i=1:length(spontRate)
    curr=spontRate(i,:);
    some=find(curr>-1);
    other=find(curr==-1);
    %     some=1:13;
    if ~isempty(some)
        data=curr(some);
        norm=data/max(data);
        tmp=find(isnan(norm));
        norm(tmp)=-1;
        spont(i,some)=norm;
        spont(i,other)=-1;
        % plot(some,norm)
        % hold on
    end
end


for i=1:13
    curr=spont(:,i);
    tmp=find(curr>-1);
    %     tmp=1:length(curr);
    data=curr(tmp);
    figure(1)
    plot(i,mean(data),'ok','markerfacecolor','k','markersize',10)
    hold on
    plot(i,median(data),'oc','markerfacecolor','c','markersize',10)
    plot([i i],[mean(data)-std(data) mean(data)+std(data)],'-k','linewidth',5)
    plot([i i],[min(data) max(data)],'-g','linewidth',2)
    text(i,2,int2str(length(tmp)))
end
axis([0 14 -0.5 2.3])



figure
for i=1:13
    curr=spontRate(:,i);
    tmp=find(curr>-1);
    %     tmp=1:length(curr);
    data=curr(tmp);
    plot(i,mean(data),'ok','markerfacecolor','k','markersize',10)
    hold on
    plot(i,median(data),'oc','markerfacecolor','c','markersize',10)
    plot([i i],[mean(data)-std(data) mean(data)+std(data)],'-k','linewidth',5)
    plot([i i],[min(data) max(data)],'-g','linewidth',2)
    %     text(i,2,int2str(length(tmp)))
end
% axis([0 14 -0.5 2.3])

