clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
load(fullfile(path1,'h_responding'))

%%

load(fullfile(path1,['h_info']))

good=find(goodones==1);
LF=cell2mat(info(good,9));

lfOFF=find(LF==-1);lfOFF=good(lfOFF);
lfON=find(LF==1);lfON=good(lfON);

load(fullfile(path1,['h_LFs']))

figure
for of=1:length(lfOFF)
    curr=LFs(:,lfOFF(of));
    subplot(2,1,1)
    plot(curr)
    hold on
    subplot(2,1,2)
    plot((curr-min(curr))/(max(curr)-min(curr)))
    hold on
end


figure
for of=1:length(lfON)
    curr=LFs(:,lfON(of));
    subplot(2,1,1)
    plot(curr)
    hold on
    subplot(2,1,2)
    plot((curr-min(curr))/(max(curr)-min(curr)))
    hold on
end


%%
close all
load(fullfile(path1,['h_info']))

good=find(goodones==1);
LF=cell2mat(info(good,9));

lfOFF=find(LF==-1);lfOFF=good(lfOFF);
lfON=find(LF==1);lfON=good(lfON);

load(fullfile(path1,['h_LFs']))

latencies=zeros(length(lfOFF)+length(lfON),3);
for of=1:length(lfOFF)
    latencies(of,1)=lfOFF(of);
%     figure
    curr=LFs(:,lfOFF(of));
    data=curr(1:210);
    dataS=[];
    for d=1:7:203
        dataS=[dataS;mean(data(d:d+6))];
    end
    [val s]=min(dataS);
%     plot(1:7:203,dataS)
%     hold on
%     plot(s,val,'o')
    latencies(of,2)=-1;
    latencies(of,3)=s*7;
end
cnt=of;

for of=1:length(lfON)
    latencies(of+cnt,1)=lfON(of);
%     figure
    curr=LFs(:,lfON(of));
    data=curr(1:210);
    dataS=[];
    for d=1:7:203
        dataS=[dataS;mean(data(d:d+6))];
    end
    [val s]=max(dataS(5:25));
%     plot(dataS)
%     hold on
%     plot(s+1,val,'o')
    latencies(of+cnt,2)=1;
    latencies(of+cnt,3)=(s+4)*7;
end

barnum=10;
figure
hist(latencies(:,3),barnum)
figure
tmp=find(latencies(:,2)==-1);
hist(latencies(tmp,3),barnum)
figure
tmp=find(latencies(:,2)==1);
hist(latencies(tmp,3),barnum)

%%
[a b]=sort(latencies(:,3));
cols='kbcgrkbcgr';
fa='kbcgrwwwww';
k=5;
ma='o*';

figure
for i=1:length(b)
    if latencies(b(i),2)==1
        type=1;
    else
        type=2;
    end
    plot(i,a(i),ma(type),'color',cols(ID(b(i),k)),'markerfacecolor',fa(ID(b(i),k)))
    hold on
end




