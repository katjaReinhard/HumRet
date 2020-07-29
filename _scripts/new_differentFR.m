clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
load(fullfile(path1,'h_responding'))
%% flash
close all



load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);


toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];

cd('R:\Users\Katja\human_retina\_scripts')

for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    currDate=toplotDate{g};
    
    cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    startF=info{ID2,4};
    stopF=info{ID2,5};
    
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
    allrate=[];allrate2=[];
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        spikes=spikes(find(spikes<6000));
        spikes=round(spikes*10);
        spikes=spikes(spikes>0);
        rate=zeros(1,60000);
        rate(spikes)=1;
        rate(spikes(diff(spikes)==0))=2;
        %         rate2=MEA_spikerates(spikes/10,40,6000);
        allrate=[allrate;rate];
        %         allrate2=[allrate2;rate2];
    end
    yfilt = filter (gaussFilter,1, mean(allrate));
    %     yfilt2 = conv (mean(allrate), gaussFilter, 'same');
figure;plot(0.1:0.1:length(yfilt)/10,yfilt*10000)    %      figure;plot(yfilt2*10000)
    %     figure;plot(mean(allrate)*10000)
    
    title(ulist(m).name,'interpreter','none')
end

%% chirp
close all



load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);


toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];

cd('R:\Users\Katja\human_retina\_scripts')

for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    currDate=toplotDate{g};
    
     cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    startF=info{ID2,4};
    stopF=info{ID2,5};
    
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
        if ~isempty(regexp(heka{h},'chirp'))
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
    allrate=[];allrate2=[];
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        spikes=spikes(find(spikes<12000));
        spikes=round(spikes*10);
        spikes=spikes(spikes>0);
        rate=zeros(1,120000);
        rate(spikes)=1;
        rate(spikes(diff(spikes)==0))=2;
        %         rate2=MEA_spikerates(spikes/10,40,6000);
        allrate=[allrate;rate];
        %         allrate2=[allrate2;rate2];
    end
    yfilt = filter (gaussFilter,1, mean(allrate));
    %     yfilt2 = conv (mean(allrate), gaussFilter, 'same');
    figure;plot(yfilt*10000)
    %      figure;plot(yfilt2*10000)
    %     figure;plot(mean(allrate)*10000)
    
    title(ulist(m).name,'interpreter','none')
end
%% DG: dynamic range
close all

 load(fullfile(path1,'newAttempt','DG','h_normPeaks'))

load(fullfile(path1,['h_info']))
rDates=info(togo,1);
rUnits=info(togo,2);

dglist={'0100um%dur_12000_speed_1';'0100um%dur_12000_speed_2';'0100um%dur_12000_speed_4';'0100um%dur_12000_speed_8';'0200um%dur_12000_speed_1';'0200um%dur_12000_speed_2';'0200um%dur_12000_speed_4';'0200um%dur_12000_speed_8';'0500um%dur_12000_speed_1';'0500um%dur_12000_speed_2';'0500um%dur_12000_speed_4';'0500um%dur_12000_speed_8';'1000um%dur_12000_speed_1%';'1000um%dur_12000_speed_2';'1000um%dur_12000_speed_4';'1000um%dur_12000_speed_8';'2000um%dur_12000_speed_1';'2000um%dur_12000_speed_2';'2000um%dur_12000_speed_4';'2000um%dur_12000_speed_8';'4000um%dur_12000_speed_1';'4000um%dur_12000_speed_2';'4000um%dur_12000_speed_4';'4000um%dur_12000_speed_8'};
spats=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
frs=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];

cd('R:\Users\Katja\human_retina\_scripts')

% for g=1:length(togo)
for g=50:70
    
    currUnit=strcat('000',num2str(rUnits{g}));
    currUnit=currUnit(end-3:end);
    currDate=rDates{g};
    
   
    startF=info{togo(g),4};
    stopF=info{togo(g),5};
    
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
    
    [~, best]=max(normPeaks(g,:,2));
    
    heka=unit{1,2}(:,1);
    file=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},dglist{best}))
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
    allrate=[];allrate2=[];
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        spikes=spikes(find(spikes<14000));
%         spikes=spikes(find(spikes>2000));
        spikes=round(spikes*10);
        spikes=spikes(spikes>0);
        rate=zeros(1,140000);
        rate(spikes)=1;
        rate(spikes(diff(spikes)==0))=2;
%         rate=rate(20000:end);
                rate2=MEA_spikerates(spikes/10,40,14000);
        allrate=[allrate;rate];
                allrate2=[allrate2;rate2(1:14000)];
    end
    yfilt = filter (gaussFilter,1, mean(allrate));
  
    currTemp=frs(best);
%     amp=max(yfilt*10000);
%     ttt=1/1000:1/1000:12;
%     sinewave=amp/2*sin(2*pi*currTemp*ttt);
%     
    figure
    subplot(1,2,1)
    plot(yfilt*10000)
     subplot(1,2,2)
    plot(mean(allrate2))
%     hold on
%     plot(20010:10:140000,sinewave+amp/2,'-r')

%     cycles=12*currTemp;
%     intervals=120000/cycles;
%     
%     for cc=20000:intervals:140000
%        curr=yfilt(cc:cc+intervals-1)*10000; 
%        ma=max(curr);
%        mi=min(curr);
%        
%     end
    
end
