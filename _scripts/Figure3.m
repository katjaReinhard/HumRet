clear
close all

%% pathes etc
path1='C:\Users\katja\OneDrive - imec\HumRet';
pathcode = 'C:\Users\katja\OneDrive - imec\HumRet\_scripts';
pathColor = 'F:\LabCode\Matlab\Plots\othercolor';
pathSave = 'C:\Users\katja\OneDrive - imec\HumanRetina\PlosOne';
pathPanel = 'C:\Users\katja\OneDrive - imec\HumanRetina\PlosOne\stim';

clrs = [166,206,227;31,120,180;178,223,138;51,160,44;251,154,153;227,26,28;...
    253,191,111;255,127,0;202,178,214;106,61,154];
clrs = clrs./255;


%% load
addpath(genpath(pathcode))
addpath(genpath(pathColor))
load(fullfile(path1,'h_responding'))
load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);
load(fullfile(path1,'EXAMPLE_INFO'))

fsize1 = 7; fsize2=9; fsize3=10;
%% choose units
toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];
order = [17 16 15 4 7 5 1 2 3 6 8 14 9 13 10];
%% axes
XX = [ 0.0500    0.2100    0.3700    0.5300    0.6900    0.805];
YY = [0.06:0.059:0.97];
H = 0.056;
W = [0.145 0.145 0.145 0.145 0.1 0.19];
breaksat = [4 11];
YY(breaksat) = YY(breaksat)+0.01;

%% start figure
FIG = figure('position',[635.6667 821 810.0000 1140],'color','w')

save_png = 1;
%% stimuli
for a = 1:6
    if a <5
        axes('position',[XX(a) 0.96 W(a) 0.03])
    elseif a == 5
        axes('position',[0.71  0.96 0.06 0.03])
    elseif a == 6
        axes('position',[0.84  0.96 0.09 0.03])
    end
    set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
    I = imread(fullfile(pathPanel,['Stim',int2str(a),'.png']));
    imagesc(I)
    % axis equal
    axis off
    box off
    if a == 3 || a==4
        axis tight
    end
    
end
%% plot example cells: flash
load(fullfile(path1,'h_flashSpikes'))


for g=1:length(order)
    numU = toplotUnit(order(g));
    dateU = toplotDate{order(g)};
    currUnit=strcat('000',num2str(numU));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,dateU);
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    
    data=flash_raster(ID2,:);
    
    axes('position',[XX(1) YY(g) W(1) H/3-0.005])
    allrate=[];
    for f=1:length(data)
        spikes=data{f};
        rate=MEA_spikerates(spikes,40,15000);
        rate=rate(200:10500);
        allrate=[allrate;rate];
        spikes=spikes(find(spikes>=200));
        spikes=spikes-200+1;
        spikes=spikes(find(spikes<=5800));
        if ~isempty(spikes)
            for s=1:length(spikes)
                plot([spikes(s) spikes(s)],[f-0.4 f+0.4],'-k')
                hold on
                
            end
        end
        
    end
    plot([1800 1800],[0 f+1],'-','color',[0.6 0.6 0.6])
    plot([3800 3800],[0 f+1],'-','color',[0.6 0.6 0.6])
    axis([0 5800 0 inf])
    axis tight
    axis off
    box off
    
    axes('position',[XX(1) YY(g)+H/3 W(1) H/3])
    %
    meanrate=mean(allrate,1); meanrate = meanrate(1:5800);
    patch([1:5800 5800 1],[meanrate 0 0],'k','edgecolor','k')
    %     plot(meanrate,'-k')
    axis([0 5800 0 inf])
    hold on
    ma=max(meanrate);
    plot([1800 1800],[0 ma],'-','color',[0.6 0.6 0.6])
    plot([3800 3800],[0 ma],'-','color',[0.6 0.6 0.6])
    axis off
    box off
end

%% bar
noresp = [1 6 8 10 11 12 14 ];
load(fullfile(path1,'h_vel6_new_Black'))

prot=read_header_field_heka('C:\Users\katja\OneDrive - imec\HumRet\data\20120224\HEKA','20120224_C1#0014_OPA1_bar_6veloc_downup_black%speed_200_gap_3000%_FW0ND4.phys','Stimulus Protocol');

cons=prot(1:end,2);
changes=find(cons==-1);
startpoints=round(prot(changes(1:7),1))-20500;
startpoints=startpoints(3:end-1);

for g=1:length(order)
    numU = toplotUnit(order(g));
    dateU = toplotDate{order(g)};
    currUnit=strcat('000',num2str(numU));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,dateU);
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    
    
    %     data=vel6_newB{ID2};
    
    %----alternative
    startF=info{ID2,4};
    stopF=info{ID2,5};
    
    
    dpath=fullfile(path1,'data',dateU);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    clear unit
    load(fullfile(upath,ulist(m).name))
    
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'bar_6veloc_downup_black%'))
            flashfiles=[flashfiles;h];
        end
    end
    tmp=find(flashfiles>startF);
    flashfiles=flashfiles(tmp);
    tmp=find(flashfiles<=stopF);
    flashfiles=flashfiles(tmp);
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(spikes,40,44100);
        %             rate=rate(start2:start2+8500);
        allrate=[allrate;rate];
    end
    %----
    data = mean(allrate);
    
    
    axes('position',[XX(2) YY(g) W(1) H-0.02])
    dat=data(20500:end-2000);
    if ~isempty(find(noresp==g))
        patch([1:length(dat) length(dat) 1],[dat 0 0],[0.4 0.4 0.4],'edgecolor',[0.4 0.4 0.4])
    else
        patch([1:length(dat) length(dat) 1],[dat 0 0],'k','edgecolor','k')
    end
    hold on
    
    for s=1:length(startpoints)
        plot([startpoints(s) startpoints(s)],[0 max(dat)],'-','color',[0.6 0.6 0.6])
    end
    bar_startpoints = startpoints;
    bar_end = length(dat);
    
    axis([0 inf 0 inf])
    axis off
    box off
end
%% chirp

prot = read_header_field_heka('C:\Users\katja\OneDrive - imec\HumRet\data\20130211c\HEKA',...
    '20130211_C1#0283_HumanPig_chirp%mean_128_amplitude_128_pause_500%_FW0ND4.phys','Stimulus Protocol');
start=prot(2,1);
start=start-200;
start2=prot(485,1);
for g=1:length(order)
    numU = toplotUnit(order(g));
    dateU = toplotDate{order(g)};
    currUnit=strcat('000',num2str(numU));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,dateU);
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    
    startF=info{ID2,4};
    stopF=info{ID2,5};
    
    
    dpath=fullfile(path1,'data',dateU);
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
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'chirp'))
            flashfiles=[flashfiles;h];
        end
    end
    tmp=find(flashfiles>startF);
    flashfiles=flashfiles(tmp);
    tmp=find(flashfiles<=stopF);
    flashfiles=flashfiles(tmp);
    
    if ~isempty(flashfiles)
        %         prot=read_header_field_heka(fullfile(dpath,'HEKA'),unit{1,2}{flashfiles(1),1},'Stimulus Protocol');
        %         start=prot(2,1);
        %         start=start-200;
        %         start2=prot(485,1);
        
        
        % === contrast
        allrate=[];
        for f=1:length(flashfiles)
            spikes=unit{1,2}{flashfiles(f),2};
            rate=MEA_spikerates(spikes,10,23000);
            rate=rate(start2:start2+8500);
            allrate=[allrate;rate];
        end
        
        meanrate=mean(allrate,1);
        axes('position',[XX(3) YY(g) W(1) H-0.02])
        patch([1:length(meanrate) length(meanrate) 1],[meanrate 0 0],'k','edgecolor','k')
        axis([0 8500 0 inf])
        axis tight
        axis off
        box off
        
        % === freq
        allrate=[];
        for f=1:length(flashfiles)
            spikes=unit{1,2}{flashfiles(f),2};
            rate=MEA_spikerates(spikes,10,13000);
            rate=rate(start:10000+start);
            allrate=[allrate;rate];
        end
        meanrate=mean(allrate,1);
        axes('position',[XX(4) YY(g) W(1) H-0.02])
        patch([1:length(meanrate) length(meanrate) 1],[meanrate 0 0],'k','edgecolor','k')
        hold on
        axis([0 10000 0 inf])
        axis tight
        axis off
        box off
        
    end
end

%% DG responses
% noresp = [1 6 8 10 11 12 14 ];
% load(fullfile(path1,'h_DGrates'))
load(fullfile(path1,'h_normPeaks'))

spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];

COLL_ID = [];
for g=1:length(order)
    numU = toplotUnit(order(g));
    dateU = toplotDate{order(g)};
    currUnit=strcat('000',num2str(numU));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,dateU);
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    COLL_ID = [COLL_ID;ID2];
    ID3=find(togo==ID2);
    if ~isempty(ID3)
        [dg id]=max(normPeaks(ID3,:,2));
        
        frnow = freq(id); spnow = spat(id);
        stimname = ['0' num2str(spnow)];
        stimname = stimname(end-3:end);
        stimname = [stimname 'um%dur_12000_speed_' int2str(frnow)];%0200um%dur_12000_speed_1%
        
        %     data=vel6_newB{ID2};
        
        %----alternative
        startF=info{ID2,4};
        stopF=info{ID2,5};
        
        
        dpath=fullfile(path1,'data',dateU);
        upath=fullfile(dpath,'units');
        ulist=dir(fullfile(upath,'*mat'));
        
        found=0;m=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
                found=1;
            end
        end
        
        clear unit
        load(fullfile(upath,ulist(m).name))
        
        
        
        heka=unit{1,2}(:,1);
        flashfiles=[];
        for h=1:length(heka)
            if ~isempty(regexp(heka{h},stimname))
                flashfiles=[flashfiles;h];
            end
        end
        tmp=find(flashfiles>startF);
        flashfiles=flashfiles(tmp);
        tmp=find(flashfiles<=stopF);
        flashfiles=flashfiles(tmp);
        allrate=[];
        for f=1:length(flashfiles)
            spikes=unit{1,2}{flashfiles(f),2};
            rate=MEA_spikerates(spikes,40,20000);
            %             rate=rate(start2:start2+8500);
            allrate=[allrate;rate];
        end
        %----
        data = mean(allrate);
        
        
        axes('position',[XX(6) YY(g) W(6) H-0.02])
        dat=data(1:15000);
        
        patch([1:length(dat) length(dat) 1],[dat 0 0],'k','edgecolor','none')
        
        hold on
        %
        plot([2000 2000],[0 max(dat)],'-','color',[0.6 0.6 0.6])
        %
        axis([0 inf 0 inf])
        axis off
        box off
        
        axes('position',[XX(5) YY(g) W(5) H-0.02])
        htmp = normPeaks(ID3,:,2);
        htmp = reshape(htmp,4,6);
        imagesc(fliplr(flipud(1-htmp)))
        colormap('gray')
        hold on
        rectangle('position',[0.5 0.5 6 4],'edgecolor','k','facecolor','none','linewidth',0.5)
        axis equal
        
        if g == 2
            axis off
            spl = [4000 2000 1000 500 200 100];
            frl = [8 4 2 1];
            hold on
            for xa = 1:6
                text(xa,5,int2str(spl(xa)),'rotation',90,'fontname','Arial','fontsize',fsize1,'HorizontalAlignment','right')
            end
            for xa = 1:4
                text(0.3,xa,int2str(frl(xa)),'fontname','Arial','fontsize',fsize1,'HorizontalAlignment','right')
            end
            
            text(3.5,8,'\mum','HorizontalAlignment','center','fontname','Arial','fontsize',fsize1)
            text(1,-0.3,'Hz','fontname','Arial','fontsize',fsize1,'HorizontalAlignment','right')
        else
            axis off
        end
    else
        %         id = 24;
    end
end
%% lines
for aa = 2:length(XX)
    axes('position',[XX(aa)-0.01 0.04 0.01 0.95])
    plot([0 0],[0 1],'-k','linewidth',1.5)
    axis tight
    axis off
    box off
end

for aa = 1:length(breaksat)
    axes('position',[0.03 YY(breaksat(aa))-0.022 0.965 0.01])
    plot([0 1],[0 0],'-k','linewidth',1.5)
    axis tight
    axis off
    box off
end

axes('position',[0.03 YY(2)-0.011 0.655 0.01])
 plot([0 1],[0 0],'-','linewidth',0.5,'color',[0.4 0.4 0.4])
    axis tight
    axis off
    box off
    
 

for aa = 3:length(YY)-1
    if isempty(find(breaksat==aa))
    axes('position',[0.03 YY(aa)-0.011 0.965 0.01])
    plot([0 1],[0 0],'-','linewidth',0.5,'color',[0.4 0.4 0.4])
    axis tight
    axis off
    box off
    end
end

axes('position',[0.03 0.945 0.965 0.01])
plot([0 1],[0 0],'-k','linewidth',1.8)
axis tight
axis off
box off
%% scale  bars
axes('position',[XX(1) YY(1)-0.015 W(1) 0.009])
axis off
box off
plot([1 2],[1 1],'-k','linewidth',3)
text(1.5,0.2,'2 seconds','horizontalalignment','center','fontsize',fsize1,'fontname','Arial')
axis([1 4 0 1])
box off
axis off
box off

axes('position',[XX(2) YY(1)-0.015 W(2) 0.009])
axis off
box off
barlist = [1 2 4 8 16];
for s = 1:length(bar_startpoints)
    plot([bar_startpoints(s) bar_startpoints(s)],[0 1],'-','color',[0.4 0.4 0.4],'linewidth',1)
    hold on
    text(bar_startpoints(s)+2000,0.4,int2str(barlist(s+1)),'horizontalalignment','center','fontsize',fsize1,'fontname','Arial')
end
    text(bar_startpoints(1)/2,0.4,[int2str(barlist(1)) 'mm/s'],'horizontalalignment','center','fontsize',fsize1,'fontname','Arial')
axis([1 bar_end 0 1])
axis off
box off

axes('position',[XX(3) YY(1)-0.015 W(3) 0.009])
axis off
box off
plot([0 2],[1 1],'-k','linewidth',3)
text(1.7,0.2,'2 seconds','horizontalalignment','center','fontsize',fsize1,'fontname','Arial')
axis([0 8 0 1])
axis off
box off

axes('position',[XX(4) YY(1)-0.015 W(4) 0.009])
axis off
box off
plot([0 2],[1 1],'-k','linewidth',3)
text(1.8,0.2,'2 seconds','horizontalalignment','center','fontsize',fsize1,'fontname','Arial')
axis([0 10 0 1])
axis off
box off

axes('position',[XX(6) YY(2)-0.015 W(6) 0.009])
axis off
box off
plot([2 4],[1 1],'-k','linewidth',3)
text(3,0.2,'2 seconds','horizontalalignment','center','fontsize',fsize1,'fontname','Arial')
axis([0 14 0 1])
axis off
box off
%% letters
ex_names = {'A_{1}','A_{2}','A_{3}','A_{4}','A_{5}',...
    'B_{1}','B_{2}','B_{3}','B_{4}','B_{5}','B_{6}','B_{7}','C_{1}','C_{2}','C_{3}'};
% ex_names = {'A1','A2','A3','A4','A5','B1','B2','B3','B4','B5','B6','B7','C1','C2','C3'};

ex_names = fliplr(ex_names);
for ii = 1:length(COLL_ID)
    inow = COLL_ID(ii);
    iddg = find(SUPER_COLL(:,1)==inow);
    pos =  YY(ii)+H-0.02;
    if ii == 4 || ii == 11
        pos = pos-0.01;
    end
    if ~isempty(iddg)
        cln = SUPER_COLL(iddg,2);
        axes('position',[0.025 pos  0.01 0.01])
        text(0,0,ex_names{ii},'fontsize',fsize3,'fontweight','bold','fontname','Arial','color',clrs(cln,:))
        set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
        box off
    else
        axes('position',[0.025 pos  0.01 0.01])
        text(0,0,ex_names{ii},'fontsize',fsize3,'fontweight','bold','fontname','Arial','color','k')
        set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
        box off
    end
end
%% save
exportgraphics(gcf,fullfile(pathSave,'NEW_Fig3_examples.png'),'Resolution',600)
exportgraphics(gcf,fullfile(pathSave,'NEW_Fig3_examples.tif'),'Resolution',600)
% saveas(gcf,fullfile(pathSave,'NEW_Fig3_examples.tif'))