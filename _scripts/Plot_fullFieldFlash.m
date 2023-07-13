%% EXAMPLE CODE: PLOTS EXAMPLE CELL RESPONSES LIKE IN FIGURE 3
% code for data from Reinhard & Münch 2021

%% VARIABLES THAT NEED TO BE CHANGED
%PATH TO DOWNLOADED SCRIPTS (including MEA_spikerates.m and read_header_field_heka.m
path2code = 'C:\Users\Lenovo2020\OneDrive - imec\HumRet\_scripts';

% PATH TO DOWNLOADED VARIABLES from folder HumanRGC_processedData
path1 = 'C:\Users\Lenovo2020\OneDrive - imec\HumRet';

% PATH TO STIMULUS INFO
% set this to any HEKA folder within a data folder (downloaded from HumanRGC_spiketimes)
% e.g. HumanRGC_spiketimes\20130211c\HEKA
path2Heka = 'C:\Users\Lenovo2020\OneDrive - imec\HumRet\data\20130211c\HEKA';
% name of any HEKA files with your desired stimulus in that folder
% FFFashBW for full field flashes as in Figure 1G
file2read = '20130211_C1#0022_HumanPig_FFFlashBW%dur_2000_gap_2000%_FW0ND4.phys';




%% load variable indicating responding units AND SELECTED CELLS TO PLOT

% rDates and rUnits will contain the date and unit number that has been
% classified as "responding"
cd(path2code)
load(fullfile(path1,'h_responding'))
load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);

% >>>>>> SELECT CELLS TO PLOT
% based on rDates and rUnits, select cells to plot
% recording date, unit number, e.g.:
toplotDate={'20120224','20120224'};
toplotUnit=[35 65];

%% plot flash responses

% read stimulus information
prot = read_header_field_heka(path2Heka, file2read,'Stimulus Protocol');
% column 1: end times of the gray values defined in column 2 and 3 (or start time of the values defined in the next row)
% column 2: 0 if background value should be used (=128), 2 if value in column 3 should be used
% column 3: gray value to be used if column 2 = 2

start=prot(2:end,1);
val = prot(2:end,4);
idx = prot(2:end,3);
tmp = find(idx==0);
val(tmp)=128;
% !! start(1) is the END TIME of val(1) or the START TIME of val(2)

figure
for g=1:length(toplotUnit) % go through chosen units
    
    % --- find IDs of chosen unit ---
    numU = toplotUnit(g);
    dateU = toplotDate{g};
    currUnit=strcat('000',num2str(numU));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,dateU);
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    
    % --- find start and end file in light level that was analyzed ---
    startF=info{ID2,4};
    stopF=info{ID2,5};
    
    % --- paths to data for this unit ---
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
    
    % --- find files that contain chosen stimulus ---
    heka=unit{1,2}(:,1);
    go_files=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlashBW'))
            go_files=[go_files;h];
        end
    end
    tmp=find(go_files>startF);
    go_files=go_files(tmp);
    tmp=find(go_files<=stopF);
    go_files=go_files(tmp);
    
    if ~isempty(go_files)
  
        % collect spike times and compute firing rate
        allrate=[];
        for f=1:length(go_files)
            spikes=unit{1,2}{go_files(f),2};
            rate=MEA_spikerates(spikes,10,23000);
            rate = rate(1:start(end));
            allrate=[allrate;rate];
        end
        meanrate=mean(allrate,1);
        
        % --- plot the firing rate ---
        subplot(ceil(length(toplotUnit)/2),2,g)
        plot(meanrate,'k')
        % alternative way to plot (like in Figure 3)
        %                 patch([1:length(meanrate) length(meanrate) 1],[meanrate 0 0],'k','edgecolor','k')
        hold on
        
        % --- plot stimulus changes ---
        for st = 1:length(start)-2
            plot([start(st) start(st)],[0 max(meanrate)],'-','color',[0.5 0.5 0.5])
            text(start(st),max(meanrate)*1.01,['begin of ',int2str(val(st+1))],'fontsize',6)
        end
        
        ylabel('firing rate (spks/sec)')
        xlabel('time (ms)')
        
     
        title([dateU,' - ',strcat('unit_',currUnit)],'interpreter','none')
        
    end
end