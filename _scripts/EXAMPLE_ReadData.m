%% EXAMPLE CODE: PLOTS EXAMPLE CELL RESPONSES LIKE IN FIGURE 3
% code for data from Reinhard & Münch 2021

%% load variable indicating responding units
path2code = ' '; %PATH TO DOWNLOADED SCRIPTS
cd(path2code)

path1 = ' ';% PATH TO DOWNLOADED VARIABLES
load(fullfile(path1,'h_responding'))
load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);

%% select cells to plot
% recording date, unit number

toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];

%% plot chirp responses

% --- read stimulus data ---
prot = read_header_field_heka('C:\Users\Lenovo2020\OneDrive - imec\HumRet\data\20130211c\HEKA',...
    '20130211_C1#0283_HumanPig_chirp%mean_128_amplitude_128_pause_500%_FW0ND4.phys','Stimulus Protocol');
start=prot(2,1);% start of frequency part
start=start-200; 
start2=prot(485,1); % start of contrast part

figure
for g=1:length(toplotUnit)
    
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
        if ~isempty(regexp(heka{h},'chirp'))
            go_files=[go_files;h];
        end
    end
    tmp=find(go_files>startF);
    go_files=go_files(tmp);
    tmp=find(go_files<=stopF);
    go_files=go_files(tmp);
    
    if ~isempty(go_files)
      
        % === contrast
        allrate=[];
        for f=1:length(go_files)
            spikes=unit{1,2}{go_files(f),2};
            rate=MEA_spikerates(spikes,10,23000);
            rate=rate(start2:start2+8500);
            allrate=[allrate;rate];
        end
        
        meanrate=mean(allrate,1);
        subplot(17,2,g*2-1)
        patch([1:length(meanrate) length(meanrate) 1],[meanrate 0 0],'k','edgecolor','k')
        axis([0 8500 0 inf])
        axis tight
        axis off
        box off
        
        % === freq
        allrate=[];
        for f=1:length(go_files)
            spikes=unit{1,2}{go_files(f),2};
            rate=MEA_spikerates(spikes,10,13000);
            rate=rate(start:10000+start);
            allrate=[allrate;rate];
        end
        meanrate=mean(allrate,1);
        subplot(17,2,g*2)
        patch([1:length(meanrate) length(meanrate) 1],[meanrate 0 0],'k','edgecolor','k')
        hold on
        axis([0 10000 0 inf])
        axis tight
        axis off
        box off
        
    end
end