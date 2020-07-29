clear
close all

path1='S:\data\Katja\CMOS_MEA\forSaad';
date='20131210c_ntk';
conf='config2';

ufolder=fullfile(path1,date,conf,'units');
ulist=dir(fullfile(ufolder,'*mat'));

if ~exist(fullfile(path1,date,conf,'pics'),'dir')
    mkdir(fullfile(path1,date,conf,'pics'))
end

stimuli={'FFFlash','DS_alldirections_black','DS_alldirections_white','chirp','bar_6veloc_downup_black','bar_6veloc_downup_white'};
stimLength=[12000 45000 45000 25000 45000 45000];
take=1:6;

for u=1:length(ulist)
    load(fullfile(ufolder,ulist(u).name))
    
    number=ulist(u).name(end-7:end-4);
    if ~exist(fullfile(path1,date,conf,'pics',number),'dir')
        mkdir(fullfile(path1,date,conf,'pics',number))
    end
    
    for t=1:length(take)
        currstim=stimuli{take(t)};
        
         if ~exist(fullfile(path1,date,conf,'pics',currstim),'dir')
        mkdir(fullfile(path1,date,conf,'pics',currstim))
    end
        
        
        files=[];
        for f=1:length(unit{1,2})
            if ~isempty(regexp(unit{1,2}{f,1},currstim))
                files=[files;f];
            end
        end
        
        allrate=[];
        for x=1:length(files)
            spikes=unit{1,2}{files(x),2};
            rate=MEA_spikerates(spikes,20,stimLength(take(t)));
            allrate=[allrate;rate];
        end
        meanrate=mean(allrate,1);
        figure
        plot(meanrate)
        title(ulist(u).name,'interpreter','none')
        set(gcf,'position',[1601 1 1600 824])
        saveas(gcf,fullfile(path1,date,conf,'pics',number,[currstim,'.emf']))
           saveas(gcf,fullfile(path1,date,conf,'pics',currstim,[number,'.emf']))
        close
    end
end
    
    
    