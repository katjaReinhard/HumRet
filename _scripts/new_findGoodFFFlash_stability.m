clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
load(fullfile(path1,'h_responding'))

%%
close all
load(fullfile(path1,['h_info']))

% toplotDate={'20120224','20130510a'};
% toplotUnit=[35 100];

toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];



rDates=info(resp,1);
rUnits=info(resp,2);

cd('R:\Users\Katja\human_retina\_scripts')


for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    startF=info{ID2,4};
    stopF=info{ID2,5};
    
    dpath=fullfile(path1,'Human',toplotDate{g});
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
    
    
    
    
    figure
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        spikes=spikes(find(spikes<6000));
        for s=1:length(spikes)
            plot([spikes(s) spikes(s)],[f-0.4 f+0.4],'-b')
            hold on
            
        end
        
        
        
    end
    title(ulist(m).name,'interpreter','none')
end

%%
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


for g=41:50
    
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
        if ~isempty(regexp(heka{h},'FFFlash'))
            file=[file;h];
        end
    end
    
    
    
    
    figure
    for f=1:length(file)
        spikes=unit{1,2}{file(f),2};
        spikes=spikes(find(spikes<6000));
        for s=1:length(spikes)
            plot([spikes(s) spikes(s)],[f-0.4 f+0.4],'-b')
            hold on
            
        end
        
        
        
    end
    title(ulist(m).name,'interpreter','none')
end