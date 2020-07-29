% clear
% close all
% path1='R:\Users\Katja\human_retina';
path1 = 'C:\Users\katja\OneDrive - imec\HumRet';

%---this part is just an ID-list of cells where I want to analyse chirp---%
load(fullfile(path1,['h_info']))
togo=[];
for i=1:length(info)
    if info{i,19}==1
        togo=[togo;i];
    end
end
togo = [togo(33);togo(1:32);togo(34:end)];

%---rDates contains dates of chosen units, rUnits contains unit number---%
rDates=info(togo,1);
rUnits=info(togo,2);

% cd('R:\Users\Katja\human_retina\_scripts')
cd('C:\Users\katja\OneDrive - imec\HumRet\_scripts')


for g=1
    
    currUnit=rUnits{g};
    currDate=rDates{g};
    
    startF=info{togo(g),4};
    stopF=info{togo(g),5};
    
    %pathes%
    dpath=fullfile(path1,'data',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    %load current unit
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    %find files with chirp stimulus
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
        prot=read_header_field_heka(fullfile(dpath,'HEKA'),unit{1,2}{flashfiles(3),1},'Stimulus Protocol');
        start=prot(2,1);
        start=start-200;
        start2=prot(485,1);
        
        int=prot(3:end-1,4);
        
        %----first part of chirp
        time=prot(2:end-2,1);
        timediff=diff(time);
        firstpart=find(timediff>50);
        if ~isempty(firstpart)
            firstpart=firstpart(1);
        else
            firstpart=length(timediff);
        end
        
        su=sum(timediff);
        timediff=round(timediff);
        roundsum=sum(timediff);
        changing=ceil(roundsum-su);
        stepping=floor(length(timediff)/changing);
        for ti=1:stepping:length(timediff)
            timediff(ti)=16;
        end
        
        int2=[]; idds=[];
        for ii=1:firstpart-1
            tmp=repmat(int(ii),round(timediff(ii)),1);
            int2=[int2;tmp];
            idds=[idds;round(timediff(ii))];
        end
        missing=8000-length(int2);
        NFFT=9000;
        Y=fft(int2,NFFT)/numel(int2);
        fS=1000/2*linspace(0,1,NFFT/2+1);
        cS=2*abs(Y(1:NFFT/2+1));
        
        figure('position',[987 1311 700 70])
        rectangle('position',[1 0 length(int2) 255],'facecolor','w','edgecolor','none')
        hold on
        plot(int2,'-k','linewidth',2)
        axis off
        box off
                axis tight

        exportgraphics(gcf,fullfile('C:\Users\katja\OneDrive - imec\HumRet\PlosOne\stim','Stim4.png'),'resolution',600)
        
        
        
        
        
        %---contrast
        int=prot(485:end-1,4);
        int=[128; int];
        
        time=prot(2:end-2,1);
        time=time-start2;
        tmp=find(time<=19000);
        time=time(tmp);
        tmp=find(time>0);
        time=time(tmp);
        timediff=diff(time);
        
        
        
        su=sum(timediff);
        timediff=round(timediff);
        roundsum=sum(timediff);
        changing=ceil(roundsum-su);
        stepping=floor(length(timediff)/changing);
        for ti=1:stepping:length(timediff)+2
            timediff(ti)=16;
        end
        tmp=find(timediff>40);
        timediff(tmp)=16;
        
        int2=[]; idds=[];
        for ii=1:length(timediff)-1
            tmp=repmat(int(ii),round(timediff(ii)),1);
            int2=[int2;tmp];
            idds=[idds;round(timediff(ii))];
        end
        int2=[repmat(128,round(time(1)),1); int2];
        int2=[int2;repmat(128,9000-length(int2),1)];
        
        figure('position',[987 1311 700 70])
                rectangle('position',[1 0 length(int2) 255],'facecolor','w','edgecolor','none')

                hold on
        plot(int2(1:8100),'-k','linewidth',2)
        axis off
        axis tight
        box off
        exportgraphics(gcf,fullfile('C:\Users\katja\OneDrive - imec\HumRet\PlosOne\stim','Stim3.png'),'resolution',600)
        
        
    end
end
