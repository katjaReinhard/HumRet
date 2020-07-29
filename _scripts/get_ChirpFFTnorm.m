function [chirpFFTnorm,chirpAmpsmoothed,smooF,togo]=get_ChirpFFTnorm(path1,info)
warning('off')
%---this part is just an ID-list of cells where I want to analyse chirp---%
% load(fullfile(path1,['h_info']))
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


chirpFFTsmoothed=zeros(length(togo),63);
chirpAmpsmoothed=zeros(length(togo),8501);


for g=1:length(togo)
    
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
        if g==1
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
            
        end
        
        allrate=[];
        for f=1:length(flashfiles)
            spikes=unit{1,2}{flashfiles(f),2};
            rate=MEA_spikerates(round(spikes),10,13000);
            rate=rate(start:10000+start);
            allrate=[allrate;rate];
        end
        meanrate=mean(allrate,1);
        
        
        Y=fft(meanrate,NFFT)/numel(meanrate);
        freq=1000/2*linspace(0,1,NFFT/2+1);
        data=2*abs(Y(1:NFFT/2+1));
        data=data';
        
        tmp1=find(freq>=0.35);
        tmp2=find(freq<=7.5);
        %         freq=freq(tmp1(1):tmp2(end));
        
        if g==1
           
             cS=cS(tmp1(1):tmp2(end)); ccS=cS/max(cS); 
        end
        
        c=data(tmp1(1):tmp2(end));
        cc=c/max(c);  ratio=cc./ccS;
        
        window=5; step=1;
        smoo=[]; smooF=[];
        for to=1:step:length(ratio)-step
            if to+window>length(ratio)
                endval=length(ratio);
            else
                endval=to+window;
            end
            me=mean(ratio(to:endval));
            smoo=[smoo;me];
            fr=mean(freq(to+tmp1(1)-1:endval+tmp1(1)-1));
            smooF=[smooF;fr];
        end
        chirpFFTsmoothed(g,:)=smoo;
        
        
        
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
        
        
        allrate=[];
        for f=1:length(flashfiles)
            spikes=unit{1,2}{flashfiles(f),2};
            rate=MEA_spikerates(spikes,40,23000);
            rate=rate(start2:start2+8500);
            allrate=[allrate;rate];
        end
        
        meanrate=mean(allrate,1);
        
        bckgr = mean(meanrate(1:500));
        meanrate = meanrate-bckgr;
        meanrate = meanrate./max(meanrate);
      
        chirpAmpsmoothed(g,:)=meanrate;
        
        
        
    end
end

%

chirpFFTnorm=[];
for c=1:size(chirpFFTsmoothed,1)
    curr=chirpFFTsmoothed(c,:);
    norm=curr/max(curr);
    chirpFFTnorm=[chirpFFTnorm;norm];
    
    
end
warning('on')