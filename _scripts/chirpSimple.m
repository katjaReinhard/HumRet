%% FFT chirp for various sigmas and various frequency resolutions (NFFT)
% this script calculates the FFT (fast fourier transform) of the chirp
% stimulus and the chirp response. 

% I'm assuming that you still use our standard HEKA file format

plotting=0;

sigmas=[5 10 20 40]; %sigmas for firing rate convolution
NFFTs=[2000 4000 6000 8000 10000 15000]; %resoultion of FFT


co='kbcgrm';

%--------this is specific for our data format-----
% load(fullfile(superpath,[Dtype,'_info']))
% % load(fullfile(superpath,[Dtype,'_chirp']))
% cd(superpath)
% savepath=fullfile(superpath,'chirp_FFTnew');
% 
% togo=[];
% for i=1:length(info)
%     if ~isempty(info{i,19})
%         if info{i,19}>0
%             togo=[togo;i];
%         end
%     end
% end
%-------------------------------



% +++ variable "togo" contains identifiers for all cells you want to
% analyze +++

chirpFFTabs=cell(length(togo),2,2);
chirpFFTsmooth=cell(length(togo),2,2);
chirpRates=cell(length(togo),2);
for i=1:length(togo)
    i
    currID=togo(i);
    
    % +++ since you don't have "info" you might have to change this part
    currd=info{currID,1}; % date of experiment
    curru=info{currID,2}; % name of the cell you want to analyze
    currstart=info{currID,4}; % first file that you want to analyze
    currstop=info{currID,5}; % last file that you want to analyze
    
    hpath=fullfile(mainpath,currd,'HEKA'); %here are the heka files
    hlist=dir(fullfile(hpath,'*phys')); %make list with all heka files
    upath=fullfile(mainpath,currd,'units');%here is the data (if it's in the heka files then you have to write that yourself)
    ulist=dir(fullfile(upath,'*mat')); %list of data files
    % ++++++++++++++++++++++++
    
    % find the cell you want to analyse
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,['unit_',curru]))
            found=1;
        end
    end
    load(fullfile(upath,ulist(m).name)) %load data file
    
    % +++ CHANGE TO YOUR DATA FORMAT: list of all stimuli that have been shwon
    hekalist=unit{1,2}(:,1);
    %+++++++++++++++++
    
    %go through list of stimuli and find all that contain "chirp"
    chirp=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'chirp%'))
            ftype=1;
            chirp=[chirp;1];
        else
            chirp=[chirp;0];
        end
    end
    %select only the files within the range you set above
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(chirp==1);
    
    %load protocol using "read_header_field_heka.m"
    prot=read_header_field_heka(hpath,hekalist{goodfiles(1)},'Stimulus Protocol');
   
        
    if ~isempty(goodfiles)
        
        sicnt=0; %try different sigmas for firing rate convolution
        for si=1:5
            sicnt=sicnt+1;
            if si<5
                sigma=sigmas(si);
                
                %calculate firing rate with the current sigma
                allrate=[];
                for f=goodfiles
                    spikes=unit{1,2}{f,2};
                    rate=MEA_spikerates(spikes,sigma,25000);
                    allrate=[allrate;rate];
                end
                meanrate=mean(allrate);
            else
                % no convolution, simply take spike times
                allrate=[];
                for f=goodfiles
                    spikes=unit{1,2}{f,2};
                    spikes=round(spikes);
                    rate=zeros(1,25000);
                    rate(spikes)=1;
                    rate(spikes(diff(spikes)==0))=2;
                    allrate=[allrate;rate];
                end
                meanrate=mean(allrate);
            end
            % save firing rates
            chirpRates{i,sicnt}=meanrate;
            
            % check with your chirp files, but first 2 lines of protocol
            % should just be background settings. otherwise change
            int=prot(3:end-1,4);
            
            %find the frequency part of the chirp
            time=prot(2:end-2,1);
            timediff=diff(time);
            firstpart=find(timediff>50);
            if ~isempty(firstpart)
                firstpart=firstpart(1);
            else
                firstpart=length(timediff);
            end
            rate=meanrate(round(time(1)):round(time(firstpart)));
            
            %make a vector that contains for each millisecond the gray
            %value of the chirp stimulus
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
            
            %go through different FFT resolutions
            ncnt=0;
            for NFFT=[2000 4000 6000 8000 10000 15000]
                ncnt=ncnt+1;
                %FFT of stimulus
                Y=fft(int2,NFFT)/numel(int2);
                fS=1000/2*linspace(0,1,NFFT/2+1);
                cS=2*abs(Y(1:NFFT/2+1));
                
                %FFT of response
                Y=fft(rate,NFFT)/numel(rate);
                f=1000/2*linspace(0,1,NFFT/2+1);
                c=2*abs(Y(1:NFFT/2+1));
                c=c';
                
                %cut frequency space (you don't want to analyze super low
                %frequencies as the FFT causes weird artefacts, nor
                %frequencies higher than your stimulus
                tmp1=find(f>0.8);
                tmp1=tmp1(1);
                tmp2=find(f<32);
                tmp2=tmp2(end);
                tmp3=find(f<=7.5);
                tmp3=tmp3(end);
                %
                %overview plots of the different sigmas and NFFTs
                if plotting==1
                    %figure of FFTs of response
                    figure(1)
                    subplot(2,6,ncnt+6)
                    plot(f(tmp1:tmp3),c(tmp1:tmp3),'color',co(sicnt))
                    hold on
                    subplot(2,6,sicnt)
                    plot(f(tmp1:tmp3),c(tmp1:tmp3),'color',co(ncnt))
                    hold on
                    
                    figure(2)
                     %figure of FFTs of stimulus
                    subplot(2,6,ncnt+6)
                    plot(fS(tmp1:tmp3),cS(tmp1:tmp3),'color',co(sicnt))
                    hold on
                    subplot(2,6,sicnt)
                    plot(fS(tmp1:tmp3),cS(tmp1:tmp3),'color',co(ncnt))
                    hold on
                end
                
                %calculate ratio of FFT response/stimulus
                c=c(tmp1:tmp3); cS=cS(tmp1:tmp3);
                cc=c/max(c); ccS=cS/max(cS); ratio=cc./ccS;
                
                chirpFFTabs{i,sicnt+1,ncnt}=c;
                chirpFFTabs{i,7,ncnt}=cS;
                chirpFFTabs{i,1,ncnt}=f(tmp1:tmp3);
                
                %do some smoothing that makes comparison easier
                window=3; step=1;
                smoo=[]; smooF=[];
                for to=1:step:length(ratio)-step
                    if to+window>length(ratio)
                        endval=length(ratio);
                    else
                        endval=to+window;
                    end
                    me=mean(ratio(to:endval));
                    smoo=[smoo;me];
                    fr=mean(f(to+tmp1-1:endval+tmp1-1));
                    smooF=[smooF;fr];
                end
                
                %plot the smoothed results for each sigma and NFFT
                if plotting==1
                    figure(3)
                    subplot(2,6,ncnt+6)
                    plot(smooF,smoo,'color',co(sicnt))
                    hold on
                    
                    subplot(2,6,sicnt)
                    plot(smooF,smoo,'color',co(ncnt))
                    hold on
                end
                
                chirpFFTsmooth{i,sicnt+1,ncnt}=smoo;
                chirpFFTsmooth{i,1,ncnt}=smooF;
            end
        end
        
        %add titles to figures and save them
        if plotting==1
            figure(1)
            for sup=1:5
                subplot(2,6,sup)
                axis([0 8 0 Inf])
                if sup<5
                    title(['sigma=',int2str(sigmas(sup)),' / all NFFTs (kbcgrm)'])
                else
                    title(['binary / all NFFTs (kbcgrm)'])
                end
            end
            for sup=1:6
                subplot(2,6,sup+6)
                axis([0 8 0 Inf])
                title(['all sigma (kbcgr) / NFFT=',int2str(NFFTs(sup))])
            end
            set(gcf,'position',[1 31 1600 794])
            saveas(gcf,fullfile(superpath,'newAttempt','chirp',[Dtype,'_',int2str(currID),'_responseFFT.bmp']))
            close
            figure(2)
            for sup=1:5
                subplot(2,6,sup)
                axis([0 8 0 Inf])
                if sup<5
                    title(['sigma=',int2str(sigmas(sup)),' / all NFFTs (kbcgrm)'])
                else
                    title(['binary / all NFFTs (kbcgrm)'])
                end
            end
            for sup=1:6
                subplot(2,6,sup+6)
                axis([0 8 0 Inf])
                title(['all sigma (kbcgrm) / NFFT=',int2str(NFFTs(sup))])
            end
            set(gcf,'position',[1 31 1600 794])
            saveas(gcf,fullfile(superpath,'newAttempt','chirp',[Dtype,'_',int2str(currID),'_stimulusFFT.bmp']))
            close
            figure(3)
            for sup=1:5
                subplot(2,6,sup)
                axis([0 8 0 Inf])
                if sup<5
                    title(['sigma=',int2str(sigmas(sup)),' / all NFFTs (kbcgrm)'])
                else
                    title(['binary / all NFFTs (kbcgrm)'])
                end
            end
            for sup=1:6
                subplot(2,6,sup+6)
                axis([0 8 0 Inf])
                title(['all sigma (kbcgrm) / NFFT=',int2str(NFFTs(sup))])
            end
            set(gcf,'position',[1 31 1600 794])
            saveas(gcf,fullfile(superpath,'newAttempt','chirp',[Dtype,'_',int2str(currID),'_respDivByStimFFT.bmp']))
            close
        end
    end
    
end
%save necessary variables
save(fullfile(superpath,'newAttempt2','chirp',[Dtype,'_chirpFFT']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates')

