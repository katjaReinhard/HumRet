%% FFT chirp for various sigmas and various frequency resolutions (NFFT)

plotting=0;

sigmas=[5 10 20 40];
NFFTs=[2000 4000 6000 8000 10000 15000];


co='kbcgrm';

load(fullfile(superpath,[Dtype,'_info']))
% load(fullfile(superpath,[Dtype,'_chirp']))
cd(superpath)
savepath=fullfile(superpath,'chirp_FFTnew');

togo=[];
for i=1:length(info)
    if ~isempty(info{i,19})
        if info{i,19}>0
            togo=[togo;i];
        end
    end
end

chirpFFTabs=cell(length(togo),2,2);
chirpFFTsmooth=cell(length(togo),2,2);
chirpRates=cell(length(togo),2);
for i=1:length(togo)
    i
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    currstart=info{currID,4};
    currstop=info{currID,5};
    
    hpath=fullfile(mainpath,currd,'HEKA');
    hlist=dir(fullfile(hpath,'*phys'));
    upath=fullfile(mainpath,currd,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,['unit_',curru]))
            found=1;
        end
    end
    load(fullfile(upath,ulist(m).name))
    
    hekalist=unit{1,2}(:,1);
    
    chirp=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'chirp%'))
            ftype=1;
            chirp=[chirp;1];
        else
            chirp=[chirp;0];
        end
    end
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(chirp==1);
    
    if ~isempty(find(Dtype=='w'))
       goodfiles=goodfiles(6:end); 
    end
    
    prot=read_header_field_heka(hpath,hekalist{goodfiles(1)},'Stimulus Protocol');
    if prot(end,1)>24000
        goodfiles=[];
    end
    
    
    if ~isempty(goodfiles)
        
        sicnt=0;
        for si=1:5
            sicnt=sicnt+1;
            if si<5
                sigma=sigmas(si);
                
                
                allrate=[];
                for f=goodfiles
                    spikes=unit{1,2}{f,2};
                    rate=MEA_spikerates(spikes,sigma,25000);
                    allrate=[allrate;rate];
                end
                meanrate=mean(allrate);
            else
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
            
            chirpRates{i,sicnt}=meanrate;
            
            
            int=prot(3:end-1,4);
            
            time=prot(2:end-2,1);
            timediff=diff(time);
            firstpart=find(timediff>50);
            if ~isempty(firstpart)
                firstpart=firstpart(1);
            else
                firstpart=length(timediff);
            end
            rate=meanrate(round(time(1)):round(time(firstpart)));
            
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
            
            ncnt=0;
            for NFFT=[2000 4000 6000 8000 10000 15000]
                ncnt=ncnt+1;
                %             NFFT=8000;
                Y=fft(int2,NFFT)/numel(int2);
                fS=1000/2*linspace(0,1,NFFT/2+1);
                cS=2*abs(Y(1:NFFT/2+1));
                
                
                
                %             NFFT=8000;
                Y=fft(rate,NFFT)/numel(rate);
                f=1000/2*linspace(0,1,NFFT/2+1);
                c=2*abs(Y(1:NFFT/2+1));
                c=c';
                
                tmp1=find(f>0.8);
                tmp1=tmp1(1);
                tmp2=find(f<32);
                tmp2=tmp2(end);
                tmp3=find(f<=7.5);
                tmp3=tmp3(end);
                %
                
                
                
                if plotting==1
                    figure(1)
                    subplot(2,6,ncnt+6)
                    plot(f(tmp1:tmp3),c(tmp1:tmp3),'color',co(sicnt))
                    hold on
                    subplot(2,6,sicnt)
                    plot(f(tmp1:tmp3),c(tmp1:tmp3),'color',co(ncnt))
                    hold on
                    
                    figure(2)
                    subplot(2,6,ncnt+6)
                    plot(fS(tmp1:tmp3),cS(tmp1:tmp3),'color',co(sicnt))
                    hold on
                    subplot(2,6,sicnt)
                    plot(fS(tmp1:tmp3),cS(tmp1:tmp3),'color',co(ncnt))
                    hold on
                end
                
                c=c(tmp1:tmp3); cS=cS(tmp1:tmp3);
                cc=c/max(c); ccS=cS/max(cS); ratio=cc./ccS;
                
                chirpFFTabs{i,sicnt+1,ncnt}=c;
                chirpFFTabs{i,7,ncnt}=cS;
                chirpFFTabs{i,1,ncnt}=f(tmp1:tmp3);
                
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
save(fullfile(superpath,'newAttempt2','chirp',[Dtype,'_chirpFFT']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates')

%% a) chirp stats (binary, NFFT=8000)
clear
close all

nfft=4; %2000,4000,6000,8000,10000,15000
nlist=[2000 4000 6000 8000 10000 15000];
superpath='f:\all_species_temp';

Dtype='wph';
cols='bgr';

colldata2=cell(3,1); colldata1a=cell(3,1); colldata1b=cell(3,1);colldata3=cell(3,1);
frs1=cell(3,1); frs2=cell(3,1); togos=cell(3,1);
figure(1)
for d=1:3
    load(fullfile(superpath,'newAttempt2','chirp',[Dtype(d),'_chirpFFT']))
    
    resp=cell2mat(reshape(chirpFFTabs(:,6,nfft),size(chirpFFTabs,1),1)')'; %non-smooth binary
   realgood=[];
for ts=1:size(chirpFFTabs,1)
    if ~isempty(chirpFFTabs{ts,1,1})
        realgood=[realgood;ts];
    end
end
togo=togo(realgood);
    stim=cell2mat(reshape(chirpFFTabs(:,7,nfft),size(chirpFFTabs,1),1)')';
    data2=cell2mat(reshape(chirpFFTsmooth(:,6,nfft),size(chirpFFTabs,1),1)')'; %smooth binary
    fr1=cell2mat(reshape(chirpFFTabs(:,1,nfft),size(chirpFFTabs,1),1));
    fr2=cell2mat(reshape(chirpFFTsmooth(:,1,nfft),size(chirpFFTabs,1),1)')';
    frs1{d}=fr1;
    frs2{d}=fr2;
    
    clear chirpFFTabs chirpFFTsmooth
    
    data1a=resp./stim;
    resp2=[]; stim2=[];
    for i=1:size(resp,1)
        resp2=[resp2;resp(i,:)/max(resp(i,:))];
        stim2=[stim2;stim(i,:)/max(stim(i,:))];
    end
    data1b=resp2./stim2;
    
    data3=[];
    for i=1:size(data2,1)
        data3=[data3;data2(i,:)/max(data2(i,:))];
    end
    
    subplot(4,4,d)
    plot(fr1(1,:),data1a,'color',[0.6 0.6 0.6])
    hold on
    subplot(4,4,d+4)
    plot(fr1(1,:),data1b,'color',[0.6 0.6 0.6])
    hold on
    subplot(4,4,d+8)
    plot(fr2(1,:),data2,'color',[0.6 0.6 0.6])
    hold on
    subplot(4,4,d+12)
    plot(fr2(1,:),data3,'color',[0.6 0.6 0.6])
    hold on
    
    subplot(4,4,d)
    plot(fr1,median(data1a),'color',cols(d),'linewidth',2)
    patch([fr1 fliplr(fr1)],[median(data1a)-std(data1a) fliplr(median(data1a)+std(data1a))],cols(d),'facealpha',0.4)
    %     plot(fr1,median(data1a)+std(data1a),'color',cols(d),'linewidth',1)
    %     plot(fr1,median(data1a)-std(data1a),'color',cols(d),'linewidth',1)
    title('non-smoothed FFTresp/FFTstim')
    axis tight
    subplot(4,4,4)
    plot(fr1,median(data1a),'color',cols(d),'linewidth',2)
    hold on
    patch([fr1 fliplr(fr1)],[median(data1a)-std(data1a) fliplr(median(data1a)+std(data1a))],cols(d),'facealpha',0.2)
    %     plot(fr1,median(data1a)+std(data1a),'color',cols(d),'linewidth',1)
    %     plot(fr1,median(data1a)-std(data1a),'color',cols(d),'linewidth',1)
    axis tight
    
    subplot(4,4,d+4)
    plot(fr1,median(data1b),'color',cols(d),'linewidth',2)
    patch([fr1 fliplr(fr1)],[median(data1b)-std(data1b) fliplr(median(data1b)+std(data1b))],cols(d),'facealpha',0.4)
    %     plot(fr1,median(data1b)+std(data1b),'color',cols(d),'linewidth',1)
    %     plot(fr1,median(data1b)-std(data1b),'color',cols(d),'linewidth',1)
    title('non-smoothed norm(FFTresp)/norm(FFTstim)')
    axis tight
    subplot(4,4,8)
    plot(fr1,median(data1b),'color',cols(d),'linewidth',2)
    hold on
    patch([fr1 fliplr(fr1)],[median(data1b)-std(data1b) fliplr(median(data1b)+std(data1b))],cols(d),'facealpha',0.2)
    %     plot(fr1,median(data1b)+std(data1b),'color',cols(d),'linewidth',1)
    %     plot(fr1,median(data1b)-std(data1b),'color',cols(d),'linewidth',1)
    axis tight
    
    subplot(4,4,d+8)
    plot(fr2,median(data2),'color',cols(d),'linewidth',2)
    patch([fr2 fliplr(fr2)],[median(data2)-std(data2) fliplr(median(data2)+std(data2))],cols(d),'facealpha',0.4)
    %     plot(fr2,median(data2)+std(data2),'color',cols(d),'linewidth',1)
    %     plot(fr2,median(data2)-std(data2),'color',cols(d),'linewidth',1)
    title('smoothed norm(FFTresp)/norm(FFTstim)')
    axis tight
    
    subplot(4,4,12)
    plot(fr2,median(data2),'color',cols(d),'linewidth',2)
    hold on
    patch([fr2 fliplr(fr2)],[median(data2)-std(data2) fliplr(median(data2)+std(data2))],cols(d),'facealpha',0.2)
    
    %     plot(fr2,median(data2)+std(data2),'color',cols(d),'linewidth',1)
    %     plot(fr2,median(data2)-std(data2),'color',cols(d),'linewidth',1)
    axis tight
    
    figure(2)
    
    plot(fr2,median(data2),'color',cols(d),'linewidth',2)
    hold on
    %     patch([fr2 fliplr(fr2)],[median(data2)-std(data2) fliplr(median(data2)+std(data2))],cols(d),'facealpha',0.2)
    axis tight
    
    figure(1)
    subplot(4,4,d+12)
    plot(fr2,median(data3),'color',cols(d),'linewidth',2)
    patch([fr2 fliplr(fr2)],[median(data3)-std(data3) fliplr(median(data3)+std(data3))],cols(d),'facealpha',0.4)
    %     plot(fr2,median(data2)+std(data2),'color',cols(d),'linewidth',1)
    %     plot(fr2,median(data2)-std(data2),'color',cols(d),'linewidth',1)
    title('norm(smoothed nnorm(FFTresp)/norm(FFTstim))')
    axis tight
    
    subplot(4,4,16)
    plot(fr2,median(data3),'color',cols(d),'linewidth',2)
    hold on
    patch([fr2 fliplr(fr2)],[median(data3)-std(data3) fliplr(median(data3)+std(data3))],cols(d),'facealpha',0.2)
    
    %     plot(fr2,median(data2)+std(data2),'color',cols(d),'linewidth',1)
    %     plot(fr2,median(data2)-std(data2),'color',cols(d),'linewidth',1)
    axis tight
    
    colldata2{d}=data2; %smooth resp/stim (mouse/pig/human)
    colldata3{d}=data3; %norm-smooth resp/stim (mouse/pig/human)
    colldata1a{d}=data1a; %non-smooth rsp/stim
    colldata1b{d}=data1b; %non-smooth norm(resp)/norm(stim)
    togos{d}=togo;
end
set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(superpath,'newAttempt2',['chirp_comparison',int2str(nlist(nfft)),'.bmp']))
% saveas(gcf,fullfile(superpath,'newAttempt2',['chirp_comparison',int2str(nlist(nfft)),'.fig']))
save(fullfile(superpath,'newAttempt2','chirp',['chirpRatios_',int2str(nlist(nfft))]),'colldata2','colldata3','colldata1a','colldata1b','frs1','frs2','togos')

figure(2)
title(['w: ',int2str(length(togos{1})),', p: ',int2str(length(togos{2})),' ,h: ',int2str(length(togos{3}))])
set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(superpath,'newAttempt2',['chirp_comparison',int2str(nlist(nfft)),'.eps']))

%% b) !!! continues -> plotting
clear
close all
nfft=4;
nlist=[2000 4000 6000 8000 10000 15000];
superpath='f:\all_species_temp';
load(fullfile(superpath,'newAttempt2','chirp',['chirpRatios_',int2str(nlist(nfft))]))

step=5;
data=colldata2;
fr=frs2{1};
dataset='smooth'; %data2
% dataset='noSmooth_noNorm'; %data1a
% dataset='noSmooth_Norm'; %data1b

pvals=zeros(ceil(size(data{1},2)/step),3); %M-P, M-H, P-H
cnt=0;
for t=1:step:size(data{1},2)
    cnt=cnt+1;
    now=t:t+step-1;
    if t+step-1 > size(data{1},2)
        now=t:size(data{1},2);
    end
    p=ranksum(mean(data{1}(:,now),2),mean(data{2}(:,now),2));
    pvals(cnt,1)=p;
    p=ranksum(mean(data{1}(:,now),2),mean(data{3}(:,now),2));
    pvals(cnt,2)=p;
    p=ranksum(mean(data{2}(:,now),2),mean(data{3}(:,now),2));
    pvals(cnt,3)=p;
end

figure
startpoint=floor(step/2);
if startpoint==0
    startpoint=1;
end
if length(fr(1,startpoint:step:end))<size(pvals,1)
    frnow=fr(1,startpoint:step:end);
    frnow=[frnow frnow(end)+diff([frnow(1) frnow(2)])];
else
    frnow=fr(1,startpoint:step:end);
end
subplot(2,2,1)
% patch([fr(1,startpoint:step:end)' fliplr(fr(1,startpoint:step:end)')],[pvals(:,1) repmat(0,length(fr(1,startpoint:step:end)'),1)],[0 0.5 0.5],'edgecolor',[0 0.5 0.5])
plot(frnow,pvals(:,1),'color',[0 0.5 0.5],'linewidth',3)
title('mouse - pig')
subplot(2,2,2)
% patch([fr(1,startpoint:step:end)' fliplr(fr(1,startpoint:step:end)')],[pvals(:,2) repmat(0,length(fr(1,startpoint:step:end)'),1)],[0.5 0 0.5],'edgecolor',[0.5 0 0.5])
plot(frnow,pvals(:,2),'color',[0.5 0 0.5],'linewidth',3)
title('mouse - human')
subplot(2,2,3)
% patch([fr(1,startpoint:step:end)' fliplr(fr(1,startpoint:step:end)')],[repmat(-0.001,length(fr(1,startpoint:step:end)'),1) fliplr(pvals(:,3)) ],[0.5 0.5 0],'edgecolor',[0.5 0.5 0])
plot(frnow,pvals(:,3),'color',[0.5 0.5 0],'linewidth',3)
title('pig - human')
for s=1:3
    subplot(2,2,s)
    hold on
%     plot([fr(1) fr(end)],[0.05 0.05],'-k')
    plot([fr(1) fr(end)],[0.01 0.01],'-k')
    axis([fr(1) fr(end) 0 0.05])
    xlabel('fr [Hz]')
    ylabel('p-value')
end
% saveas(gcf,fullfile(superpath,'newAttempt','stat_figs',['chirpPval_',dataset,'_',int2str(step),'step',int2str(nlist(nfft)),'.bmp']))
% saveas(gcf,fullfile(superpath,'newAttempt','stat_figs',['chirpPval_',dataset,'_',int2str(step),'step',int2str(nlist(nfft)),'.fig']))
% close
% save(fullfile(superpath,'newAttempt','chirp',[Dtype,'_chirpFFT']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates','varInfo')