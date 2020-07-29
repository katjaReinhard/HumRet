clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
% path1='R:\Users\Katja\human_retina';
path1='C:\Users\katja\OneDrive - imec\HumRet';
load(fullfile(path1,'h_responding'))

%% plot example cells: flash
close all



load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);
load(fullfile(path1,'h_flashSpikes'))



toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];

% cd('/media/NERFFS01/Users/Katja/human_retina/_scripts')
cd('C:\Users\katja\OneDrive - imec\HumRet\_scripts')

for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    
    data=flash_raster(ID2,:);
    
    figure
    allrate=[];
    for f=1:length(data)
        spikes=data{f};
        rate=MEA_spikerates(spikes,40,15000);
        rate=rate(200:10500);
        allrate=[allrate;rate];
        spikes=spikes(find(spikes>=200));
        spikes=spikes-200+1;
        spikes=spikes(find(spikes<=5800));
        subplot(2,1,1)
        for s=1:length(spikes)
            plot([spikes(s) spikes(s)],[f-0.4 f+0.4],'-b')
            hold on
            
        end
        %         plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
        %         hold on
    end
    plot([1800 1800],[0 f+1],'-k')
    plot([3800 3800],[0 f+1],'-k')
    axis([0 5800 0 f+1])
    title([toplotDate{g},'_',currUnit],'interpreter','none')
    subplot(2,1,2)
    meanrate=mean(allrate,1);
    plot(meanrate)
    axis([0 5800 0 inf])
    hold on
    ma=max(meanrate);
    plot([1800 1800],[0 ma],'-k')
    plot([3800 3800],[0 ma],'-k')
    
    set(gcf,'position',[1601 1 1600 824])
    name=['ex_flash_',int2str(resp(g)),'_',toplotDate{g},'_','currUnit.emf'];
%     saveas(gcf,fullfile(path1,'Figs',name))
    %     cd(fullfile(path1,'Figs'))
    %
    %     print(h,'-dpsc','test')
%     close
end

%% plot example cells: 6vel
close all



load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);
load(fullfile(path1,'h_vel6_new_Black'))

cd('R:\Users\Katja\human_retina\_scripts')

toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];

prot=read_header_field_heka('R:\Users\Katja\human_retina\Human\20120320c\HEKA','20120320_C1#0018_OPA1_bar_6veloc_downup_black%speed_200_gap_3000%_FW0ND4.phys','Stimulus Protocol');
cons=prot(1:end,2);
changes=find(cons==-1);
startpoints=round(prot(changes(1:7),1))-20500;
startpoints=startpoints(3:end-1);


for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    
    data=vel6_newB{ID2};
    
    figure
    dat=data(20500:end-2000);
    plot(dat)
    hold on
    
    for s=1:length(startpoints)
        plot([startpoints(s) startpoints(s)],[0 max(dat)],'-k')
    end
    
    axis([0 inf 0 inf])
    hold on
    
    set(gcf,'position',[1601 1 1600 824])
    name=['ex_bar_',int2str(resp(g)),'_',toplotDate{g},'_','currUnit.emf'];
    saveas(gcf,fullfile(path1,'Figs',name))
    %     cd(fullfile(path1,'Figs'))
    %
    %     print(h,'-dpsc','test')
    close
end



%% plot example cells: DS stimulus
close all



load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);
load(fullfile(path1,'h_DSBlack'))

% cd('R:\Users\Katja\human_retina\_scripts')
cd('C:\Users\KatjaR\OneDrive - imec\HumRet\_scripts')

toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];

% prot=read_header_field_heka('R:\Users\Katja\human_retina\Human\20120320c\HEKA','20120320_C1#0014_OPA1_DS_alldirections_black%speed_1000_gap_1000%_FW0ND4.phys','Stimulus Protocol');
prot=read_header_field_heka('R:\Users\Katja\human_retina\Human\20120320c\HEKA','20120320_C1#0014_OPA1_DS_alldirections_black%speed_1000_gap_1000%_FW0ND4.phys','Stimulus Protocol');

cons=prot(1:end,2);
changes=find(cons==-1);
startpoints=round(prot(changes(1:8),1))-500;


for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    
    data=DS_ratesB{ID2};
    
    figure
    dat=data(500:end-2000);
    plot(dat)
    hold on
    
    for s=1:length(startpoints)
        plot([startpoints(s) startpoints(s)],[0 max(dat)],'-k')
    end
    
    axis([0 inf 0 inf])
    hold on
    
    set(gcf,'position',[1601 1 1600 824])
    name=['ex_DSstim_',int2str(resp(g)),'_',toplotDate{g},'_','currUnit.emf'];
    saveas(gcf,fullfile(path1,'Figs',name))
    %     cd(fullfile(path1,'Figs'))
    %
    %     print(h,'-dpsc','test')
    close
end

%% plot example cells: best DG 
close all



load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);
load(fullfile(path1,'h_DGrates'))
load(fullfile(path1,'h_normPeaks'))

% cd('R:\Users\Katja\human_retina\_scripts')
cd('C:\Users\KatjaR\OneDrive - imec\HumRet\_scripts')

toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];



for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    ID3=find(togo==ID2);
    [dg id]=max(normPeaks(ID3,:,2));
    
    data=DGrates{ID2,id};
    
    figure
    dat=data(1000:14000);
    plot(dat)
    
    title(['DG number ',int2str(id)])
    axis([0 inf 0 inf])
    
%     set(gcf,'position',[1601 1 1600 824])
    name=['ex_DGbest_',int2str(resp(g)),'_',toplotDate{g},'_','currUnit.emf'];
    saveas(gcf,fullfile(path1,'Figs',name))
    %     cd(fullfile(path1,'Figs'))
    %
    %     print(h,'-dpsc','test')
    close
end
%% plot example cells: heatmap
close all



load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);
load(fullfile(path1,'h_DGrates'))
load(fullfile(path1,'h_normPeaks'))

% cd('R:\Users\Katja\human_retina\_scripts')
cd('C:\Users\KatjaR\OneDrive - imec\HumRet\_scripts')

toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];

selection = [];

for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    ID3=find(togo==ID2);
    now = normPeaks(ID3,:,2);
    now2=fliplr(reshape(now,4,6));
    
    
    figure
   imagesc(now2,[0 1])
   colormap gray
   set(gca,'plotboxaspectratio',[1 1/6*4 1],'box','off')
   axis off
  
    name=['final_ex_DGbest_',int2str(resp(g)),'_',toplotDate{g},'_','currUnit.emf'];
    saveas(gcf,fullfile(path1,'Figs',name))
    %     cd(fullfile(path1,'Figs'))
    %
    %     print(h,'-dpsc','test')
    close
end


%% plot example cells: STA
close all



load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);
load(fullfile(path1,'h_LFs'))

cd('R:\Users\Katja\human_retina\_scripts')

toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];




for g=1:length(toplotDate)
    
    currUnit=strcat('000',num2str(toplotUnit(g)));
    currUnit=currUnit(end-3:end);
    cDate = ismember(rDates,toplotDate{g});
    cUnit = ismember(rUnits(cDate),currUnit);
    ID = find(cUnit==1);
    tmp = find(cDate==1);
    ID = tmp(ID);
    ID2=resp(ID);
    data=LFs(:,ID2);
    data=data(1:245);
    dataS=[];
    for d=1:7:238
        dataS=[dataS;mean(data(d:d+6))];
    end
    
    figure
    plot(1:7:238,dataS)
    hold on
    
    
    axis([-inf inf -inf inf])
    
    
    set(gcf,'position',[1601 1 1600 824])
    name=['ex_STA_',int2str(resp(g)),'_',toplotDate{g},'_','currUnit.emf'];
    saveas(gcf,fullfile(path1,'Figs',name))
    %     cd(fullfile(path1,'Figs'))
    %
    %     print(h,'-dpsc','test')
    close
end

%% plot example cells: chirp

load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);

cd('R:\Users\Katja\human_retina\_scripts')

toplotDate={'20120224','20120224','20120224','20130124a','20130211c','20130718b','20131210a','20130510a','20130510a','20130510d','20130510d','20130510d','20130718b','20130718b','20130124a','20130124a','20130718c'};
toplotUnit=[35 65 24 38 2 70 31 28 100 13 54 5 6 3 41 53 24];


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
        prot=read_header_field_heka(fullfile(dpath,'HEKA'),unit{1,2}{flashfiles(1),1},'Stimulus Protocol');
        start=prot(2,1);
        start=start-200;
        start2=prot(485,1);
        
        
        if g==length(toplotDate)
            
            
            int=prot(2:end-1,4);
            int=[128; int];
            
            time=prot(2:end-2,1);
            time=time-start;
            tmp=find(time<=10000);
            time=time(tmp);
            timediff=diff(time);
            timediff=timediff(1:482);
            
            
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
            int2=[int2;repmat(128,10000-length(int2),1)];
            
            figure
            plot(int2)
            set(gcf,'position',[1601 1 1600 824])
            saveas(gcf,fullfile(path1,'Figs','stim_chirpFreq.emf'))
            close
            
            
            
            
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
            
            figure
            plot(int2)
            set(gcf,'position',[1601 1 1600 824])
            saveas(gcf,fullfile(path1,'Figs','stim_chirpAmp.emf'))
            close
        end
        
        
        
        figure
        allrate=[];
        for f=1:length(flashfiles)
            spikes=unit{1,2}{flashfiles(f),2};
            rate=MEA_spikerates(spikes,10,13000);
            rate=rate(start:10000+start);
            allrate=[allrate;rate];
            subplot(2,1,1)
            spikes=spikes(find(spikes>=start));
            spikes=spikes-start+1;
            spikes=spikes(find(spikes<=10000));
            for s=1:length(spikes)
                plot([spikes(s) spikes(s)],[f-0.4 f+0.4],'-b')
                hold on
                
            end
            %         plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
            %         hold on
        end
        %     plot([0 100],[f+1 f+1],'-r')
        axis([0 10000 0 f+1])
%         title(ulist(m).name,'interpreter','none')
        subplot(2,1,2)
        meanrate=mean(allrate,1);
        plot(meanrate)
        axis([0 10000 0 inf])
        
        set(gcf,'position',[1601 1 1600 824])
        name=['ex_chirpFreq_',int2str(resp(g)),'_',toplotDate{g},'_',currUnit,'.emf'];
        saveas(gcf,fullfile(path1,'Figs',name))
        close
        
        
        
        figure
        allrate=[];
        for f=1:length(flashfiles)
            spikes=unit{1,2}{flashfiles(f),2};
            rate=MEA_spikerates(spikes,10,23000);
            rate=rate(start2:start2+8500);
            allrate=[allrate;rate];
            subplot(2,1,1)
            spikes=spikes(find(spikes>=start2));
            spikes=spikes-start2+1;
                    spikes=spikes(find(spikes<=8500));
            for s=1:length(spikes)
                plot([spikes(s) spikes(s)],[f-0.4 f+0.4],'-b')
                hold on
                
            end
            %         plot(spikes,repmat(f,length(spikes),1),'o','markerfacecolor','b','markersize',1)
            %         hold on
        end
        %     plot([0 100],[f+1 f+1],'-r')
        axis([0 8500 0 f+1])
%         title(ulist(m).name,'interpreter','none')
        subplot(2,1,2)
        meanrate=mean(allrate,1);
        plot(meanrate)
        axis([0 8500 0 inf])
        
        set(gcf,'position',[1601 1 1600 824])
        name=['ex_chirpAmp_',int2str(resp(g)),'_',toplotDate{g},'_',currUnit,'.emf'];
        saveas(gcf,fullfile(path1,'Figs',name))
        close
    end
end