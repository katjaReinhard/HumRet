%% PLAN:
% 1) make matrix containing all cells with cleanliness (standardized as
% very good, usable, none), part of data to be used (start file, end file),
% ND used, polarity from Full-field (3), polarity from LF if possible (4)

% 2) check DS

% 3) check for correlation -> exclude cells

% 4) do speed (MEA_extractTemporal_new, rewritten) and spatiotemporal
% extraction (MEA_extractSpatiotemporal, rewritten, followed by
% MEA_rawSpatemp, rewritten)

% 3) analyze chirp
%% adjustable variables
clear
close all

Dtype='w';
% superpath='S:\data\Katja\allspecies';
% superpath='/gpfs01/muench/data/Katja/allspecies';
superpath='C:\Users\Katja\Desktop\all_species_temp';
%% data



switch Dtype
    
    
    case 'h'
        % Human data
        dates={'20120210a','20120224','20120320c','20120608b','20130124a','20130211a',...,
            '20130211c','20130211d','20130510a','20130510d','20130510e','20130718b','20130718c','20131210a','20131217f'};
        mainpath='S:\data\Katja\Human';
        starts=[255 245 8 91 1 97 192 97 1 1 1 173 1 1 88]; %index as in unit files
        stops=[374 365 277 280 191 190 449 285 246 164 246 429 172 86 260];
        NDs=[4 4 4 4 4 3 4 4 4 4 4 3 4 4 3];
    case 'p'
        % Pig data
        dates={'20121127','20121211a','20130613a','20130613b','20130613c','20140723d','20120213','20120319a','20120319b'};
        mainpath='S:\data\Katja\Pig';
        starts=[1 1 1 1 1 156 383 1 361];
        stops=[190 190 164 164 164 285 502 270 630];
        NDs=[4 4 4 4 4 4 3 3 3];
        
    case 'm'
        % Mouse data OPA
        dates={'20120126','20120127','20120221','20120222'};
        mainpath='S:\data\Katja\OPA';
        starts=[363 363 363 363];
        stops=[543 543 543 543];
        NDs=[4 4 4 4];
    case 'w'
        % Mouse data C3H
        dates={'20150114b','20150128b','20150203b','20150324b','20150424a','20150507a','20150507b'};
        %         mainpath='S:\data\Katja\WT';
        mainpath='C:\Users\Katja\Desktop\all_species_temp';
        starts=[412 412 412 412 412 412 412];
        stops=[577 577 577 577 577 577 577];
        NDs=[4 4 4 4 4 4 4];
end
%% 1) make info variable: first time, only basic info
% col1) date
% col2) unit number
% col3) quality: 2=great, 1=ok, 0=not usable
info=cell(10,3);
scnt=0; goodones=[];
for d=1:length(dates)
    currd=dates{d};
    upath=fullfile(mainpath,currd,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    for u=1:length(ulist)
        scnt=scnt+1;
        info{scnt,1}=currd;
        info{scnt,2}=ulist(u).name(end-7:end-4);
        load(fullfile(upath,ulist(u).name))
        if size(unit{1,1},1)>5
            cle=unit{1,1}{6,2};
            if ~isnumeric(cle)
                cle=str2num(cle);
            end
        else
            tmp=regexp(ulist(u).name,'_unit');
            tmp2=regexp(ulist(u).name,'-');
            tmp3=find(tmp2<tmp);
            tmp2=tmp2(tmp3);
            cle=str2num(ulist(u).name(tmp2+1:tmp-1));
        end
        
        if cle<10 %old
            if cle==1
                qua=2;
            elseif cle==-1 || cle==-2
                qua=1;
            else
                qua=0;
            end
        else %new
            if cle>85
                qua=2;
            elseif cle>70
                qua=1;
            else
                qua=0;
            end
        end
        info{scnt,3}=qua;
        %         disp([int2str(scnt),'/',int2str(qua)])
        if qua>0
            goodones=[goodones;1];
        else
            goodones=[goodones;0];
        end
    end
end

save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 2) add more info: start and stop file, ND
% col4) starts
% col5) stops
% col6) ND
load(fullfile(superpath,[Dtype,'_info']))

for i=1:length(info)
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(dates{m},info{i,1}))
            found=1;
        end
    end
    info{i,4}=starts(m);
    info{i,5}=stops(m);
    info{i,6}=NDs(m);
end

save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 3) polarity judgement from full-field
% col7) #quick files
% col8) pol from flash: 0=none, -1=off, 1=on, 2=onoff

load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
if ~exist(fullfile(superpath,'polarity_FFFlash'),'dir')
    mkdir(fullfile(superpath,'polarity_FFFlash'))
end
flash_rates=cell(length(info),1);

togo=find(goodones>0);
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
    
    quicks=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'FFFlash'))
            ftype=1;
            quicks=[quicks;1];
        elseif ~isempty(regexp(hekalist{f},'quick'))
            ftype=2;
            quicks=[quicks;1];
        else
            quicks=[quicks;0];
        end
    end
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(quicks==1);
    if Dtype=='m'
        goodfiles=goodfiles(11:end);
    end
    info{currID,7}=length(goodfiles);
    
    if ftype==1 %FFFlash
        first=goodfiles(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        cols=[0.6 0.3 0.6 0.9 0.6 0.6];
        if prot(6,1)<12000
            tt=12000;
            tim=[0; round(prot(2:6,1)); tt];
        else
            tt=14000;
            tim=[0; round(prot(2:6,1)); tt];
        end
        order=[1 4 2 3];
    else %quick
        first=goodfiles(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        timing=find(prot(:,4)==50);
        timing=round(prot(timing,1));
        timing=timing-2500;
        di=diff(timing);
        di=median(di);
        for t=2:length(timing)+1
            timing(t)=timing(t-1)+di;
        end
        cols=[0.6 0.9 0.6 0.3 0.6];
        tt=75000;
        tim=[0 500 2500 4500 9500 11500 13500];
        order=[2 3 1 4];
    end
    
    allrate=[];
    for f=goodfiles
        spikes=unit{1,2}{f,2};
        rate=MEA_spikerates(spikes,40,tt);
        allrate=[allrate;rate];
    end
    meanrate=mean(allrate);
    
    if ftype==1
        mrate=meanrate;
    else
        mrate=[];
        for t=1:length(timing)-1
            curr=meanrate(timing(t):timing(t+1));
            mrate=[mrate;curr];
        end
    end
    flash_rates{currID,1}=mrate;
    
    figure
    for t=1:length(tim)-1
        rectangle('position',[tim(t),0,tim(t+1)-tim(t),max(mrate)+2],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
        hold on
    end
    backgr=mean(mrate(1:500));
    stdBack=std(mrate(1:500));
    maxvals=[];
    for t=2:length(tim)-1
        curr=mrate(tim(t)+50:tim(t)+350);
        curr=curr-backgr;
        plot([tim(t)+350 tim(t)+350],[0 max(mrate)+2],'--b')
        maxval=max(curr);
        maxvals=[maxvals;maxval];
    end
    maxvals=maxvals(order);
    if ~isempty(regexp(Dtype,'h'))
        tmp=find(maxvals>3*stdBack);
    else
        tmp=find(maxvals>6*stdBack);
    end
    tmp2=find(maxvals(tmp)>=2);
    tmp=tmp(tmp2);
    off=0; on=0;
    if ~isempty(tmp)
        if ~isempty(find(tmp==1)) || ~isempty(find(tmp==2));
            off=1;
        end
        if ~isempty(find(tmp==3)) || ~isempty(find(tmp==4));
            on=1;
        end
        if off==1 && on==0
            pol=-1;
        elseif off==0 && on==1
            pol=1;
        else
            pol=2;
        end
    else
        pol=0;
    end
    if pol~=0
        plot(mrate,'-r','linewidth',2)
    else
        plot(mrate,'--k')
    end
    title([currd,'_',curru,' // pol: ',num2str(pol)],'interpreter','none')
    axis tight
    set(gcf,'position',[-46 387 867 420])
    choice = questdlg('Detected Polarity', ...
        'Automatic', ...
        'OK','Nope','OK');
    
    switch choice
        case 'Nope'
            
            prompt = {'correct polarity'};
            dlg_title = 'Polarity';
            num_lines = 1;
            def = {num2str(1-pol)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            oldpol=pol;
            pol=str2num(answer{1});
            
            if pol~=0
                plot(mrate,'-r','linewidth',2)
            else
                plot(mrate,'--k')
            end
            title([currd,'_',curru,' // pol: ',num2str(oldpol), ', corrected: ',num2str(pol)],'interpreter','none')
        case 'OK'
    end
    
    info{currID,8}=pol;
    saveas(gcf,fullfile(superpath,'polarity_FFFlash',[Dtype,'_',int2str(i),'_',currd,'_',curru,'.bmp']))
    close
    
end

save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
save(fullfile(superpath,[Dtype,'_flashrates.mat']),'flash_rates')
%% 3b) corrections human
togo=find(goodones>0);
ids=[241 650 836 841 168 577 762 767];
pols=[-1 2 -1 -1 0 0 -1 -1];
for i=1:length(ids)
    info{togo(ids(i)),8}=pols(i);
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 3c) corrections mouse
togo=find(goodones>0);
ids=[3 40 19];
pols=[1 1 1];
for i=1:length(ids)
    info{togo(ids(i)),8}=pols(i);
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 3d) corrections pig
togo=find(goodones>0);
ids=[148 206 285 317];
pols=[2 1 1 -1];
for i=1:length(ids)
    info{togo(ids(i)),8}=pols(i);
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 4) LF calculations and plotting
% col9) polarity from LF

apath=mainpath;
savepath=fullfile(superpath,'LFs');
if ~exist(savepath,'dir')
    mkdir(savepath)
end

for dd=1:length(dates)
    date=dates{dd};
    dpath=fullfile(apath,date);
    hekapath=fullfile(dpath,'HEKA');
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    load(fullfile(upath,ulist(1).name));
    heka=unit{1,2}(:,1);
    hekaNames=heka;
    
    found=0;mm=starts(dd);
    while found==0
        mm=mm+1;
        if ~isempty(regexp(hekaNames{mm},'FFFlicker'))
            found=1;
            ftype=1;
            codeWord='FFFlicker';
            maxr=256;
        elseif ~isempty(regexp(hekaNames{mm},'HCseq'))
            found=1;
            ftype=2;
            codeWord='HCseq';
            maxr=60;
        end
    end
    
    % select heka files
    hekas2consider=starts(dd):stops(dd);
    file_list=[];
    for i=hekas2consider
        if ~isempty(regexp(hekaNames{i},codeWord,'once'))
            file_list=[file_list i];
        end
    end
    
    % get protocol
    gotit=0;tester=2;
    while gotit==0
        try
            protocol=read_header_field_heka(hekapath, hekaNames{file_list(tester)}, 'Stimulus Protocol');
            gotit=1;
        catch
            tester=tester+1;
        end
    end
    temp=find(protocol(:,4)>1);
    protocol=protocol(temp,:);
    protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
    correctedProtocols=zeros(size(protocol,1),2,length(file_list));
    
    tmp=(0:size(correctedProtocols,1)-1)'*(7.85/3600);
    ndFullList=zeros(1,length(file_list));
    
    for i=1:length(file_list)
        hekaFile=hekaNames{file_list(i)};
        try
            protocol=read_header_field_heka(hekapath, hekaFile, 'Stimulus Protocol');
            tester=find(protocol(:,4)>1);
        catch
            tmp3=regexp(hekaFile,codeWord);
            tmp2=regexp(hekaFile(tmp3:end),'_');
            code2=hekaFile(tmp3:tmp2(1)+tmp3-2);
            
            found=0;m=0;
            while found==0
                m=m+1;
                if ~isempty(regexp(templist(m).name,code2))
                    found=1;
                end
            end
            protocol=read_header_field_heka('C:\Users\Katja\Desktop\Basel\protocol_temp', templist(m).name, 'Stimulus Protocol');
            
        end
        
        temp=find(protocol(:,4)>1);
        protocol=protocol(temp,:);
        protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
        protocol(:,1)=protocol(:,1)+tmp(1:size(protocol,1));
        if size(protocol,1)<size(correctedProtocols,1)
            protocol=[protocol;protocol(end,1)+16 -1 2 30 0];
        end
        correctedProtocols(1:length(protocol),:,i)=round(protocol(:,[1,4]));
        
        a=regexp(hekaFile,'ND');
        a=a+2;
        a=hekaFile(a);
        if length(a)==2
            ndFullList(i)=str2num(a(1))+str2num(a(2));
        else
            ndFullList(i)=str2num(a(1));
        end
    end
    clear protocol
    
    startTimes=reshape(correctedProtocols(1,1,:),1,length(file_list));
    endTimes=reshape(correctedProtocols(end,1,:),1,length(file_list));
    maxLength=max(endTimes-startTimes);
    flicker=zeros(maxLength,length(file_list));
    
    
    for i=1:length(file_list)
        
        flips=correctedProtocols(:,1,i);
        pxls=correctedProtocols(2:end,2,i);
        pxls=(pxls-(maxr/2))/(maxr/2)+maxr;
        
        flips=flips-flips(1)+1;
        
        if flips(end)<0
            flips=flips(1:end-1);
            pxls=pxls(1:end-1);
        end
        tmp=zeros(flips(end),1);
        tmp(flips(2:end),1)=pxls;
        
        % fill stimulus down to ms with brightness values where spikes
        % happen
        droppedMax=max(diff(flips));
        for j=1:droppedMax+1
            subFrame=flips(2:end)-j;
            subFrame(subFrame<1)=[];
            subFrame(find(tmp(subFrame,1)))=[];
            tmp(subFrame,1)=tmp(subFrame+1,1);
        end
        tmp=tmp-256;
        flicker(1:length(tmp),i)=tmp;
    end
    
    save(fullfile(savepath,[Dtype,'_protocols_',date]),'ndFullList','dpath','date','flicker','startTimes','maxLength','correctedProtocols','codeWord','file_list')
    
    
    
    per=1800; % frames to take for calculations
    filter_length=500;
    
    if exist('analyzed_cells','var')
        units=[]; IDs=[];
        for q=1:length(analyzed_cells)
            if ~isempty(regexp(analyzed_cells{q,1},date))
                units=[units;str2num(analyzed_cells{q,2})];
                IDs=[IDs;q];
            end
        end
    else
        units=[]; IDs=[];
        for q=1:length(ulist)
            units=[units;str2num(ulist(q).name(end-7:end-4))];
            IDs=[IDs;q];
        end
    end
    
    
    
    LinearFilter=zeros(filter_length,length(file_list),length(units));
    SpikeCount=zeros(length(file_list),length(units));
    if length(unique(startTimes))>1
        for i=1:size(correctedProtocols,3)
            correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
        end
    else
        correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
    end
    
    ulist=dir(fullfile(dpath,'units','*mat'));
    
    for cnt=1:length(units)
        unum=['000',int2str(units(cnt))];
        unum=unum(end-3:end);
        
        found=0; m=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(ulist(m).name,['unit_',unum]))
                found=1;
            end
        end
        load(fullfile(dpath,'units',ulist(m).name));
        
        
        %     colla=[]; colla2=[];
        for i=1:length(file_list)
            if file_list(i)<=size(unit{1,2},1)
                spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
                
                if ~isempty(spikes)
                    spikes=spikes-startTimes(i)+1;
                    startPoint=correctedProtocols(1,1,i)+1;
                    endPoint=correctedProtocols(per,1,i);
                    spikesTMP=spikes(spikes>startPoint+500&spikes<endPoint); %correct for filter length
                    
                    n=zeros(length(spikesTMP),filter_length);
                    %             for k=1:filter_length
                    for k=1:500
                        n(:,k)=flicker(spikesTMP-k+1,i); %flicker contains value for each ms
                        %                 n(:,k+500)=flicker(spikesTMP+k+1,i); %flicker contains value for each ms
                    end
                    LinearFilter(1:filter_length,i,cnt)=sum(n);
                    %             LinearFilter(1:filter_length,i,cnt)=mean(n);
                    SpikeCount(i,cnt)=length(spikesTMP);
                end
            end
        end
        %     figure
        %     subplot(2,1,1)
        %     hist(colla,50)
        %     subplot(2,1,2)
        %     hist(colla2,50)
    end
    
    save(fullfile(savepath,[Dtype,'_LF_',date]),'ndFullList','LinearFilter','SpikeCount','IDs','units','per')
    
end

load(fullfile(superpath,[Dtype,'_info']))
list=dir(fullfile(superpath,'LFs',[Dtype,'_LF_*']));
load(fullfile(superpath,'LFs',list(1).name))

LFs=zeros(500,length(units));


cnt=0;
for i=1:length(list)
    load(fullfile(superpath,'LFs',list(i).name))
    
    for j=1:length(IDs)
        cnt=cnt+1;
        names{cnt}=[info{cnt,1} '_', info{cnt,2}];
        LFs(:,cnt)=mean(LinearFilter(:,:,j),2);
    end
    
end


save(fullfile(superpath,'LFs',[Dtype,'_LFs']),'LFs','ndFullList','names')

togo=find(goodones>0);
for l=1:length(togo)
    currID=togo(l);
    curr=LFs(:,currID);
    figure
    rectangle('position',[50,min(curr)-2,200,max(curr)-min(curr)+4],'facecolor',[0.6 0.6 0.6])
    hold on
    plot(curr)
    
    axis tight
    set(gcf,'position',[-46 387 867 420])
    
    
    %     prompt = {'what is the polarity?'};
    %     dlg_title = 'LF';
    %     num_lines = 1;
    %     def = {'0'};
    %     answer = inputdlg(prompt,dlg_title,num_lines,def);
    %     info{currID,9}=str2num(answer{1});
    
    choice = questdlg('what is the polarity?', ...
        'LF', ...
        'ON','OFF','none','ON');
    switch choice
        case 'ON'
            info{currID,9}=1;
        case 'OFF'
            info{currID,9}=-1;
        case 'none'
            info{currID,9}=0;
    end
    
    title([names{currID},' // pol: ',num2str(info{currID,9})],'interpreter','none')
    saveas(gcf,fullfile(superpath,'LFs',[Dtype,'_',int2str(currID),'_',names{currID},'.bmp']))
    close
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 4b) corrections human
togo=find(goodones>0);
ids=[206 971 1061];
pols=[1 -1 -1];
for i=1:length(ids)
    info{ids(i),9}=pols(i);
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 4c) data loss!
savepath=fullfile(superpath,'LFs');
apath=mainpath; allunits=0; allIDs=0;
for dd=1:length(dates)
    dd
    date=dates{dd};
    dpath=fullfile(apath,date);
    hekapath=fullfile(dpath,'HEKA');
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    load(fullfile(upath,ulist(1).name));
    heka=unit{1,2}(:,1);
    hekaNames=heka;
    
    found=0;mm=starts(dd);
    while found==0
        mm=mm+1;
        if ~isempty(regexp(hekaNames{mm},'FFFlicker'))
            found=1;
            ftype=1;
            codeWord='FFFlicker';
            maxr=256;
        elseif ~isempty(regexp(hekaNames{mm},'HCseq'))
            found=1;
            ftype=2;
            codeWord='HCseq';
            maxr=60;
        end
    end
    
    % select heka files
    hekas2consider=starts(dd):stops(dd);
    file_list=[];
    for i=hekas2consider
        if ~isempty(regexp(hekaNames{i},codeWord,'once'))
            file_list=[file_list i];
        end
    end
    
    % get protocol
    gotit=0;tester=2;
    while gotit==0
        try
            protocol=read_header_field_heka(hekapath, hekaNames{file_list(tester)}, 'Stimulus Protocol');
            gotit=1;
        catch
            tester=tester+1;
        end
    end
    temp=find(protocol(:,4)>1);
    protocol=protocol(temp,:);
    protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
    correctedProtocols=zeros(size(protocol,1),2,length(file_list));
    
    tmp=(0:size(correctedProtocols,1)-1)'*(7.85/3600);
    ndFullList=zeros(1,length(file_list));
    
    for i=1:length(file_list)
        hekaFile=hekaNames{file_list(i)};
        try
            protocol=read_header_field_heka(hekapath, hekaFile, 'Stimulus Protocol');
            tester=find(protocol(:,4)>1);
        catch
            tmp3=regexp(hekaFile,codeWord);
            tmp2=regexp(hekaFile(tmp3:end),'_');
            code2=hekaFile(tmp3:tmp2(1)+tmp3-2);
            
            found=0;m=0;
            while found==0
                m=m+1;
                if ~isempty(regexp(templist(m).name,code2))
                    found=1;
                end
            end
            protocol=read_header_field_heka('C:\Users\Katja\Desktop\Basel\protocol_temp', templist(m).name, 'Stimulus Protocol');
            
        end
        
        temp=find(protocol(:,4)>1);
        protocol=protocol(temp,:);
        protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
        protocol(:,1)=protocol(:,1)+tmp(1:size(protocol,1));
        if size(protocol,1)<size(correctedProtocols,1)
            protocol=[protocol;protocol(end,1)+16 -1 2 30 0];
        end
        correctedProtocols(1:length(protocol),:,i)=round(protocol(:,[1,4]));
        
        a=regexp(hekaFile,'ND');
        a=a+2;
        a=hekaFile(a);
        if length(a)==2
            ndFullList(i)=str2num(a(1))+str2num(a(2));
        else
            ndFullList(i)=str2num(a(1));
        end
    end
    clear protocol
    
    startTimes=reshape(correctedProtocols(1,1,:),1,length(file_list));
    endTimes=reshape(correctedProtocols(end,1,:),1,length(file_list));
    maxLength=max(endTimes-startTimes);
    flicker=zeros(maxLength,length(file_list));
    
    
    for i=1:length(file_list)
        
        flips=correctedProtocols(:,1,i);
        pxls=correctedProtocols(2:end,2,i);
        pxls=(pxls-(maxr/2))/(maxr/2)+maxr;
        
        flips=flips-flips(1)+1;
        
        if flips(end)<0
            flips=flips(1:end-1);
            pxls=pxls(1:end-1);
        end
        tmp=zeros(flips(end),1);
        tmp(flips(2:end),1)=pxls;
        
        % fill stimulus down to ms with brightness values where spikes
        % happen
        droppedMax=max(diff(flips));
        for j=1:droppedMax+1
            subFrame=flips(2:end)-j;
            subFrame(subFrame<1)=[];
            subFrame(find(tmp(subFrame,1)))=[];
            tmp(subFrame,1)=tmp(subFrame+1,1);
        end
        tmp=tmp-256;
        flicker(1:length(tmp),i)=tmp;
    end
    
    save(fullfile(savepath,[Dtype,'_protocols_',date]),'ndFullList','dpath','date','flicker','startTimes','maxLength','correctedProtocols','codeWord','file_list')
    
    
    
    per=1800; % frames to take for calculations
    filter_length=500;
    
    if exist('analyzed_cells','var')
        units=[]; IDs=[];
        for q=1:length(analyzed_cells)
            if ~isempty(regexp(analyzed_cells{q,1},date))
                units=[units;str2num(analyzed_cells{q,2})];
                IDs=[IDs;q];
            end
        end
    else
        units=[]; IDs=[];
        for q=1:length(ulist)
            units=[units;str2num(ulist(q).name(end-7:end-4))];
            IDs=[IDs;q];
        end
    end
    
    allunits=allunits+length(units);
    allIDs=allIDs+length(IDs);
    
    LinearFilter=zeros(filter_length,length(file_list),length(units));
    SpikeCount=zeros(length(file_list),length(units));
    if length(unique(startTimes))>1
        for i=1:size(correctedProtocols,3)
            correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
        end
    else
        correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
    end
    
    ulist=dir(fullfile(dpath,'units','*mat'));
    
    for cnt=1:length(units)
        unum=['000',int2str(units(cnt))];
        unum=unum(end-3:end);
        
        found=0; m=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(ulist(m).name,['unit_',unum]))
                found=1;
            end
        end
        load(fullfile(dpath,'units',ulist(m).name));
        
        
        %     colla=[]; colla2=[];
        for i=1:length(file_list)
            if file_list(i)<=size(unit{1,2},1)
                spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
                
                if ~isempty(spikes)
                    spikes=spikes-startTimes(i)+1;
                    startPoint=correctedProtocols(1,1,i)+1;
                    endPoint=correctedProtocols(per,1,i);
                    spikesTMP=spikes(spikes>startPoint+500&spikes<endPoint); %correct for filter length
                    
                    n=zeros(length(spikesTMP),filter_length);
                    %             for k=1:filter_length
                    for k=1:500
                        n(:,k)=flicker(spikesTMP-k+1,i); %flicker contains value for each ms
                        %                 n(:,k+500)=flicker(spikesTMP+k+1,i); %flicker contains value for each ms
                    end
                    LinearFilter(1:filter_length,i,cnt)=sum(n);
                    %             LinearFilter(1:filter_length,i,cnt)=mean(n);
                    SpikeCount(i,cnt)=length(spikesTMP);
                end
            end
        end
        %     figure
        %     subplot(2,1,1)
        %     hist(colla,50)
        %     subplot(2,1,2)
        %     hist(colla2,50)
    end
    
    save(fullfile(savepath,[Dtype,'_LF_',date]),'ndFullList','LinearFilter','SpikeCount','IDs','units','per')
    
end

load(fullfile(superpath,[Dtype,'_info']))
list=dir(fullfile(superpath,'LFs',[Dtype,'_LF_*']));
load(fullfile(superpath,'LFs',list(1).name))

LFs=zeros(500,length(units));


cnt=0;
for i=1:length(list)
    load(fullfile(superpath,'LFs',list(i).name))
    
    for j=1:length(IDs)
        cnt=cnt+1;
        names{cnt}=[info{cnt,1} '_', info{cnt,2}];
        if info{cnt,9}~=0
            LFs(:,cnt)=mean(LinearFilter(:,:,j),2);
        else
            LFs(:,cnt)=zeros(500,1);
        end
    end
    
end


save(fullfile(superpath,'LFs',[Dtype,'_LFs']),'LFs','ndFullList','names')

% togo=find(goodones>0);
% for l=1:length(togo)
%     currID=togo(l);
%     curr=LFs(:,currID);
%     figure
%     rectangle('position',[50,min(curr)-2,200,max(curr)-min(curr)+4],'facecolor',[0.6 0.6 0.6])
%     hold on
%     plot(curr)
%
%     axis tight
%     set(gcf,'position',[-46 387 867 420])
%
%
%
%     title([names{currID},' // pol: ',num2str(info{currID,9})],'interpreter','none')
%     saveas(gcf,fullfile(superpath,'LFs',[Dtype,'_',int2str(currID),'_',names{currID},'.bmp']))
%     close
% end
%% -> 3) plotting afterwards due to data loss => do later for human and
load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
togo=find(goodones>0);
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
    
    quicks=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'FFFlash'))
            ftype=1;
            quicks=[quicks;1];
        elseif ~isempty(regexp(hekalist{f},'quick'))
            ftype=2;
            quicks=[quicks;1];
        else
            quicks=[quicks;0];
        end
    end
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(quicks==1);
    if Dtype=='m'
        goodfiles=goodfiles(11:end);
    end
    
    
    if ftype==1 %FFFlash
        first=goodfiles(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        cols=[0.6 0.3 0.6 0.9 0.6 0.6];
        if prot(6,1)<12000
            tt=12000;
            tim=[0; round(prot(2:6,1)); tt];
        else
            tt=14000;
            tim=[0; round(prot(2:6,1)); tt];
        end
        order=[1 4 2 3];
    else %quick
        first=goodfiles(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        timing=find(prot(:,4)==50);
        timing=round(prot(timing,1));
        timing=timing-2500;
        di=diff(timing);
        di=median(di);
        for t=2:length(timing)+1
            timing(t)=timing(t-1)+di;
        end
        cols=[0.6 0.9 0.6 0.3 0.6];
        tt=75000;
        tim=[0 500 2500 4500 9500 11500 13500];
        order=[2 3 1 4];
    end
    
    allrate=[];
    for f=goodfiles
        spikes=unit{1,2}{f,2};
        rate=MEA_spikerates(spikes,40,tt);
        allrate=[allrate;rate];
    end
    meanrate=mean(allrate);
    
    if ftype==1
        mrate=meanrate;
    else
        mrate=[];
        for t=1:length(timing)-1
            curr=meanrate(timing(t):timing(t+1));
            mrate=[mrate;curr];
        end
    end
    
    figure
    for t=1:length(tim)-1
        rectangle('position',[tim(t),0,tim(t+1)-tim(t),max(mrate)+2],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
        hold on
    end
    
    pol=info{currID,8};
    if pol~=0
        plot(mrate,'-r','linewidth',2)
    else
        plot(mrate,'--k')
    end
    title([currd,'_',curru,' // pol: ',num2str(pol)],'interpreter','none')
    axis tight
    set(gcf,'position',[-46 387 867 420])
    
    saveas(gcf,fullfile(superpath,'polarity_FFFlash',[Dtype,'_',int2str(i),'_',currd,'_',curru,'.bmp']))
    close
    
end
%% 5) check discrepancies between LF and full-field flash
load(fullfile(superpath,[Dtype,'_info']))

for i=1:size(info,1)
    currF=info{i,8};
    currL=info{i,9};
    currD=info{i,1};
    currU=info{i,2};
    
    gofor=0;
    
    if ~isempty(currF) || ~isempty(currL)
        if currF==1 && currL==-1
            gofor=1;
        end
        if currF==-1 && currL==1
            gofor=1;
        end
        if currF==2 && currL==1
            gofor=1;
        end
        
        if gofor==1
            winopen(fullfile(superpath,'LFs',[Dtype,'_',int2str(i),'_',currD,'_',currU,'.bmp']))
            
            choice = questdlg('want to correct LF?', ...
                'LF', ...
                'ON','OFF','none','ON');
            switch choice
                case 'ON'
                    info{i,9}=1;
                case 'OFF'
                    info{i,9}=-1;
                case 'none'
                    info{i,9}=0;
            end
        end
    end
end

save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 6) compute DS
% col10) DS index black
% col11) DS index white
load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
if ~exist(fullfile(superpath,'DS'),'dir')
    mkdir(fullfile(superpath,'DS'))
end
DS_ratesB=cell(length(info),1);
DS_ratesW=cell(length(info),1);

togo=find(goodones>0);
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
    
    black=[]; white=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DS_alldirections_black'))
            ftype=1;
            black=[black;1];
        else
            black=[black;0];
        end
        if ~isempty(regexp(hekalist{f},'DS_alldirections_white'))
            ftype=2;
            white=[white;1];
        else
            white=[white;0];
        end
    end
    goodfiles=currstart:currstop;
    goodfilesB=goodfiles(black==1);
    goodfilesW=goodfiles(white==1);
    
    first=goodfilesB(1);
    prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
    cons=prot(1:end,2);
    changes=find(cons==-1);
    startpoints=round(prot(changes(1:8),1));
    endpoints=round(prot(changes(2:9)-1,1));
    cols=[0.1 0.6 0.9 0.6 0.9 0.6 0.9 0.6 0.9 0.1];
    changes=[0; startpoints; endpoints(end); round(prot(end,1)+1000)];
    tt=round(prot(end,1)+1000);
    
    allrateB=[];
    for f=goodfilesB
        spikes=unit{1,2}{f,2};
        rate=MEA_spikerates(spikes,40,tt);
        allrateB=[allrateB;rate];
    end
    meanrateB=mean(allrateB);
    DS_ratesB{currID,1}=meanrateB;
    
    allrateW=[];
    for f=goodfilesW
        spikes=unit{1,2}{f,2};
        rate=MEA_spikerates(spikes,40,tt);
        allrateW=[allrateW;rate];
    end
    meanrateW=mean(allrateW);
    DS_ratesW{currID,1}=meanrateW;
    
    figure
    set(gcf,'position',[1601 1 1600 824])
    subplot(2,2,1)
    for c=1:length(changes)-1
        rectangle('position',[changes(c),0,changes(c+1)-changes(c),max(meanrateB)+2],'facecolor',[cols(c) cols(c) cols(c)],'edgecolor',[cols(c) cols(c) cols(c)])
        hold on
    end
    plot(meanrateB)
    axis tight
    
    subplot(2,2,2)
    for c=1:length(changes)-1
        rectangle('position',[changes(c),0,changes(c+1)-changes(c),max(meanrateW)+2],'facecolor',[cols(c) cols(c) cols(c)],'edgecolor',[cols(c) cols(c) cols(c)])
        hold on
    end
    plot(meanrateW)
    axis tight
    
    backB=mean(meanrateB(1:startpoints(1)));
    newrateB=meanrateB-backB;
    backW=mean(meanrateW(1:startpoints(1)));
    newrateW=meanrateW-backW;
    
    startpoints=[startpoints; endpoints(end)];
    vals=[];
    for fli=1:length(startpoints)-1
        curr=max(newrateB(startpoints(fli):startpoints(fli+1)));
        vals=[vals;curr];
    end
    normv=vals/max(vals);
    xcoord=(normv(1)-normv(2)+normv(8)/sqrt(2)+normv(5)/sqrt(2)-normv(6)/sqrt(2)-normv(7)/sqrt(2))/(sum(normv));
    ycoord=(normv(4)-normv(3)+normv(6)/sqrt(2)+normv(8)/sqrt(2)-normv(5)/sqrt(2)-normv(7)/sqrt(2))/(sum(normv));
    index=sqrt(xcoord^2+ycoord^2);
    info{currID,10}=index;
    %
    subplot(2,2,3)
    plot([0 normv(1)],[0 0],'-ok')
    hold on
    plot([0 -normv(2)],[0 0],'-ok')
    plot([0 0],[0 -normv(3)],'-ok')
    plot([0 0],[0 normv(4)],'-ok')
    plot([0 normv(5)/sqrt(2)],[0 -normv(5)/sqrt(2)],'-ok')
    plot([0 -normv(6)/sqrt(2)],[0 normv(5)/sqrt(6)],'-ok')
    plot([0 -normv(7)/sqrt(2)],[0 -normv(7)/sqrt(2)],'-ok')
    plot([0 normv(8)/sqrt(2)],[0 normv(8)/sqrt(2)],'-ok')
    plot([0 xcoord],[0 ycoord],'-r','linewidth',3)
    axis([-1 1 -1 1])
    rectangle('position',[-0.25,-0.25,0.5,0.5],'curvature',[1,1],'edgecolor','g')
    rectangle('position',[-0.35,-0.35,0.7,0.7],'curvature',[1,1],'edgecolor','g')
    rectangle('position',[-0.5,-0.5,1,1],'curvature',[1,1],'edgecolor','g')
    axis([-1 1 -1 1])
    
    
    vals=[];
    for fli=1:length(startpoints)-1
        curr=max(newrateW(startpoints(fli):startpoints(fli+1)));
        vals=[vals;curr];
    end
    normv=vals/max(vals);
    xcoord=(normv(1)-normv(2)+normv(8)/sqrt(2)+normv(5)/sqrt(2)-normv(6)/sqrt(2)-normv(7)/sqrt(2))/(sum(normv));
    ycoord=(normv(4)-normv(3)+normv(6)/sqrt(2)+normv(8)/sqrt(2)-normv(5)/sqrt(2)-normv(7)/sqrt(2))/(sum(normv));
    index=sqrt(xcoord^2+ycoord^2);
    info{currID,11}=index;
    
    subplot(2,2,4)
    plot([0 normv(1)],[0 0],'-ok')
    hold on
    plot([0 -normv(2)],[0 0],'-ok')
    plot([0 0],[0 -normv(3)],'-ok')
    plot([0 0],[0 normv(4)],'-ok')
    plot([0 normv(5)/sqrt(2)],[0 -normv(5)/sqrt(2)],'-ok')
    plot([0 -normv(6)/sqrt(2)],[0 normv(5)/sqrt(6)],'-ok')
    plot([0 -normv(7)/sqrt(2)],[0 -normv(7)/sqrt(2)],'-ok')
    plot([0 normv(8)/sqrt(2)],[0 normv(8)/sqrt(2)],'-ok')
    plot([0 xcoord],[0 ycoord],'-r','linewidth',3)
    axis([-1 1 -1 1])
    rectangle('position',[-0.25,-0.25,0.5,0.5],'curvature',[1,1],'edgecolor','g')
    rectangle('position',[-0.35,-0.35,0.7,0.7],'curvature',[1,1],'edgecolor','g')
    rectangle('position',[-0.5,-0.5,1,1],'curvature',[1,1],'edgecolor','g')
    axis([-1 1 -1 1])
    
    saveas(gcf,fullfile(superpath,'DS',[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
    close
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
save(fullfile(superpath,[Dtype,'_DSblack.mat']),'DS_ratesB')
save(fullfile(superpath,[Dtype,'_DSwhite.mat']),'DS_ratesW')
%% 6b) check DS
% col12) DS yes/no

load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)

togo=find(goodones>0);
for i=1:length(togo)
    curID=togo(i);
    currB=info{curID,10};
    currW=info{curID,11};
    currd=info{curID,1};
    curru=info{curID,2};
    
    if currB>=0.25 || currW>=0.25
        winopen(fullfile(superpath,'DS',[Dtype,'_',int2str(curID),'_',currd,'_',curru,'.bmp']))
        
        choice = questdlg('is it DS?', ...
            'DS', ...
            'yes','no','yes');
        switch choice
            case 'yes'
                info{curID,12}=1;
            case 'no'
                info{curID,12}=0;
        end
    else
        info{i,12}=0;
    end
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 7) get spont activity
load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)

spontRates=cell(3,1);
togo=find(goodones>0);
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
    
    spont=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'spont'))
            ftype=1;
            spont=[spont;1];
        else
            spont=[spont;0];
        end
        
    end
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(spont==1);
    
    allrate=[];
    for f=goodfiles
        spikes=unit{1,2}{f,2};
        rate=MEA_spikerates(spikes,40,60000);
        allrate=[allrate;rate];
    end
    if size(allrate,1)>1
        meanrate=mean(allrate);
    else
        meanrate=allrate;
    end
    spontRates{currID,1}=meanrate;
end
save(fullfile(superpath,[Dtype,'_spont']),'spontRates')
%% 8) check for correlation

load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_DSblack']))
load(fullfile(superpath,[Dtype,'_DSwhite']))
load(fullfile(superpath,[Dtype,'_flashrates']))
load(fullfile(superpath,[Dtype,'_spont']))
load(fullfile(superpath,'LFs',[Dtype,'_LFs']))

cd(superpath)
if ~exist(fullfile(superpath,'correlation'),'dir')
    mkdir(fullfile(superpath,'correlation'))
end

togo=find(goodones>0);
corrs=zeros(size(info,1),5,2);
cellIDs=zeros(size(info,1),2);
for i=1:length(togo)
    currID=togo(i);
    
    currd=info{currID,1};
    curru=info{currID,2};
    
    thisdate=[];
    for f=1:length(info)
        if ~isempty(regexp(info{f,1},currd)) && info{f,3}>0
            thisdate=[thisdate;f];
        end
    end
    
    thisdate=setdiff(thisdate,currID);
    
    
    fl1=flash_rates{currID}(500:end-1000);
    lf1=LFs(1:300,currID);
    dsb1=DS_ratesB{currID}(500:end-1000);
    dsw1=DS_ratesW{currID}(500:end-1000);
    sp1=spontRates{currID}(500:end-1000);
    
    for c=1:length(thisdate)
        currc=thisdate(c);
        cellIDs(currID,c)=currc;
        %flash
        fl2=flash_rates{currc}(500:end-1000);
        if ~isempty([find(fl1>0) find(fl2>0)])
            [rho1,pval]=corr(fl1',fl2');
        else
            rho1=0;
        end
        
        %LF
        lf2=LFs(1:300,currc);
        if ~isempty([find(lf1'>0) find(lf2'>0)])
            [rho2,pval]=corr(lf1,lf2);
        else
            rho2=0;
        end
        
        %DSblack
        dsb2=DS_ratesB{currc}(500:end-1000);
        if ~isempty([find(dsb1>0) find(dsb2>0)])
            [rho3,pval]=corr(dsb1',dsb2');
        else
            rho3=0;
        end
        
        %DSwhite
        dsw2=DS_ratesW{currc}(500:end-1000);
        if ~isempty([find(dsw1>0) find(dsw2>0)])
            [rho4,pval]=corr(dsw1',dsw2');
        else
            rho4=0;
        end
        
        %spont
        sp2=spontRates{currc}(500:end-1000);
        if ~isempty([find(sp1>0) find(sp1>0)])
            [rho5,pval]=corr(sp1',sp2');
        else
            rho5=0;
        end
        
        corrs(currID,1,c)=rho1;
        corrs(currID,2,c)=rho2;
        corrs(currID,3,c)=rho3;
        corrs(currID,4,c)=rho4;
        corrs(currID,5,c)=rho5;
    end
end

save(fullfile(superpath,'correlation',[Dtype,'_corrs']),'corrs','cellIDs')
%% 8b) search for high correlations
clc
load(fullfile(superpath,'correlation',[Dtype,'_corrs']))
load(fullfile(superpath,[Dtype,'_info']))
% save(fullfile(superpath,[Dtype,'_info_keep']),'info','goodones') %safety copy

medcorrs=[]; stcorrs=[];
for i=1:size(corrs,2)
    curr=abs(reshape(corrs(:,i,:),size(corrs,1),size(corrs,3)));
    tmp=find(curr>0);
    curr=curr(tmp);
    tmp=find(~isnan(curr));
    curr=curr(tmp);
    m=median(curr);
    st=std(curr);
    %     figure
    %     hist(curr)
    
    medcorrs=[medcorrs m];
    stcorrs=[stcorrs st];
end

tocheck=cell(size(corrs,1),1);
cvals=cell(size(corrs,1),1);
for i=1:size(corrs,1)
    currdata=reshape(corrs(i,:,:),5,size(corrs,3));
    
    %         clear test
    %         for j=1:5
    %             test{j}=find(abs(currdata(j,:))>medcorrs(j)+2*stcorrs(j));
    %         end
    %
    %         candidates=intersect(test{1},test{2});
    %         candidates=intersect(candidates,test{3});
    %         candidates=intersect(candidates,test{4});
    %         candidates=intersect(candidates,test{5});
    
    
    candidates=find(currdata(5,:)>0.6);
    newcands=candidates;
    for cc=1:length(candidates)
        cnt=0;
        for jj=[1 3 4]
            if ~isempty(find(currdata(jj,candidates(cc))>medcorrs(jj)+2*stcorrs(jj)))
                cnt=cnt+1;
            end
        end
        if cnt<3 || isempty(find(currdata(2,candidates(cc))>0.5))
            newcands=setdiff(newcands,candidates(cc));
        elseif info{cellIDs(i,candidates(cc)),9}==-1 && info{i,9}==1
            newcands=setdiff(newcands,candidates(cc));
        elseif info{cellIDs(i,candidates(cc)),9}==1 && info{i,9}==-1
            newcands=setdiff(newcands,candidates(cc));
        end
    end
    candidates=newcands;
    tocheck{i}=cellIDs(i,candidates);
    for cc=1:length(candidates)
        cvals{i,cc}=currdata(:,candidates(cc));
    end
end

for t=1:length(tocheck)
    if ~isempty(tocheck{t})
        t
        cands=tocheck{t};
        allcands=[t cands];
        a=0;
        while a<length(allcands)
            a=a+1;
            qual=info{allcands(a),3};
            if qual<0
                allcands=setdiff(allcands,allcands(a));
            end
        end
        keepcands=allcands;
        while length(allcands)>1
            checkpoints=[0 0];
            %check LF
            first=info{allcands(1),9};
            second=info{allcands(2),9};
            if first~=0
                checkpoints(1,1)=1;
            end
            if second~=0
                checkpoints(1,2)=1;
            end
            
            %check flash
            first=info{allcands(1),8};
            second=info{allcands(2),8};
            if first~=0
                checkpoints(1,1)= checkpoints(1,1)+1;
            end
            if second~=0
                checkpoints(1,2)=checkpoints(1,2)+1;
            end
            if checkpoints(1)>checkpoints(2)
                allcands=setdiff(allcands,allcands(2));
            elseif checkpoints(1)<checkpoints(2)
                allcands=setdiff(allcands,allcands(1));
            else
                allcands=setdiff(allcands,allcands(2));
            end
        end
        lost=setdiff(keepcands,allcands);
        if ~isempty(lost)
            lost
            for l=1:length(lost)
                info{lost(l),3}=-1; %-1 = out because "copy" of other unit
            end
        end
        
    end
end
goodones=[];
for i=1:size(info,1)
    if info{i,3}>0
        goodones=[goodones;1];
    else
        goodones=[goodones;0];
    end
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 9) compute speed
% col13) #of files for 6vel
% col14) #maxid for 6vel black (speed for which it's best)
% col15) #maxid for 6vel white (optimal speed)


load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
if ~exist(fullfile(superpath,'6vel'),'dir')
    mkdir(fullfile(superpath,'6vel'))
end
names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[0.2 1 2 4 8 16];

vel6_ratesB=cell(length(info),1);
vel6_ratesW=cell(length(info),1);
togo=find(goodones>0);
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
    
    black=[]; white=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'bar_6veloc_downup_black'))
            ftype=1;
            black=[black;1];
        else
            black=[black;0];
        end
        if ~isempty(regexp(hekalist{f},'bar_6veloc_downup_white'))
            ftype=2;
            white=[white;1];
        else
            white=[white;0];
        end
    end
    goodfiles=currstart:currstop;
    goodfilesB=goodfiles(black==1);
    goodfilesW=goodfiles(white==1);
    
    if ~isempty(goodfilesB)
        info{currID,13}=length(goodfilesB);
        first=goodfilesB(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        cons=prot(1:end,2);
        changes=find(cons==-1);
        startpoints=round(prot(changes(1:7),1));
        tt=round(prot(end,1)+1000);
        
        allrateB=[];
        for f=goodfilesB
            spikes=unit{1,2}{f,2};
            rate=MEA_spikerates(spikes,40,tt);
            allrateB=[allrateB;rate];
        end
        meanrateB=mean(allrateB);
        vel6_ratesB{currID,1}=meanrateB;
        
        allrateW=[];
        for f=goodfilesW
            spikes=unit{1,2}{f,2};
            rate=MEA_spikerates(spikes,40,tt);
            allrateW=[allrateW;rate];
        end
        meanrateW=mean(allrateW);
        vel6_ratesW{currID,1}=meanrateW;
        
        
        clear data3
        avRate=meanrateB;
        avbackground=mean(avRate(1:startpoints(1)));
        for j=1:6
            [data3(j) data4(j)]=max(avRate(startpoints(j):startpoints(j+1)));
        end
        [maxval maxid]=max(data3);
        data3=data3-avbackground;
        
        
        total=sum(data3); sofar=0; nn=0;
        if total>0
            while sofar<total/2
                nn=nn+1;
                sofar=sofar+data3(nn);
            end
            if data3(nn)==sofar && sum(data3(nn+1:end))==0
                position=speeds(nn);
                value=total;
                cutoff=total;
            elseif data3(nn)==sofar && nn==1
                position=speeds(nn);
                value=total/2;
                cutoff=total/2;
            else
                cutoff=total/2;
                a=nn-1; b=nn;
                before=sum(data3(1:a));
                after=sum(data3(1:b));
                slope=(after-before)/(speeds(b)-speeds(a));
                position=(slope*speeds(b)-after+total/2)/slope;
                c=find(speeds<position);
                c=c(end);
                cc=c+1;
                slope2=(data3(cc)-data3(c))/(speeds(cc)-speeds(c));
                value=data3(cc)-slope2*speeds(cc)+position*slope2;
            end
            summed_coll(1)=data3(1);
            for rr=2:length(data3)
                summed_coll(rr)=summed_coll(rr-1)+data3(rr);
            end
        else
            position=0;
            value=0;
            summed_coll=repmat(0,length(data3),1);
            cutoff=total/2;
        end
        
        figure
        set(gcf,'position',[1 31 1600 794])
        subplot(2,4,[3 4])
        plot(avRate,'-k')
        hold on
        for j=1:6
            plot([startpoints(j) startpoints(j)],[0 maxval+1],'-r')
        end
        axis([0 max(startpoints)+50 0 maxval+1])
        subplot(2,4,1)
        plot(speeds,summed_coll,'-d','color','k','markerfacecolor','k','markersize',3);
        hold on
        plot(position,cutoff,'*k','markersize',5)
        xlabel('speed (mm/s)')
        ylabel('summed response')
        axis([0 max(speeds)+1 0 max(summed_coll)+1])
        title(['Black: total summed activity: ',num2str(total),' (',num2str(total/avbackground),')'])
        subplot(2,4,2)
        plot(speeds,data3,'-dk','markerfacecolor','black','markersize',7);
        hold on
        plot(position,value,'*k','markersize',7);
        xlabel('speed (mm/s)')
        ylabel('peak spike frequency')
        axis([0 max(speeds)+1 0 max(data3)+1])
        title(['50% of activity after ',num2str(position),' mm/s'])
        info{currID,14}=position;
        
        keep_position=position;
        keep_total=total;
        keep_summed=summed_coll;
        keep_data=data3;
        keep_value=value;
        keep_cutoff=cutoff;
        keep_background=avbackground;
        keep_maxval=maxval;
        keep_maxid=maxid;
        
        
        
        clear data3
        avRate=meanrateW;
        avbackground=mean(avRate(1:startpoints(1)));
        for j=1:6
            [data3(j) data4(j)]=max(avRate(startpoints(j):startpoints(j+1)));
        end
        [maxval maxid]=max(data3);
        data3=data3-avbackground;
        
        
        total=sum(data3); sofar=0; nn=0;
        if total>0
            while sofar<total/2
                nn=nn+1;
                sofar=sofar+data3(nn);
            end
            if data3(nn)==sofar && sum(data3(nn+1:end))==0
                position=speeds(nn);
                value=total;
                cutoff=total;
            elseif data3(nn)==sofar && nn==1
                position=speeds(nn);
                value=total/2;
                cutoff=total/2;
            else
                cutoff=total/2;
                a=nn-1; b=nn;
                before=sum(data3(1:a));
                after=sum(data3(1:b));
                slope=(after-before)/(speeds(b)-speeds(a));
                position=(slope*speeds(b)-after+total/2)/slope;
                c=find(speeds<position);
                c=c(end);
                cc=c+1;
                slope2=(data3(cc)-data3(c))/(speeds(cc)-speeds(c));
                value=data3(cc)-slope2*speeds(cc)+position*slope2;
            end
            summed_coll(1)=data3(1);
            for rr=2:length(data3)
                summed_coll(rr)=summed_coll(rr-1)+data3(rr);
            end
        else
            position=0;
            value=0;
            summed_coll=repmat(0,length(data3),1);
            cutoff=total/2;
        end
        
        subplot(2,4,[7 8])
        plot(avRate,'-k')
        hold on
        for j=1:6
            plot([startpoints(j) startpoints(j)],[0 maxval+1],'-r')
        end
        axis([0 max(startpoints)+50 0 maxval+1])
        subplot(2,4,5)
        plot(speeds,summed_coll,'-d','color','k','markerfacecolor','k','markersize',3);
        hold on
        plot(position,cutoff,'*k','markersize',7)
        xlabel('speed (mm/s)')
        ylabel('summed response')
        axis([0 max(speeds)+1 0 max(summed_coll)+1])
        title(['White: total summed activity: ',num2str(total),' (',num2str(total/avbackground),')'])
        subplot(2,4,6)
        plot(speeds,data3,'-dk','markerfacecolor','black','markersize',7);
        hold on
        plot(position,value,'*k','markersize',5);
        xlabel('speed (mm/s)')
        ylabel('peak spike frequency')
        axis([0 max(speeds)+1 0 max(data3)+1])
        title(['50% of activity after ',num2str(position),' mm/s'])
        info{currID,15}=position;
        
        
        
        if keep_background==0 && avbackground==0
            testval=keep_total>total;
        else
            testval=keep_total/keep_background>total/avbackground;
        end
        if testval==0
            testval=2;
        end
        
        %     choice = questdlg('which one is better?', ...
        %         '6vel', ...
        %         '1','2','none',int2str(testval));
        %     switch choice
        %         case '1'
        %             info{currID,16}=1;
        %             subplot(2,4,1)
        %             plot(speeds,keep_summed,'-d','color',[0 0.6 0.2],'markerfacecolor',[0 0.6 0.2],'markersize',3);
        %             plot(keep_position,keep_cutoff,'*k','markersize',12)
        %         case '2'
        %             info{currID,16}=2;
        %             subplot(2,4,5)
        %             plot(speeds,summed_coll,'-d','color',[0 0.6 0.2],'markerfacecolor',[0 0.6 0.2],'markersize',3,'linewidth',3);
        %             hold on
        %             plot(position,cutoff,'*k','markersize',12)
        %         case 'none'
        %             info{currID,16}=0;
        %             info{currID,14}=0;
        %             info{currID,15}=0;
        %     end
        saveas(gcf,fullfile(superpath,'6vel',[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
        close
    else
        info{currID,14}=0;
        info{currID,15}=0;
    end
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
save(fullfile(superpath,[Dtype,'_vel6Black']),'vel6_ratesB')
save(fullfile(superpath,[Dtype,'_vel6White']),'vel6_ratesW')
%% 9b) check whether real responses
load(fullfile(superpath,[Dtype,'_info']))

for i=1:size(info,1)
    if ~isempty(info{i,14})
        if info{i,14}>0 || info{i,15}>0
            currd=info{i,1};
            
            curru=info{i,2};
            winopen(fullfile(superpath,'6vel',[Dtype,'_',int2str(i),'_',currd,'_',curru,'.bmp']))
            str={'black','white'};
            [s,v] = listdlg('PromptString','Select responses',...
                'SelectionMode','multiple',...
                'ListString',str);
            
            if isempty(find(s==1))
                info{i,14}=0;
            end
            if isempty(find(s==2))
                info{i,15}=0;
            end
        end
    end
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 9c) compute speed only for judged units (better background) => do this again but with subtracting
% background + 2STD
% col16) #maxid for 6vel black (speed for which it's best)
% col17) #maxid for 6vel white (optimal speed)


load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
if ~exist(fullfile(superpath,'6vel_new'),'dir')
    mkdir(fullfile(superpath,'6vel_new'))
end
names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[0.2 1 2 4 8 16];

vel6_newB=cell(length(info),1);
vel6_newW=cell(length(info),1);
togo=[];
for i=1:length(info)
    if ~isempty(info{i,14})
        if info{i,14}>0 || info{i,15}>0
            togo=[togo;i];
        end
    end
end
for i=1:length(togo)
    i
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    currstart=info{currID,4};
    currstop=info{currID,5};
    goodcond=[];
    if ~isempty(info{currID,14})
        if info{currID,14}>0
            goodcond=[goodcond 1];
        end
        if info{currID,15}>0
            goodcond=[goodcond 2];
        end
    end
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
    
    black=[]; white=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'bar_6veloc_downup_black'))
            ftype=1;
            black=[black;1];
        else
            black=[black;0];
        end
        if ~isempty(regexp(hekalist{f},'bar_6veloc_downup_white'))
            ftype=2;
            white=[white;1];
        else
            white=[white;0];
        end
    end
    goodfiles=currstart:currstop;
    goodfilesB=goodfiles(black==1);
    goodfilesW=goodfiles(white==1);
    
    if ~isempty(goodfilesB)
        first=goodfilesB(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        cons=prot(1:end,2);
        changes=find(cons==-1);
        startpoints=round(prot(changes(1:7),1));
        tt=round(prot(end,1)+1000);
        
        allrateB=[];
        for f=goodfilesB
            spikes=unit{1,2}{f,2};
            rate=MEA_spikerates(spikes,40,tt);
            allrateB=[allrateB;rate];
        end
        meanrateB=mean(allrateB);
        vel6_newB{currID,1}=meanrateB;
        
        allrateW=[];
        for f=goodfilesW
            spikes=unit{1,2}{f,2};
            rate=MEA_spikerates(spikes,40,tt);
            allrateW=[allrateW;rate];
        end
        meanrateW=mean(allrateW);
        vel6_newW{currID,1}=meanrateW;
        
        if find(goodcond==1)
            clear data3
            avRate=meanrateB;
            avbackground=mean(avRate(200:startpoints(1)));
            
            for j=1:6
                [data3(j) data4(j)]=max(avRate(startpoints(j):startpoints(j+1)));
            end
            [maxval maxid]=max(data3);
            data3=data3-avbackground;
            data3=(data3-min(data3))/(max(data3)-min(data3));
            
            
            total=sum(data3); sofar=0; nn=0;
            if total>0
                while sofar<total/2
                    nn=nn+1;
                    sofar=sofar+data3(nn);
                end
                if data3(nn)==sofar && sum(data3(nn+1:end))==0
                    position=speeds(nn);
                    value=total;
                    cutoff=total;
                elseif data3(nn)==sofar && nn==1
                    position=speeds(nn);
                    value=total/2;
                    cutoff=total/2;
                else
                    cutoff=total/2;
                    a=nn-1; b=nn;
                    before=sum(data3(1:a));
                    after=sum(data3(1:b));
                    slope=(after-before)/(speeds(b)-speeds(a));
                    position=(slope*speeds(b)-after+total/2)/slope;
                    c=find(speeds<position);
                    c=c(end);
                    cc=c+1;
                    slope2=(data3(cc)-data3(c))/(speeds(cc)-speeds(c));
                    value=data3(cc)-slope2*speeds(cc)+position*slope2;
                end
                summed_coll(1)=data3(1);
                for rr=2:length(data3)
                    summed_coll(rr)=summed_coll(rr-1)+data3(rr);
                end
            else
                position=0;
                value=0;
                summed_coll=repmat(0,length(data3),1);
                cutoff=total/2;
            end
            
            figure
            set(gcf,'position',[1 31 1600 794])
            subplot(2,4,[3 4])
            plot(avRate,'-k')
            hold on
            for j=1:6
                plot([startpoints(j) startpoints(j)],[0 maxval+1],'-r')
            end
            axis([0 max(startpoints)+50 0 maxval+1])
            subplot(2,4,1)
            plot(speeds,summed_coll,'-d','color','k','markerfacecolor','k','markersize',3);
            hold on
            plot(position,cutoff,'*k','markersize',5)
            xlabel('speed (mm/s)')
            ylabel('summed response')
            axis([0 max(speeds)+1 0 max(summed_coll)+1])
            title(['Black: total summed activity: ',num2str(total),' (',num2str(total/avbackground),')'])
            subplot(2,4,2)
            plot(speeds,data3,'-dk','markerfacecolor','black','markersize',7);
            hold on
            plot(position,value,'*k','markersize',7);
            xlabel('speed (mm/s)')
            ylabel('peak spike frequency')
            axis([0 max(speeds)+1 0 1.1])
            title(['50% of activity after ',num2str(position),' mm/s'])
            info{currID,16}=position;
            
            keep_position=position;
            keep_total=total;
            keep_summed=summed_coll;
            keep_data=data3;
            keep_value=value;
            keep_cutoff=cutoff;
            keep_background=avbackground;
            keep_maxval=maxval;
            keep_maxid=maxid;
        else
            info{currID,16}=0;
        end
        
        
        if find(goodcond==2)
            clear data3
            avRate=meanrateW;
            avbackground=mean(avRate(200:startpoints(1)));
            for j=1:6
                [data3(j) data4(j)]=max(avRate(startpoints(j):startpoints(j+1)));
            end
            [maxval maxid]=max(data3);
            data3=data3-avbackground;
            data3=(data3-min(data3))/(max(data3)-min(data3));
            
            total=sum(data3); sofar=0; nn=0;
            if total>0
                while sofar<total/2
                    nn=nn+1;
                    sofar=sofar+data3(nn);
                end
                if data3(nn)==sofar && sum(data3(nn+1:end))==0
                    position=speeds(nn);
                    value=total;
                    cutoff=total;
                elseif data3(nn)==sofar && nn==1
                    position=speeds(nn);
                    value=total/2;
                    cutoff=total/2;
                else
                    cutoff=total/2;
                    a=nn-1; b=nn;
                    before=sum(data3(1:a));
                    after=sum(data3(1:b));
                    slope=(after-before)/(speeds(b)-speeds(a));
                    position=(slope*speeds(b)-after+total/2)/slope;
                    c=find(speeds<position);
                    c=c(end);
                    cc=c+1;
                    slope2=(data3(cc)-data3(c))/(speeds(cc)-speeds(c));
                    value=data3(cc)-slope2*speeds(cc)+position*slope2;
                end
                summed_coll(1)=data3(1);
                for rr=2:length(data3)
                    summed_coll(rr)=summed_coll(rr-1)+data3(rr);
                end
            else
                position=0;
                value=0;
                summed_coll=repmat(0,length(data3),1);
                cutoff=total/2;
            end
            
            subplot(2,4,[7 8])
            plot(avRate,'-k')
            hold on
            for j=1:6
                plot([startpoints(j) startpoints(j)],[0 maxval+1],'-r')
            end
            axis([0 max(startpoints)+50 0 maxval+1])
            subplot(2,4,5)
            plot(speeds,summed_coll,'-d','color','k','markerfacecolor','k','markersize',3);
            hold on
            plot(position,cutoff,'*k','markersize',7)
            xlabel('speed (mm/s)')
            ylabel('summed response')
            axis([0 max(speeds)+1 0 max(summed_coll)+1])
            title(['White: total summed activity: ',num2str(total),' (',num2str(total/avbackground),')'])
            subplot(2,4,6)
            plot(speeds,data3,'-dk','markerfacecolor','black','markersize',7);
            hold on
            plot(position,value,'*k','markersize',5);
            xlabel('speed (mm/s)')
            ylabel('peak spike frequency')
            axis([0 max(speeds)+1 0 1.1])
            title(['50% of activity after ',num2str(position),' mm/s'])
            info{currID,17}=position;
        else
            info{currID,17}=0;
        end
        
        
        
        
        saveas(gcf,fullfile(superpath,'6vel_new',[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
        close
    else
        info{currID,16}=0;
        info{currID,17}=0;
    end
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
save(fullfile(superpath,[Dtype,'_vel6_new_Black']),'vel6_newB')
save(fullfile(superpath,[Dtype,'_vel6_new_White']),'vel6_newW')
%% 9d) PLOT histograms for speed tuning
clear
close all
% superpath='S:\data\Katja\allstuff';
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');


names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[0.2; 1; 2; 4; 8; 16];
types='wph';
cols='bgr';
titles={'mouse','pig','human'};
figure(1)
figure(2)
for t=1:length(types)
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,16})
            if info{i,16}>0
                goodB=[goodB;i];
            end
            if info{i,17}>0
                goodW=[goodW;i];
            end
        end
    end
    
    black=cell2mat(info(goodB,16));
    if t==1
        mblack=black;
    elseif t==2
        pblack=black;
    else
        hblack=black;
    end
    white=cell2mat(info(goodW,17));
    if t==1
        mwhite=white;
    elseif t==2
        pwhite=white;
    else
        hwhite=white;
    end
    black=[black;-1;17];
    white=[white;-1;17];
    
    figure(1)
    subplot(3,2,t*2-1)
    hist(black,38)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
    med=median(black(1:end-2));
    st=std(black(1:end-2));
    conf = 1.96*st/sqrt(length(black)-2);
    title([titles{t}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(black)-2),')'])
    hold on
    plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
    plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
    yy=get(gca,'ylim');
    plot([med-conf med-conf],[0 yy(2)+1],'--k')
    plot([med+conf med+conf],[0 yy(2)+1],'--k')
    %     set(gca,'xtick',speeds)
    % set(gca,'xticklabel',speeds)
    xlabel('speed (mm/s)')
    axis([0 16 0 yy(2)+1])
    
    subplot(3,2,[2 4])
    c=cdfplot(black(1:end-2))
    set(c,'color',cols(t),'linestyle','-','linewidth',2)
    hold on
    xlabel('speed (mm/s)')
    ylabel('probability')
    title('black bar')
    
    figure(2)
    subplot(3,2,t*2-1)
    hist(white,38)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
    med=median(white(1:end-2));
    st=std(white(1:end-2));
    conf = 1.96*st/sqrt(length(black)-2);
    title([titles{t}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(white)-2),')'])
    hold on
    plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
    plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
    xlabel('speed (mm/s)')
    yy=get(gca,'ylim');
    plot([med-conf med-conf],[0 yy(2)+1],'--k')
    plot([med+conf med+conf],[0 yy(2)+1],'--k')
    axis([0 16 0 yy(2)+1])
    %      set(gca,'xtick',speeds)
    %     set(gca,'xticklabel',speeds)
    
    subplot(3,2,[2 4])
    c=cdfplot(white(1:end-2))
    set(c,'color',cols(t),'linestyle','-','linewidth',2)
    hold on
    xlabel('speed (mm/s)')
    ylabel('probability')
    title('white bar')
    %      set(gca,'xtick',speeds)
    % set(gca,'xticklabel',speeds)
end
figure(1)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'speed_black.bmp'))
saveas(gcf,fullfile(savepath,'speed_black.fig'))
close
figure(2)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'speed_white.bmp'))
saveas(gcf,fullfile(savepath,'speed_white.fig'))
close
%% 10) compute spatiotemporal
% col18) fourier peak
load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
if ~exist(fullfile(superpath,'DG'),'dir')
    mkdir(fullfile(superpath,'DG'))
end

DGrates=cell(length(info),1);
DGfouriers=cell(length(info),1);
stimInfo=cell(length(info),1);
togo=find(goodones>0);
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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    thiscell=[];
    fouriers=[];
    rates=[];
    clear f2take
    for g=1:length(DGlist)
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            rates=[rates;meanrate];
            
            stimFreq=freqs(g);
            
            if max(meanrate)>0
                meanF=mean(meanrate);
                meanFR=round(meanF*10)/10;
                % fourier
                NFFT=2^nextpow2(numel(meanrate));
                Y=fft(meanrate,NFFT)/numel(meanrate);
                f=1000/2*linspace(0,1,NFFT/2+1);
                if ~exist('f2take','var') && ~isempty(find(f>0))
                    f2take=f;
                end
                c=2*abs(Y(1:NFFT/2+1));
                if g>1 && size(fouriers,2)>length(c)
                    fouriers=fouriers(:,1:length(c));
                end
                fouriers=[fouriers;c];
                
                ampPeak=max(c(4:end));
                ampPeakR=round(ampPeak);
                posF=find(c==ampPeak);
                posF=posF(1);
                freq=f(posF);
                findfreq1=find(f<stimFreq);
                findfreq1=findfreq1(end)-1;
                if findfreq1<=0
                    findfreq1=4;
                end
                findfreq2=max(c(findfreq1:findfreq1+5)); %recorded amp at stimfreq
                pos1=find(c==findfreq2);
                pos1=f(pos1);
            else
                findfreq2=0;
                pos1=0;
                if g==1
                    fouriers=zeros(1,10000);
                else
                    fouriers=[fouriers;zeros(1,size(fouriers,2))];
                end
            end
            thiscell=[thiscell;findfreq2 pos1];
            size(rates,1)
        end
    end
    
    
    if gofor==1
        normthis=thiscell(:,1)/max(thiscell(:,1));
        
        
        
        zi=reshape(normthis,length(unique(freqs)),length(unique(spats)));
        zi=zi';
        zi=flipud(zi);
        figure
        set(gcf,'position',[1 31 1600 794])
        subplot(3,2,[2 4])
        imagesc(1:length(unique(freqs)),1:length(unique(spats)),zi)
        set(gca,'YTick',1:length(unique(spats)))
        set(gca,'XTick',1:length(unique(freqs)))
        set(gca,'XTickLabel',unique(freqs))
        xlabel('Hz')
        set(gca,'YTickLabel',flipud(unique(spats)))
        ylabel('period (um)')
        
        [listed ids]=sort(normthis,'descend');
        
        if exist('f2take','var')
            tis={'best','second','third'};
            for l=1:3
                subplot(3,2,l*2-1)
                plot(f2take,fouriers(ids(l),:))
                hold on
                plot(thiscell(ids(l),2),thiscell(ids(l),1),'*r')
                axis([0 16 0 max(fouriers(ids(l),:))+1])
                title(tis{l})
            end
            
            subplot(3,2,6)
            plot(rates(ids(1),:))
            title('best')
        end
        
        for r=1:size(rates,1)
            DGrates{currID,r}=rates(r,:);
            DGfouriers{currID,r}=fouriers(r,:);
        end
        stimInfo{currID,1}=freqs;
        stimInfo{currID,2}=spats;
        
        
        choice = questdlg('does the cell respond?', ...
            'DG', ...
            'yes','no','yes');
        switch choice
            case 'yes'
                info{currID,18}=thiscell(ids(1),1);
            case 'OFF'
                info{currID,18}=0;
        end
        
        saveas(gcf,fullfile(superpath,'DG',[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
        close
    else
        info{currID,18}=0;
    end
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
save(fullfile(superpath,[Dtype,'_DG']),'stimInfo','DGfouriers','DGrates')
%% 10) !!!!!!!!!!!!!!!!!!!!!!!!work on this for human data!!
% col18) fourier peak
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_DG']))
cd(superpath)


DGrates=cell(length(info),1);
togo=find(goodones>0);
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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    
    rates=[];
    for g=1:length(DGlist)
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            rates=[rates;meanrate];
            
            stimFreq=freqs(g);
            
            
        end
    end
    
    
    if gofor==1
        
        for r=1:size(rates,1)
            DGrates{currID,r}=rates(r,:);
        end
        
        
    end
end
save(fullfile(superpath,[Dtype,'_DG']),'stimInfo','DGfouriers','DGrates')
%% 10b) list with good DG units
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_DG']))
goodDG=[];
for i=1:length(info)
    if ~isempty(info{i,18})
        if info{i,18}>0
            goodDG=[goodDG;1];
            if Dtype=='p' && ~isempty(regexp(info{i,1},'20140723d'))
                goodDG=[goodDG(1:end-1);0];
            end
        else
            goodDG=[goodDG;0];
        end
    else
        goodDG=[goodDG;0];
    end
end
save(fullfile(superpath,[Dtype,'_goodDG']),'goodDG','stimInfo')
%% 10c) forgot to save all peaks....redo it for the judged units
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)


DGpeaks=cell(length(info),1);
togo=find(goodDG==1);
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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    thiscell=[];
    fouriers=[];
    rates=[];
    clear f2take
    for g=1:length(DGlist)
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            
            stimFreq=freqs(g);
            
            if max(meanrate)>0
                meanF=mean(meanrate);
                meanFR=round(meanF*10)/10;
                % fourier
                NFFT=2^nextpow2(numel(meanrate));
                Y=fft(meanrate,NFFT)/numel(meanrate);
                f=1000/2*linspace(0,1,NFFT/2+1);
                if ~exist('f2take','var') && ~isempty(find(f>0))
                    f2take=f;
                end
                c=2*abs(Y(1:NFFT/2+1));
                
                
                ampPeak=max(c(4:end));
                ampPeakR=round(ampPeak);
                posF=find(c==ampPeak);
                posF=posF(1);
                freq=f(posF);
                findfreq1=find(f<stimFreq);
                findfreq1=findfreq1(end)-1;
                if findfreq1<=0
                    findfreq1=4;
                end
                findfreq2=max(c(findfreq1:findfreq1+5)); %recorded amp at stimfreq
                pos1=find(c==findfreq2);
                pos1=f(pos1);
            else
                findfreq2=0;
                pos1=0;
                
            end
            thiscell=[thiscell;findfreq2 pos1]; %max at freq, corresponding freq (exactly)
            size(rates,1)
        end
    end
    DGpeaks{currID}=thiscell;
end
save(fullfile(superpath,[Dtype,'_DGpeaks']),'DGpeaks')
%% --10d) PLOT spatiotemporal
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaks']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            allDataTemp{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    
    figure(1)
    data=sumSpat;
    subplot(2,2,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,2,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,2,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,2,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,2,s)
    axis([0 7 -inf inf])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq.fig'))
close
figure(2)
for s=1:4
    subplot(2,2,s)
    axis([0 5 -inf inf])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq.fig'))
close
%% 11) chirp
% col19) chirp yes/no
load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
if ~exist(fullfile(superpath,'chirp'),'dir')
    mkdir(fullfile(superpath,'chirp'))
end

window=500;
dT=1/60*1000;
f0=0.05;
k=1;
t=0:dT:8000;
t=t/1000;
f1=f0+k*t;
f2=[];
for i=1:length(f1)
    if mod(i,3)==0
        curr=repmat(f1(i),16,1);
    else
        curr=repmat(f1(i),17,1);
    end
    f2=[f2;curr];
end

togo=find(goodones>0);
chirpFreq=cell(length(info),1);
chirpContr=cell(length(info),1);
chirpMod=cell(length(info),2);
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
    
    
    if ~isempty(goodfiles)
        
        allrate=[];
        for f=goodfiles
            spikes=unit{1,2}{f,2};
            rate=MEA_spikerates(spikes,10,25000);
            allrate=[allrate;rate];
        end
        meanrate=mean(allrate);
        allrate=[];
        for f=goodfiles
            spikes=unit{1,2}{f,2};
            rate=MEA_spikerates(spikes,40,25000);
            allrate=[allrate;rate];
        end
        meanrate2=mean(allrate);
        
        figure
        set(gcf,'position',[203 177 1148 630])
        subplot(2,1,1)
        plot(meanrate)
        subplot(2,1,2)
        plot(meanrate2)
        choice = questdlg('does the cell respond?', ...
            'DG', ...
            'yes','no','yes');
        switch choice
            case 'yes'
                info{currID,19}=1;
            case 'OFF'
                info{currID,19}=0;
        end
        saveas(gcf,fullfile(superpath,'chirp',[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
        close
        
        chirpFreq{currID,1}=meanrate(1:11500);
        chirpContr{currID,1}=meanrate(11500:20500);
        
        if ~isempty(regexp(choice,'yes'))
            backR=mean(meanrate(500:1000));
            keepmean=meanrate;
            meanrate=meanrate-backR;
            absrate=abs(meanrate);
            cnt=1;
            clear param freqrange
            for fc=2500:50:11500-window
                param(cnt)=sum(absrate(fc:fc+window))/window; %modulation strength
                if fc>=3000 && length(f2)>=fc-2999+window
                    freqrange(cnt)=mean(f2(fc-2999:fc-2999+window));
                end
                cnt=cnt+1;
            end
            chirpMod{currID,1}=param;
            chirpMod{currID,2}=freqrange;
        end
        
    end
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
save(fullfile(superpath,[Dtype,'_chirp']),'chirpMod','chirpFreq','chirpContr')
%% 11b) correct for different herz (exclude) and calculated modulation for contrast part => do for m and h
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_chirp']))
cd(superpath)

window=500;
dT=1/60*1000;
f0=0.05;
k=1;
t=0:dT:8000;
t=t/1000;
f1=f0+k*t;
f2=[];
for i=1:length(f1)
    if mod(i,3)==0
        curr=repmat(f1(i),16,1);
    else
        curr=repmat(f1(i),17,1);
    end
    f2=[f2;curr];
end

togo=[];
for i=1:length(info)
    if ~isempty(info{i,19})
        if info{i,19}>0
            togo=[togo;i];
        end
    end
end

chirpModContr=cell(length(info),1);
counter=0;
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
    
    prot=read_header_field_heka(hpath,hekalist{goodfiles(1)},'Stimulus Protocol');
    if prot(end,1)>24000
        goodfiles=[];
    end
    
    
    if ~isempty(goodfiles)
        %         counter=counter+1
        allrate=[];
        for f=goodfiles
            spikes=unit{1,2}{f,2};
            rate=MEA_spikerates(spikes,40,25000);
            allrate=[allrate;rate];
        end
        meanrate=mean(allrate);
        
        backR=mean(meanrate(500:1000));
        keepmean=meanrate;
        meanrate=meanrate-backR;
        absrate=abs(meanrate);
        cnt=1;
        clear param freqrange
        for fc=2500:50:11500-window
            param(cnt)=sum(absrate(fc:fc+window))/window; %modulation strength
            if fc>=3000 && length(f2)>=fc-2999+window
                freqrange(cnt)=mean(f2(fc-2999:fc-2999+window));
            end
            cnt=cnt+1;
        end
        chirpMod{currID,1}=param;
        chirpMod{currID,2}=freqrange;
        
        cnt=1;
        clear param freqrange
        for fc=11500:50:20500-window
            param(cnt)=sum(absrate(fc:fc+window))/window; %modulation strength
            cnt=cnt+1;
        end
        chirpModContr{currID,1}=param;
    else
        info{currID,19}=0;
        chirpMod{currID,1}={};
        chirpMod{currID,2}={};
    end
    
end

save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
save(fullfile(superpath,[Dtype,'_chirp']),'chirpMod','chirpFreq','chirpContr','chirpModContr')
%% 11c) extract max. modulation strength for certain frequency ranges
% col20) max for 1-2 Hz
% col21) max for 3.8-4.2 Hz
% col22) max for 5.8-6.2 Hz
% col23) max for 7.3-7.7 Hz

load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_chirp']))
cd(superpath)

startF=[1 3.8 5.8 7.3];
stopF=[2 4.2 6.2 7.7];

togo=[];
for i=1:length(info)
    if ~isempty(info{i,19})
        if info{i,19}==1
            togo=[togo;i];
        end
    end
end

for i=1:length(togo)
    currID=togo(i);
    
    param=chirpMod{currID,1};
    freqrange=chirpMod{currID,2};
    param=param(11:end);
    freqrange=freqrange(11:end);
    
    tmp1=find(freqrange<1);
    tmp1=tmp1(end);
    tmp2=length(freqrange);
    toconsider=param(tmp1:tmp2);
    backnew=mean(param(1:5));
    stnew=std(param(1:5));
    
    norm=param(tmp1:tmp2)/max(param(tmp1:tmp2));
    newfreqr=freqrange(tmp1:tmp2);
    
    for j=1:length(startF)
        tmp1=find(newfreqr<startF(j));
        tmp1=tmp1(end);
        tmp2=find(newfreqr<=stopF(j));
        tmp2=tmp2(end);
        info{currID,19+j}=max(norm(tmp1:tmp2));
    end
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% 11d) PLOT

clear
close all
% superpath='S:\data\Katja\allstuff';
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

types='wph';
cols='bgr';
titles={'mouse','pig','human'};
texts={'1-2','3.8-4.2','5.8-6.2','7.3-7.7'};

figure(1)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_chirp']))
    togo=[];
    for i=1:length(info)
        if ~isempty(info{i,19})
            if info{i,19}==1
                togo=[togo;i];
            end
        end
    end
    
    figure(1)
    coll=[];
    for i=1:length(togo)
        currID=togo(i);
        currdata=chirpMod{currID,1};
        currFr=chirpMod{currID,2}(7:end);
        currdata=currdata(7:length(currFr)+6)/max(currdata(7:length(currFr)+6));
        coll=[coll;currdata];
        subplot(3,3,(tt-1)*3+1)
        plot(currFr,currdata,'-','color',[0.6 0.6 0.6])
        hold on
    end
    plot(currFr,median(coll,1),'-','color',cols(tt),'linewidth',3)
    fill([currFr fliplr(currFr)],[median(coll)-std(coll,0,1) fliplr(median(coll)+std(coll,0,1))],cols(tt),'facealpha',0.6)
    axis tight
    %     xx=get(gca,'xlim');
    %     set(gca,'xtick',0:xx(2)/9:xx(2))
    %     stepping=length(currFr)/10;
    %     set(gca,'xticklabel',[round(currFr(stepping/2:length(currFr)/10:length(currFr))*10)/10 xx(2)])
    xlabel('Hz')
    title('freq. modulation')
    subplot(3,3,7)
    plot(currFr,median(coll,1),'-','color',cols(tt),'linewidth',3)
    hold on
    axis tight
    %     xx=get(gca,'xlim');
    %     set(gca,'xtick',1:xx(2)/10:xx(2))
    %     set(gca,'xticklabel',round(currFr(1:length(currFr)/10:length(currFr))*10)/10)
    
    alldata=cell2mat(info(togo,20:23));
    subplot(3,3,(tt-1)*3+2)
    plot(median(alldata),'s-','color',cols(tt),'linewidth',2)
    hold on
    fill([1:length(median(alldata)) length(median(alldata)):-1:1],[median(alldata)-std(alldata,0,1) fliplr(median(alldata)+std(alldata,0,1))],cols(tt),'facealpha',0.3)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',texts)
    xlabel('Hz')
    set(gca,'ytick',0:0.2:1.2)
    set(gca,'ytick',0:0.2:1.2)
    title([titles{tt},': freq. modulation (total: ',int2str(size(alldata,1)),')'])
    axis([1 4 0 1.2])
    
    subplot(3,3,8)
    plot(median(alldata),'s-','color',cols(tt),'linewidth',2)
    hold on
    axis([1 4 0 1.2])
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',texts)
    xlabel('Hz')
    set(gca,'ytick',0:0.2:1.2)
    set(gca,'ytick',0:0.2:1.2)
    
    
    
    subplot(3,3,(tt-1)*3+3)
    coll=[];
    for i=1:length(togo)
        currID=togo(i);
        currdata=chirpModContr{currID,1};
        if ~isempty(find(currdata)>0)
            currdata=currdata/max(currdata);
            coll=[coll;currdata];
            plot(currdata,'-','color',[0.6 0.6 0.6])
        end
        hold on
    end
    plot(median(coll,1),'-','color',cols(tt),'linewidth',3)
    fill([1:length(median(coll)) length(median(coll)):-1:1],[median(coll)-std(coll,0,1) fliplr(median(coll)+std(coll,0,1))],cols(tt),'facealpha',0.6)
    title('contrast modulation')
    axis tight
    
    subplot(3,3,9)
    plot(median(coll,1),'-','color',cols(tt),'linewidth',3)
    hold on
    axis tight
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'chirp_modulation.bmp'))
saveas(gcf,fullfile(savepath,'chirp_modulation.fig'))
close
%% read out some values
clear
close all

Dtype='w';
% superpath='S:\data\Katja\allstuff';
superpath='C:\Users\Katja\Desktop\all_species_temp';
load(fullfile(superpath,[Dtype,'_info']))

good=length(find(goodones==1));
fullfield=cell2mat(info(:,8));
fullfield=length(find(fullfield~=0));
LF=cell2mat(info(:,9));
LF=length(find(LF~=0));
DS=cell2mat(info(:,12));
DS=length(find(DS==1));

vel61=cell2mat(info(:,16));
vel61=find(vel61>0);
vel62=cell2mat(info(:,17));
vel62=find(vel62>0);
vel6=length(unique([vel61;vel62]));
DG=cell2mat(info(:,18));
DG=length(find(DG>0));
chirp=cell2mat(info(:,19));
chirp=length(find(chirp>0));
%% check for ON-OFF bias
clear
close all
% superpath='S:\data\Katja\allstuff';
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    
    good=length(find(goodones==1));
    fullfield=cell2mat(info(:,8));
    
    on=length(find(fullfield==1));
    off=length(find(fullfield==-1));
    all=length(find(fullfield~=0));
    ON=100/all*on;
    OFF=100/all*off;
    
    %per exp
    figure(1)
    switch Dtype
        
        case 'h'
            % Human data
            dates={'20120210','20120224','20120320c','20120608b','20130124a','20130211a',...,
                '20130211c','20130211d','20130510a','20130510d','20130510e','20130718b','20130718c','20131210a','20131217f'};
        case 'p'
            % Pig data
            dates={'20121127','20121211a','20130613a','20130613b','20130613c','20140723d','20120213','20120319a','20120319b'};
        case 'm'
            % Mouse data
            dates={'20120126','20120127','20120221','20120222'};
    end
    onoff=zeros(length(dates),3);
    for d=1:length(dates)
        curr=dates{d};
        ids=[];
        good=find(goodones==1);
        for i=1:length(good)
            if ~isempty(regexp(info{good(i),1},curr))
                ids=[ids;good(i)];
            end
        end
        fullfield=cell2mat(info(ids,8));
        fullfield2=cell2mat(info(ids,9));
        on=length(find(fullfield==1));
        off=length(find(fullfield==-1));
        ono=length(find(fullfield==2));
        tmp1=find(fullfield~=0);
        tmp2=find(fullfield2~=0);
        all=[tmp1;tmp2];
        all=length(unique(all));
        %         all=length(find(fullfield~=0));
        % ON=100/all*on;
        % OFF=100/all*off;
        onoff(d,1)=on;
        onoff(d,2)=off;
        onoff(d,3)=all;
        onoff(d,4)=ono;
    end
    if tt<3
        subplot(2,2,tt)
    else
        subplot(2,2,[3 4])
    end
    data=[100./onoff(:,3).*onoff(:,1) 100./onoff(:,3).*onoff(:,2) 100./onoff(:,3).*onoff(:,4); 100./sum(onoff(:,3)).*sum(onoff(:,1)) 100./sum(onoff(:,3)).*sum(onoff(:,2)) 100./sum(onoff(:,3)).*sum(onoff(:,4))];
    data(:,4)=100-data(:,1)-data(:,2)-data(:,3);
    for da=1:size(data,1)-1
        if data(da,1)>0
            rectangle('position',[da-0.35,0,0.7,data(da,1)],'facecolor','c')
        end
        hold on
        if data(da,4)>0
            rectangle('position',[da-0.35,data(da,1),0.7,data(da,4)],'facecolor','w')
        end
        if data(da,3)>0
            rectangle('position',[da-0.35,data(da,1)+data(da,4),0.7,data(da,3)],'facecolor','m')
        end
        if data(da,2)>0
            rectangle('position',[da-0.35,data(da,1)+data(da,3)+data(da,4),0.7,data(da,2)],'facecolor','k')
        end
        
        
    end
    rectangle('position',[da+2-0.35,0,0.7,data(da+1,1)],'facecolor','c')
    hold on
    rectangle('position',[da+2-0.35,data(da+1,1),0.7,data(da+1,4)],'facecolor','w')
    rectangle('position',[da+2-0.35,data(da+1,1)+data(da+1,4),0.7,data(da+1,3)],'facecolor','m')
    rectangle('position',[da+2-0.35,data(da+1,1)+data(da+1,3)+data(da+1,4),0.7,data(da+1,2)],'facecolor','k')%
    %     bar(data,'stacked')
    if tt==3
        rectangle('position',[da+4,37,1,6],'facecolor','c')
        rectangle('position',[da+4,47,1,6],'facecolor','m')
        rectangle('position',[da+4,57,1,6],'facecolor','k')
        text(da+5.2,40,'on')
        text(da+5.2,50,'onoff')
        text(da+5.2,60,'off')
        axis([0 da+7 0 100])
    else
        axis([0 da+3 0 100])
    end
    dates=[dates ' ' 'all'];
    set(gca,'xtick',1:length(dates))
    set(gca,'xticklabel',dates,'fontsize',7)
    title(titles{tt})
    ylabel('% of cells with LF and/or flash response')
    
    
    
    figure(2)
    switch Dtype
        
        case 'h'
            % Human data
            dates={'20120210','20120224','20120320c','20120608b','20130124a','20130211a',...,
                '20130211c','20130211d','20130510a','20130510d','20130510e','20130718b','20130718c','20131210a','20131217f'};
        case 'p'
            % Pig data
            dates={'20121127','20121211a','20130613a','20130613b','20130613c','20140723d','20120213','20120319a','20120319b'};
        case 'm'
            % Mouse data
            dates={'20120126','20120127','20120221','20120222'};
    end
    onoff=zeros(length(dates),3);
    for d=1:length(dates)
        curr=dates{d};
        ids=[];
        good=find(goodones==1);
        for i=1:length(good)
            if ~isempty(regexp(info{good(i),1},curr))
                ids=[ids;good(i)];
            end
        end
        fullfield=cell2mat(info(ids,9));
        fullfield2=cell2mat(info(ids,8));
        on=length(find(fullfield==1));
        off=length(find(fullfield==-1));
        tmp1=find(fullfield~=0);
        tmp2=find(fullfield2~=0);
        all=[tmp1;tmp2];
        all=length(unique(all));
        %         all=length(find(fullfield~=0));
        % ON=100/all*on;
        % OFF=100/all*off;
        onoff(d,1)=on;
        onoff(d,2)=off;
        onoff(d,3)=all;
    end
    if tt<3
        subplot(2,2,tt)
    else
        subplot(2,2,[3 4])
    end
    data=[100./onoff(:,3).*onoff(:,1) 100./onoff(:,3).*onoff(:,2) ; 100./sum(onoff(:,3)).*sum(onoff(:,1)) 100./sum(onoff(:,3)).*sum(onoff(:,2)) ];
    data(:,3)=100-data(:,1)-data(:,2);
    for da=1:size(data,1)-1
        if data(da,1)>0
            rectangle('position',[da-0.35,0,0.7,data(da,1)],'facecolor','c')
        end
        hold on
        if data(da,3)>0
            rectangle('position',[da-0.35,data(da,1),0.7,data(da,3)],'facecolor','w')
        end
        if data(da,2)>0
            rectangle('position',[da-0.35,data(da,1)+data(da,3),0.7,data(da,2)],'facecolor','k')
        end
    end
    rectangle('position',[da+2-0.35,0,0.7,data(da+1,1)],'facecolor','c')
    hold on
    rectangle('position',[da+2-0.35,data(da+1,1),0.7,data(da+1,3)],'facecolor','w')
    rectangle('position',[da+2-0.35,data(da+1,1)+data(da+1,3),0.7,data(da+1,2)],'facecolor','k')
    %     bar(data,'stacked')
    if tt==3
        rectangle('position',[da+4,37,1,6],'facecolor','c')
        %         rectangle('position',[da+4,47,1,6],'facecolor','m')
        rectangle('position',[da+4,57,1,6],'facecolor','k')
        text(da+5.2,40,'on')
        %          text(da+5.2,50,'onoff')
        text(da+5.2,60,'off')
        axis([0 da+7 0 100])
    else
        axis([0 da+3 0 100])
    end
    dates=[dates ' ' 'all'];
    set(gca,'xtick',1:length(dates))
    set(gca,'xticklabel',dates,'fontsize',7)
    
    title(titles{tt})
    ylabel('% of cells with LF and/or flash response')
end
figure(1)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'polarity_flash.bmp'))
saveas(gcf,fullfile(savepath,'polarity_flash.fig'))
figure(2)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'polarity_LF.bmp'))
saveas(gcf,fullfile(savepath,'polarity_LF.fig'))
close all
%% ??datewise bootstrapping for 6vel
clear
close all
% superpath='S:\data\Katja\allstuff';
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

id=17;
border=3;
types='wph';
cols='bgr';
titles={'mouse','pig','human'};


figure
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    goods=[];
    for q=1:length(info)
        if goodones(q)==1 && ~isempty(info{q,id})
            if info{q,id}>0
                goods=[goods;q];
            end
        end
    end
    datas=cell2mat(info(goods,id));
    dates=info(goods,1);
    dats=unique(dates);
    
    prevlength=[0 0 0];
    for d=1:length(dats) %go through dates of species 1
        curr=[];
        for i=1:length(dates)
            if ~isempty(regexp(dates{i},dats{d}))
                curr=[curr;i];
            end
        end
        
        data=datas(curr);
        
        for t2=tt+1:length(types) %go through other species
            Dtype=types(t2);
            load(fullfile(superpath,[Dtype,'_info']))
            goods=[];
            for q=1:length(info)
                if goodones(q)==1 && ~isempty(info{q,id})
                    if info{q,id}>0
                        goods=[goods;q];
                    end
                end
            end
            datas2=cell2mat(info(goods,id));
            dates2=info(goods,1);
            dats2=unique(dates2);
            
            for d2=1:length(dats2) %go through dates of other species
                curr2=[];
                for i=1:length(dates2)
                    if ~isempty(regexp(dates2{i},dats2{d2}))
                        curr2=[curr2;i];
                    end
                end
                
                data2=datas2(curr2);
                
                me1=median(data);
                me2=median(data2);
                if length(curr2)>=border && length(data)>=border
                    p=ranksum(data,data2);
                    
                    collstats{prevlength(t2)+d2,tt+t2-2}=p;
                end
                medians{d,tt}=me1;
                medians{d2,t2}=me2;
                
            end
            prevlength(t2)=prevlength(t2)+d2;
        end
    end
end

col=[0 0.6 0.6;0.6 0 0.6;0.6 0.6 0];
for i=1:3
    data=cell2mat(collstats(:,i));
    plot(i,data,'o','color',col(i,:))
    hold on
    plot(i,median(data),'s','markerfacecolor',col(i,:),'markeredgecolor',col(i,:),'markersize',12)
    plot([i i],[median(data)-std(data) median(data)+std(data)],'-','color',col(i,:))
end
plot([0 4],[0.05 0.05],'-k')
axis([0 4 -inf inf])
%% "normal" bootstrapping for 6vel
clear
% close all
% superpath='S:\data\Katja\allstuff';
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

id=16;
samplesize=10;
numsamples=100;
border=5;
types='wph';
cols='bgr';
titles={'mouse','pig','human'};


figure
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    goods=[];
    for q=1:length(info)
        if goodones(q)==1 && ~isempty(info{q,id})
            if info{q,id}>0
                goods=[goods;q];
            end
        end
    end
    datas=cell2mat(info(goods,id));
    
    prevlength=[0 0 0];
    for d=1:numsamples
        now=randi(length(datas),samplesize,1);
        now=unique(now);
        
        data=datas(now);
        
        
        for t2=tt+1:length(types) %go through other species
            Dtype=types(t2);
            load(fullfile(superpath,[Dtype,'_info']))
            goods=[];
            for q=1:length(info)
                if goodones(q)==1 && ~isempty(info{q,id})
                    if info{q,id}>0
                        goods=[goods;q];
                    end
                end
            end
            datas2=cell2mat(info(goods,id));
            
            for d2=1:numsamples %go through dates of other species
                now=randi(length(datas2),samplesize,1);
                now=unique(now);
                
                data2=datas2(now);
                
                me1=median(data);
                me2=median(data2);
                if length(now)>=border && length(data)>=border
                    p=ranksum(data,data2);
                    
                    collstats{prevlength(t2)+d2,tt+t2-2}=p;
                end
                medians{d,tt}=me1;
                medians{d2,t2}=me2;
                
            end
            prevlength(t2)=prevlength(t2)+d2;
        end
    end
end

col=[0 0.6 0.6;0.6 0 0.6;0.6 0.6 0];
for i=1:3
    data=cell2mat(collstats(:,i));
    plot(i,data,'o','color',col(i,:))
    hold on
    plot(i,median(data),'s','markerfacecolor',col(i,:),'markeredgecolor',col(i,:),'markersize',12)
    plot([i i],[median(data)-std(data) median(data)+std(data)],'-','color',col(i,:))
end
plot([0 4],[0.05 0.05],'-k')
axis([0 4 -inf inf])
%% 10e) recalculate DG values: sum of first 4 harmonics
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)


DGpeaks=cell(length(info),1);
togo=find(goodDG==1);
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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    thiscell=[];
    fouriers=[];
    rates=[];
    clear f2take
    for g=1:length(DGlist)
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            
            stimFreq=freqs(g);
            
            if max(meanrate)>0
                meanF=mean(meanrate);
                meanFR=round(meanF*10)/10;
                % fourier
                NFFT=2^nextpow2(numel(meanrate));
                Y=fft(meanrate,NFFT)/numel(meanrate);
                f=1000/2*linspace(0,1,NFFT/2+1);
                if ~exist('f2take','var') && ~isempty(find(f>0))
                    f2take=f;
                end
                c=2*abs(Y(1:NFFT/2+1));
                
                
                ampPeak=max(c(4:end));
                ampPeakR=round(ampPeak);
                posF=find(c==ampPeak);
                posF=posF(1);
                freq=f(posF);
                
                collpeaks=0;
                testfreq=stimFreq:stimFreq:40;
                testfreq=testfreq(1:4);
                for te=1:length(testfreq)
                    stimFreqN=testfreq(te);
                    findfreq1=find(f<stimFreqN);
                    findfreq1=findfreq1(end)-1;
                    if findfreq1<=0
                        findfreq1=4;
                    end
                    
                    findfreq2=max(c(findfreq1:findfreq1+5)); %recorded amp at stimfreq
                    if te==1
                        pos1=find(c==findfreq2);
                        pos1=f(pos1);
                    end
                    collpeaks=collpeaks+findfreq2;
                end
            else
                findfreq2=0;
                pos1=0;
                
            end
            thiscell=[thiscell;collpeaks pos1]; %max at freq, corresponding freq (exactly)
            size(rates,1)
        end
    end
    DGpeaks{currID}=thiscell;
end
save(fullfile(superpath,[Dtype,'_DGpeaks4']),'DGpeaks')
%% 10f) recalculate DG values: STD of firing rate, baseline subtracted
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)


DGpeaks=cell(length(info),1);
togo=find(goodDG==1);
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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    thiscell=[];
    fouriers=[];
    rates=[];
    clear f2take
    for g=1:length(DGlist)
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            
            stimFreq=freqs(g);
            
            if max(meanrate)>0
                backgr=std(meanrate(500:1000));
                stim=std(meanrate(2000:14000));
                findfreq2=stim-backgr;
                pos1=0;
            else
                findfreq2=0;
                pos1=0;
                
            end
            thiscell=[thiscell;findfreq2 pos1]; %max at freq, corresponding freq (exactly)
            size(rates,1)
        end
    end
    DGpeaks{currID}=thiscell;
end
save(fullfile(superpath,[Dtype,'_DGpeaksSTD']),'DGpeaks')
%% --10g) PLOT spatiotemporal again for sum of 4 first harmonics
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaks4']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            allDataTemp{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    
    
    figure(1)
    data=sumSpat;
    subplot(2,2,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,2,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,2,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,2,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,2,s)
    axis([0 7 -inf inf])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq_sum4harm.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq_sum4harm.fig'))
close
figure(2)
for s=1:4
    subplot(2,2,s)
    axis([0 5 -inf inf])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq_sum4harm.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq_sum4harm.fig'))
close
%% --10h) PLOT spatiotemporal again for sum of STD
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaksSTD']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            allDataTemp{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    
    figure(1)
    data=sumSpat;
    subplot(2,2,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,2,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,2,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,2,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,2,s)
    axis([0 7 -0.5 2])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq_STD.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq_STD.fig'))
close
figure(2)
for s=1:4
    subplot(2,2,s)
    axis([0 5 -0.5 2])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq_STD.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq_STD.fig'))
close
%% --10i) PLOT spatiotemporal again for sum of 4 first harmonics: not normalized
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaks4']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            allDataTemp{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    figure(1)
    data=sumSpat;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,3,s)
    axis([0 7 -inf inf])
end
subplot(2,3,5)
axis([0 7 0 1])

set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq_sum4harm_nonNorm.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq_sum4harm_nonNorm.fig'))
close
figure(2)
for s=1:4
    subplot(2,3,s)
    axis([0 5 -inf inf])
end
subplot(2,3,5)
axis([0 5 0 1])

set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq_sum4harm_nonNorm.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq_sum4harm_nonNorm.fig'))
close
%% --10j) PLOT spatiotemporal again for sum of STD: not normalized
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaksSTD']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            allDataTemp{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    
    figure(1)
    data=sumSpat;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,3,s)
    axis([0 7 -inf inf])
end
subplot(2,3,5)
axis([0 7 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq_STD_nonNorm.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq_STD_nonNorm.fig'))
close
figure(2)
for s=1:4
    subplot(2,3,s)
    axis([0 5 -inf inf])
end
subplot(2,3,5)
axis([0 5 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq_STD_nonNorm.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq_STD_nonNorm.fig'))
close
%% --10k) PLOT spatiotemporal again for "normal calculation" (first harm): not normalized
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaks']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            allDataTemp{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    
    figure(1)
    data=sumSpat;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,3,s)
    axis([0 7 -inf inf])
end
subplot(2,3,5)
axis([0 7 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq_firstharm_nonNorm.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq_firstharm_nonNorm.fig'))
close
figure(2)
for s=1:4
    subplot(2,3,s)
    axis([0 5 -inf inf])
end
subplot(2,3,5)
axis([0 5 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq_firstharm_nonNorm.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq_firstharm_nonNorm.fig'))
close
%% --10l) PLOT spatiotemporal again for "normal calculation" (first harm), speed controlled: not normalized
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');

tempmat=[4 5 6; 3 4 5; 2 3 4; 1 2 3];
spatmat=[4 3 2 1];


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaks']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    spats=spats(2:end-1);
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        %spatial
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            tmp=tmp(spatmat(s));
            ids=tmp; %all temp freq to consider
            su=sum(curr(ids));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        %temporal
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            ids=tmp; %all spat freq to consider
            ids=ids(tempmat(t,:));
            allDataTemp{t,1}=curr(ids);
            su=sum(curr(ids));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    
    figure(1)
    data=sumSpat;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,3,s)
    axis([0 7 -inf inf])
end
subplot(2,3,5)
axis([0 7 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq_firstharm_nonNorm_speedCorr.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq_firstharm_nonNorm_speedCorr.fig'))
close
figure(2)
for s=1:4
    subplot(2,3,s)
    axis([0 5 -inf inf])
end
subplot(2,3,5)
axis([0 5 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq_firstharm_nonNorm_speedCorr.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq_firstharm_nonNorm_speedCorr.fig'))
close
%% --10m) PLOT spatiotemporal again for sum 4 harm, speed controlled: not normalized
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');

tempmat=[4 5 6; 3 4 5; 2 3 4; 1 2 3];
spatmat=[4 3 2 1];


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaks4']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    spats=spats(2:end-1);
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        %spatial
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            tmp=tmp(spatmat(s));
            ids=tmp; %all temp freq to consider
            su=sum(curr(ids));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        %temporal
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            ids=tmp; %all spat freq to consider
            ids=ids(tempmat(t,:));
            allDataTemp{t,1}=curr(ids);
            su=sum(curr(ids));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    
    figure(1)
    data=sumSpat;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,3,s)
    axis([0 7 -inf inf])
end
subplot(2,3,5)
axis([0 7 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq_sum4harm_nonNorm_speedCorr.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq_sum4harm_nonNorm_speedCorr.fig'))
close
figure(2)
for s=1:4
    subplot(2,3,s)
    axis([0 5 -inf inf])
end
subplot(2,3,5)
axis([0 5 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq_sum4harm_nonNorm_speedCorr.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq_sum4harm_nonNorm_speedCorr.fig'))
close
%% --10n) PLOT spatiotemporal again for STD, speed controlled: not normalized
clear
close all
superpath='S:\data\Katja\allstuff';
savepath=fullfile(superpath,'overview_pics');

tempmat=[4 5 6; 3 4 5; 2 3 4; 1 2 3];
spatmat=[4 3 2 1];


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure(1)
figure(2)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_DGpeaksSTD']))
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGpeaks{m})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    spats=spats(2:end-1);
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    sumSpat=[]; sumTemp=[]; sumSpe=[];
    
    togo=find(goodDG==1);
    for i=1:length(togo)
        currID=togo(i);
        curr=DGpeaks{currID}(:,1);
        
        %spatial
        allsum=[];
        for s=1:length(spats)
            tmp=find(DGstim(:,2)==spats(s));
            tmp=tmp(spatmat(s));
            ids=tmp; %all temp freq to consider
            su=sum(curr(ids));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpat=[sumSpat;allsum];
        
        %temporal
        clear allDataTemp
        allsum=[];
        for t=1:length(temps)
            tmp=find(DGstim(:,1)==temps(t));
            ids=tmp; %all spat freq to consider
            ids=ids(tempmat(t,:));
            allDataTemp{t,1}=curr(ids);
            su=sum(curr(ids));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumTemp=[sumTemp;allsum];
        
        allsum=[];
        clear allDataSpe
        for t=1:length(spe)
            tmp=find(spes==spe(t));
            allDataSpe{t,1}=curr(tmp);
            su=sum(curr(tmp));
            allsum=[allsum su];
        end
        %         allsum=allsum/max(allsum);
        sumSpe=[sumSpe;allsum];
    end
    
    
    
    figure(1)
    data=sumSpat;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    set(gca,'xtick',1:6)
    set(gca,'xticklabel',spats)
    xlabel('spatial freq. (um)')
    
    
    figure(2)
    data=sumTemp;
    subplot(2,3,tt)
    plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    title([titles{tt},' (total: ',int2str(size(data,1)),')'])
    subplot(2,3,4)
    plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    for p=1:size(data,2)
        plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
        plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
    end
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
    subplot(2,3,5)
    meddata=median(data);
    plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
    hold on
    
    
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',temps)
    xlabel('temporal freq. (Hz)')
    
    
end

figure(1)
for s=1:4
    subplot(2,3,s)
    axis([0 7 -inf inf])
end
subplot(2,3,5)
axis([0 7 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'spatial_freq_STD_nonNorm_speedCorr.bmp'))
saveas(gcf,fullfile(savepath,'spatial_freq_STD_nonNorm_speedCorr.fig'))
close
figure(2)
for s=1:4
    subplot(2,3,s)
    axis([0 5 -inf inf])
end
subplot(2,3,5)
axis([0 5 0 1])
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'temporal_freq_STD_nonNorm_speedCorr.bmp'))
saveas(gcf,fullfile(savepath,'temporal_freq_STD_nonNorm_speedCorr.fig'))
close
%% 10o) PLOT spatiotemporal for all three conditions, not normalized, max (instead of sum)
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

conditions={'DGpeaks','DGpeaks4','DGpeaksSTD'};
titleadd={'firstharm','sum4harm','STD'};


for co=1:length(conditions)
    types='wph';
    cols='bgr';
    titles={'mouse','pig','human'};
    
    figure(1)
    figure(2)
    for tt=1:length(types)
        Dtype=types(tt);
        load(fullfile(superpath,[Dtype,'_info']))
        load(fullfile(superpath,[Dtype,'_goodDG']))
        load(fullfile(superpath,[Dtype,'_',conditions{co}]))
        
        found=0;m=0;
        while found==0
            m=m+1;
            if ~isempty(DGpeaks{m})
                found=1;
            end
        end
        DGstim=[stimInfo{m,1} stimInfo{m,2}];
        spats=unique(DGstim(:,2));
        temps=unique(DGstim(:,1));
        spes=[];
        for ss=1:length(DGstim)
            noww=DGstim(ss,1)*DGstim(ss,2);
            spes=[spes;noww];
        end
        spes=spes/1000;
        spe=unique(spes);
        
        sumSpat=[]; sumTemp=[]; sumSpe=[];
        
        togo=find(goodDG==1);
        for i=1:length(togo)
            currID=togo(i);
            curr=DGpeaks{currID}(:,1);
            
            allsum=[];
            for s=1:length(spats)
                tmp=find(DGstim(:,2)==spats(s));
                su=max((curr(tmp)));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpat=[sumSpat;allsum];
            
            clear allDataTemp
            allsum=[];
            for t=1:length(temps)
                tmp=find(DGstim(:,1)==temps(t));
                allDataTemp{t,1}=curr(tmp);
                su=max(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumTemp=[sumTemp;allsum];
            
            allsum=[];
            clear allDataSpe
            for t=1:length(spe)
                tmp=find(spes==spe(t));
                allDataSpe{t,1}=curr(tmp);
                su=max(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpe=[sumSpe;allsum];
        end
        
        
        
        figure(1)
        data=sumSpat;
        subplot(2,3,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,3,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        subplot(2,3,5)
        meddata=median(data);
        plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        
        figure(2)
        data=sumTemp;
        subplot(2,3,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,3,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        
        subplot(2,3,5)
        meddata=median(data);
        plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        save(fullfile(superpath,[types(tt),'_',conditions{co},'_max_all']),'sumSpat','sumTemp')
        
    end
    
    figure(1)
    for s=1:4
        subplot(2,3,s)
        axis([0 7 -inf inf])
    end
    subplot(2,3,5)
    axis([0 7 0 1])
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['spatial_freq_',titleadd{co},'_max.bmp']))
    saveas(gcf,fullfile(savepath,['spatial_freq_',titleadd{co},'_max.fig']))
    close
    figure(2)
    for s=1:4
        subplot(2,3,s)
        axis([0 5 -inf inf])
    end
    subplot(2,3,5)
    axis([0 5 0 1])
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['temporal_freq_',titleadd{co},'_max.bmp']))
    saveas(gcf,fullfile(savepath,['temporal_freq_',titleadd{co},'_max.fig']))
    close
    
end
%% 10p) PLOT spatiotemporal for all three conditions, speed controlled, not normalized, max (instead of sum)
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

tempmat=[4 5 6; 3 4 5; 2 3 4; 1 2 3];
spatmat=[4 3 2 1];


conditions={'DGpeaks','DGpeaks4','DGpeaksSTD'};
titleadd={'firstharm','sum4harm','STD'};


for co=1:length(conditions)
    types='wph';
    cols='bgr';
    titles={'mouse','pig','human'};
    
    figure(1)
    figure(2)
    for tt=1:length(types)
        Dtype=types(tt);
        load(fullfile(superpath,[Dtype,'_info']))
        load(fullfile(superpath,[Dtype,'_goodDG']))
        load(fullfile(superpath,[Dtype,'_',conditions{co}]))
        
        found=0;m=0;
        while found==0
            m=m+1;
            if ~isempty(DGpeaks{m})
                found=1;
            end
        end
        DGstim=[stimInfo{m,1} stimInfo{m,2}];
        spats=unique(DGstim(:,2));
        spats=spats(2:end-1);
        temps=unique(DGstim(:,1));
        spes=[];
        for ss=1:length(DGstim)
            noww=DGstim(ss,1)*DGstim(ss,2);
            spes=[spes;noww];
        end
        spes=spes/1000;
        spe=unique(spes);
        
        sumSpat=[]; sumTemp=[]; sumSpe=[];
        
        togo=find(goodDG==1);
        for i=1:length(togo)
            currID=togo(i);
            curr=DGpeaks{currID}(:,1);
            
            %spatial
            allsum=[];
            for s=1:length(spats)
                tmp=find(DGstim(:,2)==spats(s));
                tmp=tmp(spatmat(s));
                ids=tmp; %all temp freq to consider
                su=max(curr(ids));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpat=[sumSpat;allsum];
            
            %temporal
            clear allDataTemp
            allsum=[];
            for t=1:length(temps)
                tmp=find(DGstim(:,1)==temps(t));
                ids=tmp; %all spat freq to consider
                ids=ids(tempmat(t,:));
                allDataTemp{t,1}=curr(ids);
                su=max(curr(ids));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumTemp=[sumTemp;allsum];
            
            allsum=[];
            clear allDataSpe
            for t=1:length(spe)
                tmp=find(spes==spe(t));
                allDataSpe{t,1}=curr(tmp);
                su=max(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpe=[sumSpe;allsum];
        end
        
        
        
        figure(1)
        data=sumSpat;
        subplot(2,3,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,3,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        subplot(2,3,5)
        meddata=median(data);
        plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        
        figure(2)
        data=sumTemp;
        subplot(2,3,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,3,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        
        subplot(2,3,5)
        meddata=median(data);
        plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        save(fullfile(superpath,[types(tt),'_',conditions{co},'_max_speedC']),'sumSpat','sumTemp')
        
        
    end
    
    figure(1)
    for s=1:4
        subplot(2,3,s)
        axis([0 7 -inf inf])
    end
    subplot(2,3,5)
    axis([0 7 0 1])
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['spatial_freq_',titleadd{co},'_speedCorr_max.bmp']))
    saveas(gcf,fullfile(savepath,['spatial_freq_',titleadd{co},'_speedCorr_max.fig']))
    
    close
    figure(2)
    for s=1:4
        subplot(2,3,s)
        axis([0 5 -inf inf])
    end
    subplot(2,3,5)
    axis([0 5 0 1])
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['temporal_freq_',titleadd{co},'_speedCorr_max.bmp']))
    saveas(gcf,fullfile(savepath,['temporal_freq_',titleadd{co},'_speedCorr_max.fig']))
    
    close
    
    
end
%% 10q) PLOT spatiotemporal for all three conditions, not normalized, sum (replaces 10i, 10j, k)
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

conditions={'DGpeaks','DGpeaks4','DGpeaksSTD'};
titleadd={'firstharm','sum4harm','STD'};


for co=1:length(conditions)
    types='wph';
    cols='bgr';
    titles={'mouse','pig','human'};
    
    figure(1)
    figure(2)
    for tt=1:length(types)
        Dtype=types(tt);
        load(fullfile(superpath,[Dtype,'_info']))
        load(fullfile(superpath,[Dtype,'_goodDG']))
        load(fullfile(superpath,[Dtype,'_',conditions{co}]))
        
        found=0;m=0;
        while found==0
            m=m+1;
            if ~isempty(DGpeaks{m})
                found=1;
            end
        end
        DGstim=[stimInfo{m,1} stimInfo{m,2}];
        spats=unique(DGstim(:,2));
        temps=unique(DGstim(:,1));
        spes=[];
        for ss=1:length(DGstim)
            noww=DGstim(ss,1)*DGstim(ss,2);
            spes=[spes;noww];
        end
        spes=spes/1000;
        spe=unique(spes);
        
        sumSpat=[]; sumTemp=[]; sumSpe=[];
        
        togo=find(goodDG==1);
        for i=1:length(togo)
            currID=togo(i);
            curr=DGpeaks{currID}(:,1);
            
            allsum=[];
            for s=1:length(spats)
                tmp=find(DGstim(:,2)==spats(s));
                su=sum((curr(tmp)));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpat=[sumSpat;allsum];
            
            clear allDataTemp
            allsum=[];
            for t=1:length(temps)
                tmp=find(DGstim(:,1)==temps(t));
                allDataTemp{t,1}=curr(tmp);
                su=sum(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumTemp=[sumTemp;allsum];
            
            allsum=[];
            clear allDataSpe
            for t=1:length(spe)
                tmp=find(spes==spe(t));
                allDataSpe{t,1}=curr(tmp);
                su=sum(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpe=[sumSpe;allsum];
        end
        
        
        
        figure(1)
        data=sumSpat;
        subplot(2,3,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,3,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        subplot(2,3,5)
        meddata=median(data);
        plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        
        figure(2)
        data=sumTemp;
        subplot(2,3,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,3,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        
        subplot(2,3,5)
        meddata=median(data);
        plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        save(fullfile(superpath,[types(tt),'_',conditions{co},'_sum_all']),'sumSpat','sumTemp')
        
    end
    
    figure(1)
    for s=1:4
        subplot(2,3,s)
        axis([0 7 -inf inf])
    end
    subplot(2,3,5)
    axis([0 7 0 1])
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['spatial_freq_',titleadd{co},'_sum.bmp']))
    saveas(gcf,fullfile(savepath,['spatial_freq_',titleadd{co},'_sum.fig']))
    close
    figure(2)
    for s=1:4
        subplot(2,3,s)
        axis([0 5 -inf inf])
    end
    subplot(2,3,5)
    axis([0 5 0 1])
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['temporal_freq_',titleadd{co},'_sum.bmp']))
    saveas(gcf,fullfile(savepath,['temporal_freq_',titleadd{co},'_sum.fig']))
    close
    
end
%% 10r) PLOT spatiotemporal for all three conditions, speed controlled, not normalized, sum, (replaces 10l, 10m, 10n)
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

tempmat=[4 5 6; 3 4 5; 2 3 4; 1 2 3];
spatmat=[4 3 2 1];


conditions={'DGpeaks','DGpeaks4','DGpeaksSTD'};
titleadd={'firstharm','sum4harm','STD'};


for co=1:length(conditions)
    types='wph';
    cols='bgr';
    titles={'mouse','pig','human'};
    
    figure(1)
    figure(2)
    for tt=1:length(types)
        Dtype=types(tt);
        load(fullfile(superpath,[Dtype,'_info']))
        load(fullfile(superpath,[Dtype,'_goodDG']))
        load(fullfile(superpath,[Dtype,'_',conditions{co}]))
        
        found=0;m=0;
        while found==0
            m=m+1;
            if ~isempty(DGpeaks{m})
                found=1;
            end
        end
        DGstim=[stimInfo{m,1} stimInfo{m,2}];
        spats=unique(DGstim(:,2));
        spats=spats(2:end-1);
        temps=unique(DGstim(:,1));
        spes=[];
        for ss=1:length(DGstim)
            noww=DGstim(ss,1)*DGstim(ss,2);
            spes=[spes;noww];
        end
        spes=spes/1000;
        spe=unique(spes);
        
        sumSpat=[]; sumTemp=[]; sumSpe=[];
        
        togo=find(goodDG==1);
        for i=1:length(togo)
            currID=togo(i);
            curr=DGpeaks{currID}(:,1);
            
            %spatial
            allsum=[];
            for s=1:length(spats)
                tmp=find(DGstim(:,2)==spats(s));
                tmp=tmp(spatmat(s));
                ids=tmp; %all temp freq to consider
                su=sum(curr(ids));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpat=[sumSpat;allsum];
            
            %temporal
            clear allDataTemp
            allsum=[];
            for t=1:length(temps)
                tmp=find(DGstim(:,1)==temps(t));
                ids=tmp; %all spat freq to consider
                ids=ids(tempmat(t,:));
                allDataTemp{t,1}=curr(ids);
                su=sum(curr(ids));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumTemp=[sumTemp;allsum];
            
            allsum=[];
            clear allDataSpe
            for t=1:length(spe)
                tmp=find(spes==spe(t));
                allDataSpe{t,1}=curr(tmp);
                su=sum(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpe=[sumSpe;allsum];
        end
        
        
        
        figure(1)
        data=sumSpat;
        subplot(2,3,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,3,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        subplot(2,3,5)
        meddata=median(data);
        plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        
        figure(2)
        data=sumTemp;
        subplot(2,3,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,3,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        
        subplot(2,3,5)
        meddata=median(data);
        plot(meddata/max(meddata),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        save(fullfile(superpath,[types(tt),'_',conditions{co},'_sum_speedC']),'sumSpat','sumTemp')
        
    end
    
    figure(1)
    for s=1:4
        subplot(2,3,s)
        axis([0 7 -inf inf])
    end
    subplot(2,3,5)
    axis([0 7 0 1])
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['spatial_freq_',titleadd{co},'_speedCorr_sum.bmp']))
    saveas(gcf,fullfile(savepath,['spatial_freq_',titleadd{co},'_speedCorr_sum.fig']))
    
    close
    figure(2)
    for s=1:4
        subplot(2,3,s)
        axis([0 5 -inf inf])
    end
    subplot(2,3,5)
    axis([0 5 0 1])
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['temporal_freq_',titleadd{co},'_speedCorr_sum.bmp']))
    saveas(gcf,fullfile(savepath,['temporal_freq_',titleadd{co},'_speedCorr_sum.fig']))
    
    close
    
end
%% 10s) PLOT spatiotemporal for all three conditions, normalized, sum
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

conditions={'DGpeaks','DGpeaks4','DGpeaksSTD'};
titleadd={'firstharm','sum4harm','STD'};


for co=1:length(conditions)
    types='wph';
    cols='bgr';
    titles={'mouse','pig','human'};
    
    figure(1)
    figure(2)
    for tt=1:length(types)
        Dtype=types(tt);
        load(fullfile(superpath,[Dtype,'_info']))
        load(fullfile(superpath,[Dtype,'_goodDG']))
        load(fullfile(superpath,[Dtype,'_',conditions{co}]))
        
        found=0;m=0;
        while found==0
            m=m+1;
            if ~isempty(DGpeaks{m})
                found=1;
            end
        end
        DGstim=[stimInfo{m,1} stimInfo{m,2}];
        spats=unique(DGstim(:,2));
        temps=unique(DGstim(:,1));
        spes=[];
        for ss=1:length(DGstim)
            noww=DGstim(ss,1)*DGstim(ss,2);
            spes=[spes;noww];
        end
        spes=spes/1000;
        spe=unique(spes);
        
        sumSpat=[]; sumTemp=[]; sumSpe=[];
        
        togo=find(goodDG==1);
        for i=1:length(togo)
            currID=togo(i);
            curr=DGpeaks{currID}(:,1);
            curr=curr/max(curr);
            
            allsum=[];
            for s=1:length(spats)
                tmp=find(DGstim(:,2)==spats(s));
                su=sum((curr(tmp)));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpat=[sumSpat;allsum];
            
            clear allDataTemp
            allsum=[];
            for t=1:length(temps)
                tmp=find(DGstim(:,1)==temps(t));
                allDataTemp{t,1}=curr(tmp);
                su=sum(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumTemp=[sumTemp;allsum];
            
            allsum=[];
            clear allDataSpe
            for t=1:length(spe)
                tmp=find(spes==spe(t));
                allDataSpe{t,1}=curr(tmp);
                su=sum(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpe=[sumSpe;allsum];
        end
        
        
        
        figure(1)
        data=sumSpat;
        subplot(2,2,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,2,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        
        
        figure(2)
        data=sumTemp;
        subplot(2,2,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,2,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        
    end
    
    figure(1)
    for s=1:4
        subplot(2,2,s)
        axis([0 7 -inf inf])
    end
    
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['spatial_norm_',titleadd{co},'_sum.bmp']))
    saveas(gcf,fullfile(savepath,['spatial_norm_',titleadd{co},'_sum.fig']))
    close
    figure(2)
    for s=1:4
        subplot(2,2,s)
        axis([0 5 -inf inf])
    end
    
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['temporal_norm_',titleadd{co},'_sum.bmp']))
    saveas(gcf,fullfile(savepath,['temporal_norm_',titleadd{co},'_sum.fig']))
    close
    
end
%% 10t) PLOT spatiotemporal for all three conditions, normalized, max
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

conditions={'DGpeaks','DGpeaks4','DGpeaksSTD'};
titleadd={'firstharm','sum4harm','STD'};


for co=1:length(conditions)
    types='wph';
    cols='bgr';
    titles={'mouse','pig','human'};
    
    figure(1)
    figure(2)
    for tt=1:length(types)
        Dtype=types(tt);
        load(fullfile(superpath,[Dtype,'_info']))
        load(fullfile(superpath,[Dtype,'_goodDG']))
        load(fullfile(superpath,[Dtype,'_',conditions{co}]))
        
        found=0;m=0;
        while found==0
            m=m+1;
            if ~isempty(DGpeaks{m})
                found=1;
            end
        end
        DGstim=[stimInfo{m,1} stimInfo{m,2}];
        spats=unique(DGstim(:,2));
        temps=unique(DGstim(:,1));
        spes=[];
        for ss=1:length(DGstim)
            noww=DGstim(ss,1)*DGstim(ss,2);
            spes=[spes;noww];
        end
        spes=spes/1000;
        spe=unique(spes);
        
        sumSpat=[]; sumTemp=[]; sumSpe=[];
        
        togo=find(goodDG==1);
        for i=1:length(togo)
            currID=togo(i);
            curr=DGpeaks{currID}(:,1);
            curr=curr/max(curr);
            
            allsum=[];
            for s=1:length(spats)
                tmp=find(DGstim(:,2)==spats(s));
                su=max((curr(tmp)));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpat=[sumSpat;allsum];
            
            clear allDataTemp
            allsum=[];
            for t=1:length(temps)
                tmp=find(DGstim(:,1)==temps(t));
                allDataTemp{t,1}=curr(tmp);
                su=max(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumTemp=[sumTemp;allsum];
            
            allsum=[];
            clear allDataSpe
            for t=1:length(spe)
                tmp=find(spes==spe(t));
                allDataSpe{t,1}=curr(tmp);
                su=max(curr(tmp));
                allsum=[allsum su];
            end
            %         allsum=allsum/max(allsum);
            sumSpe=[sumSpe;allsum];
        end
        
        
        
        figure(1)
        data=sumSpat;
        subplot(2,2,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %         fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,2,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spats)
        xlabel('spatial freq. (um)')
        
        
        
        figure(2)
        data=sumTemp;
        subplot(2,2,tt)
        plot(median(data),'o','markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        %     fill([1:length(median(data)) length(median(data)):-1:1],[median(data)-std(data) fliplr(median(data)+std(data))],cols(tt),'facealpha',0.5)
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        title([titles{tt},' (total: ',int2str(size(data,1)),')'])
        subplot(2,2,4)
        plot(median(data),'-o','color',cols(tt),'markerfacecolor',cols(tt),'markeredgecolor',cols(tt))
        hold on
        for p=1:size(data,2)
            plot([p p],[median(data(:,p))-1.96*std(data(:,p))/sqrt(length(data(:,p))) median(data(:,p))+1.96*std(data(:,p))/sqrt(length(data(:,p)))],'-','color',cols(tt),'linewidth',3)
            plot([p p],[median(data(:,p))-std(data(:,p)) median(data(:,p))+std(data(:,p))],'-','color',cols(tt),'linewidth',1)
        end
        
        
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',temps)
        xlabel('temporal freq. (Hz)')
        
        
    end
    
    figure(1)
    for s=1:4
        subplot(2,2,s)
        axis([0 7 -inf inf])
    end
    
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['spatial_norm_',titleadd{co},'_max.bmp']))
    saveas(gcf,fullfile(savepath,['spatial_norm_',titleadd{co},'_max.fig']))
    close
    figure(2)
    for s=1:4
        subplot(2,2,s)
        axis([0 5 -inf inf])
    end
    
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,['temporal_norm_',titleadd{co},'_max.bmp']))
    saveas(gcf,fullfile(savepath,['temporal_norm_',titleadd{co},'_max.fig']))
    close
    
end
%% to be worked on
type='m';
superpath='S:\data\Katja\allstuff';

load(fullfile(superpath,[type,'_DGpeaks_sum_all']))

sumspeedS=sumSpat;
sumspeedT=sumTemp;
%%
load(fullfile(superpath,[type,'_DGpeaks_max_all']))

sumallS=sumSpat;
sumallT=sumTemp;
% sumallS=sumallS(:,2:5);
%%
figure

for i=1:6
    subplot(2,6,i)
    plot(sumspeedS(:,i),sumallS(:,i),'.')
end
for i=1:4
    subplot(2,6,i+6)
    plot(sumspeedT(:,i),sumallT(:,i),'.')
end
%% compare LF and flash
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');
flash=8; LF=9;

types='wph';
cols='bgr';
titles={'mouse','pig','human'};
list=[1 -1 0];
li={'LF on','LF off','LF none'};
list2=[1 -1 2 0];
li2={'Flash on','Flash off','Flash on-off','Flash none'};

figure(1)
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    togo=find(cell2mat(info(:,3))>0);
    flashdata=info(togo,flash);
    for f=1:length(flashdata)
        if isempty(flashdata{f})
            flashdata{f}=0;
        end
    end %0=none, 1=on, -1=off, 2=onoff
    flashdata=cell2mat(flashdata);
    LFdata=info(togo,LF);
    for f=1:length(LFdata)
        if isempty(LFdata{f})
            LFdata{f}=0;
        end
    end
    LFdata=cell2mat(LFdata);
    
    coll=zeros(3,4);
    for l=1:length(list)
        curr=find(LFdata==list(l));
        for o=1:length(list2)
            curr2=find(flashdata(curr)==list2(o));
            coll(l,o)=length(curr2);
        end
    end
    subplot(2,2,tt)
    bar(coll)
    set(gca,'xticklabel',li)
    legend(li2)
    title(titles{tt})
    
    if tt==3
        axis([0.5 3.5 0 100])
    end
    
    
    %     dates=info(togo,1);
    %     udates=unique(dates);
    %     coll=zeros(3,4,length(udates));
    %     keepLF=LFdata;
    %     keepFlash=flashdata;
    %     for u=1:length(udates)
    %         currdates=[];
    %         for c=1:length(dates)
    %             if ~isempty(regexp(dates{c},udates{u}))
    %                 currdates=[currdates c];
    %             end
    %         end
    %         LFdata=keepLF(currdates);
    %         flashdata=keepFlash(currdates);
    %         for l=1:length(list)
    %             curr=find(LFdata==list(l));
    %             for o=1:length(list2)
    %                 curr2=find(flashdata(curr)==list2(o));
    %                 coll(l,o,u)=length(curr2);
    %             end
    %         end
    %     end
    %     subplot(2,3,tt+3)
    %     bar(coll)
    %     set(gca,'xticklabel',li)
    %     legend(li2)
    %     title(titles{tt})
    
    
end
%% 12a) more on full-field
% col7) #quick files
% col8) pol from flash: 0=none, -1=off, 1=on, 2=onoff

load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
if ~exist(fullfile(superpath,'FFFlash'),'dir')
    mkdir(fullfile(superpath,'FFFlash'))
end
togo=find(goodones>0);
tmp=cell2mat(info(togo,8));
tmp2=find(tmp~=0);
togo=togo(tmp2);
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
    
    quicks=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'FFFlash'))
            ftype=1;
            quicks=[quicks;1];
        elseif ~isempty(regexp(hekalist{f},'quick'))
            ftype=2;
            quicks=[quicks;1];
        else
            quicks=[quicks;0];
        end
    end
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(quicks==1);
    if Dtype=='m'
        goodfiles=goodfiles(11:end);
    else
        goodfiles=goodfiles(5:end);
    end
    
    if ftype==1 %FFFlash
        first=goodfiles(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        cols=[0.6 0.4 0.6 0.9 0.6 0.6];
        if prot(6,1)<12000
            tt=12000;
            tim=[0; round(prot(2:6,1)); tt];
        else
            tt=14000;
            tim=[0; round(prot(2:6,1)); tt];
        end
        order=[1 4 2 3];
    else %quick
        first=goodfiles(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        timing=find(prot(:,4)==50);
        timing=round(prot(timing,1));
        timing=timing-2500;
        di=diff(timing);
        di=median(di);
        for t=2:length(timing)+1
            timing(t)=timing(t-1)+di;
        end
        cols=[0.6 0.9 0.6 0.4 0.6];
        tt=75000;
        tim=[0 500 2500 4500 9500 11500 13500];
        order=[2 3 1 4];
    end
    
    figure
    subplot(3,1,1)
    if ftype==1
        for t=1:length(tim)-1
            rectangle('position',[tim(t),0,tim(t+1)-tim(t),length(goodfiles)+1],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
            hold on
        end
    else
        for t=1:length(tim)-1
            rectangle('position',[tim(t),0,tim(t+1)-tim(t),length(goodfiles)*5+1],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
            hold on
        end
    end
    allrate=[]; counting=0; allrate2=[];
    for f=goodfiles
        spikes=unit{1,2}{f,2};
        rate=MEA_spikerates(spikes,40,tt);
        allrate=[allrate;rate];
        rate2=MEA_spikerates(spikes,80,tt);
        allrate2=[allrate2;rate2];
        if ftype==1
            counting=counting+1;
            for ss=1:length(spikes)
                subplot(3,1,1)
                plot([spikes(ss) spikes(ss)],[counting-0.4 counting+0.4],'-k')
                hold on
            end
        else
            for t=1:length(timing)-1
                counting=counting+1;
                tmp=find(spikes>=timing(t));
                tmp=tmp(1);
                if ~isempty(tmp)
                    tmp2=find(spikes<timing(t+1));
                    if ~isempty(tmp2)
                        tmp2=tmp2(end);
                        spi=spikes(tmp:tmp2);
                        for ss=1:length(spi)
                            subplot(3,1,1)
                            plot([spi(ss) spi(ss)],[counting-0.4 counting+0.4],'-k')
                            hold on
                        end
                        
                    end
                end
            end
        end
    end
    axis([0 tt 0 inf])
    meanrate=mean(allrate);
    meanrate2=mean(allrate2);
    
    if ftype==1
        mrate=meanrate;
        mrate2=meanrate2;
    else
        mrate=[];
        for t=1:length(timing)-1
            curr=meanrate(timing(t):timing(t+1));
            mrate=[mrate;curr];
        end
        mrate2=[];
        for t=1:length(timing)-1
            curr=meanrate2(timing(t):timing(t+1));
            mrate2=[mrate2;curr];
        end
    end
    
    subplot(3,1,2)
    for t=1:length(tim)-1
        rectangle('position',[tim(t),0,tim(t+1)-tim(t),max(mrate)+2],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
        hold on
    end
    plot(mrate,'-r','linewidth',2)
    axis([0 tt 0 max(mrate)+2])
    
    subplot(3,1,3)
    backgr=mean(mrate2(300:500));
    newrate=mrate2-backgr;
    testing=find(mrate2==0);
    di=diff(testing);
    tmp=find(di>1);
    if ~isempty(tmp)
        tmp=tmp(end)+1;
        gogo=testing(tmp);
    else
        gogo=testing(1);
    end
    newrate=newrate(300:gogo-400);
    if ~isempty(newrate)
        for t=1:length(tim)-1
            rectangle('position',[tim(t),min(newrate)-1,tim(t+1)-tim(t),max(newrate)-min(newrate)+3],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
            hold on
        end
        plot(300:length(newrate)+299,newrate,'-m','linewidth',2)
        title('baseline subtracted')
        axis([0 tt min(newrate)-1 max(newrate)+2])
    end
    
    if info{currID,8}==1
        po='ON';
    elseif info{currID,8}==-1
        po='OFF';
    elseif info{currID,8}==2
        po='ON-OFF';
    end
    if info{currID,9}==1
        po2='ON';
    elseif info{currID,9}==-1
        po2='OFF';
    elseif info{currID,9}==0
        po2='none';
    end
    subplot(3,1,1)
    title([currd,'_',curru,' // manual pol: ',po, ' // LF: ',po2],'interpreter','none')
    set(gcf,'position',[1601 1 1600 824])
    
    
    
    saveas(gcf,fullfile(superpath,'FFFlash',[Dtype,'_',int2str(i),'_',currd,'_',curru,'.bmp']))
    close
    
end
%% 12b) remove some cells form the flash group
% human
% badcells=[7 9 10 19 24 32 33 35 39 40 42 43 44 45 46 47 50 53 56 57 58 60 61 62 63 64 66 67 69 70 71 72 74 76 77 78 79 80 82 87 88 90 92 94 101 110 113 115 116 118 119 120 121 128 129 131 133 136 138 140 141 143 145 146 149 151 153 154 158 168 173 176 177 178 179 181 186 188 191 192 194 195 197 198 201 204 205 208 209 210 211 212 214 215 216 217 218];

%pig
% badcells=[18 19 21 26 48 51 52 54 57 65 69 75 77 79 81 82 83 85 90 91 173 183 185 190 194 202 211 215 218 219 222 223 224 227 229 230 231 232];

%mouse
badcells=[7 8 9 28 40 42 58 65];

load(fullfile(superpath,[Dtype,'_info']))
save(fullfile(superpath,[Dtype,'_info_before11']),'info','goodones')
togo=find(goodones>0);
tmp=cell2mat(info(togo,8));
tmp2=find(tmp~=0);
togo=togo(tmp2);

badIDs=togo(badcells);
for b=1:length(badIDs)
    info{badIDs(b),8}=0;
end
save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
%% read out some values: new
clear
close all

Dtype='w';
superpath='C:\Users\Katja\Desktop\all_species_temp';
load(fullfile(superpath,[Dtype,'_info']))

good=find(goodones==1);
fullfield=cell2mat(info(good,8));
fullfield=length(find(fullfield~=0));
LF=cell2mat(info(good,9));
LF=length(find(LF~=0));
DS=cell2mat(info(good,12));
DS=length(find(DS==1));

vel61=cell2mat(info(good,16));
vel61=find(vel61>0);
vel62=cell2mat(info(good,17));
vel62=find(vel62>0);
vel6=length(unique([vel61;vel62]));
DG=cell2mat(info(good,18));
DG=length(find(DG>0));
chirp=cell2mat(info(good,19));
chirp=length(find(chirp>0));
%% 12c) more on full-field: plot again after corrections
% col7) #quick files
% col8) pol from flash: 0=none, -1=off, 1=on, 2=onoff

load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
if ~exist(fullfile(superpath,'FFFlash_new'),'dir')
    mkdir(fullfile(superpath,'FFFlash_new'))
end

togo=find(goodones>0);
tmp=cell2mat(info(togo,8));
tmp2=find(tmp~=0);
togo=togo(tmp2);
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
    
    quicks=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'FFFlash'))
            ftype=1;
            quicks=[quicks;1];
        elseif ~isempty(regexp(hekalist{f},'quick'))
            ftype=2;
            %             currd
            quicks=[quicks;1];
        else
            quicks=[quicks;0];
        end
    end
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(quicks==1);
    if Dtype=='m'
        goodfiles=goodfiles(11:end);
    else
        goodfiles=goodfiles(5:end);
    end
    
    if ftype==1 %FFFlash
        first=goodfiles(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        cols=[0.6 0.4 0.6 0.9 0.6 0.6];
        if prot(6,1)<12000
            tt=12000;
            tim=[0; round(prot(2:6,1)); tt];
        else
            tt=14000;
            tim=[0; round(prot(2:6,1)); tt];
        end
        order=[1 4 2 3];
    else %quick
        first=goodfiles(1);
        prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
        timing=find(prot(:,4)==50);
        timing=round(prot(timing,1));
        timing=timing-2500;
        di=diff(timing);
        di=median(di);
        for t=2:length(timing)+1
            timing(t)=timing(t-1)+di;
        end
        cols=[0.6 0.9 0.6 0.4 0.6];
        tt=75000;
        tim=[0 500 2500 4500 9500 11500 13500];
        order=[2 3 1 4];
    end
    
    figure
    subplot(3,1,1)
    if ftype==1
        for t=1:length(tim)-1
            rectangle('position',[tim(t),0,tim(t+1)-tim(t),length(goodfiles)+1],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
            hold on
        end
    else
        for t=1:length(tim)-1
            rectangle('position',[tim(t),0,tim(t+1)-tim(t),length(goodfiles)*5+1],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
            hold on
        end
    end
    allrate=[]; counting=0; allrate2=[];
    for f=goodfiles
        spikes=unit{1,2}{f,2};
        
        rate=MEA_spikerates(spikes,40,tt);
        allrate=[allrate;rate];
        rate2=MEA_spikerates(spikes,80,tt);
        allrate2=[allrate2;rate2];
        if ftype==1
            counting=counting+1;
            all_spikes{i,counting}=spikes;
            for ss=1:length(spikes)
                subplot(3,1,1)
                plot([spikes(ss) spikes(ss)],[counting-0.4 counting+0.4],'-k')
                hold on
            end
        else
            for t=1:length(timing)-1
                counting=counting+1;
                tmp=find(spikes>=timing(t));
                tmp=tmp(1);
                if ~isempty(tmp)
                    tmp2=find(spikes<timing(t+1));
                    if ~isempty(tmp2)
                        tmp2=tmp2(end);
                        spi=spikes(tmp:tmp2);
                        all_spikes{i,counting}=spi;
                        for ss=1:length(spi)
                            subplot(3,1,1)
                            plot([spi(ss) spi(ss)],[counting-0.4 counting+0.4],'-k')
                            hold on
                        end
                        
                    end
                end
            end
        end
    end
    axis([0 tt 0 inf])
    meanrate=mean(allrate);
    flash_rates{i}=meanrate;
    meanrate2=mean(allrate2);
    
    if ftype==1
        mrate=meanrate;
        mrate2=meanrate2;
    else
        mrate=[];
        for t=1:length(timing)-1
            curr=meanrate(timing(t):timing(t+1));
            mrate=[mrate;curr];
        end
        mrate2=[];
        for t=1:length(timing)-1
            curr=meanrate2(timing(t):timing(t+1));
            mrate2=[mrate2;curr];
        end
    end
    
    subplot(3,1,2)
    for t=1:length(tim)-1
        rectangle('position',[tim(t),0,tim(t+1)-tim(t),max(mrate)+2],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
        hold on
    end
    plot(mrate,'-r','linewidth',2)
    axis([0 tt 0 max(mrate)+2])
    
    subplot(3,1,3)
    backgr=mean(mrate2(300:500));
    newrate=mrate2-backgr;
    testing=find(mrate2==0);
    di=diff(testing);
    tmp=find(di>1);
    if ~isempty(tmp)
        tmp=tmp(end)+1;
        gogo=testing(tmp);
    else
        gogo=testing(1);
    end
    newrate=newrate(300:gogo-400);
    if ~isempty(newrate)
        for t=1:length(tim)-1
            rectangle('position',[tim(t),min(newrate)-1,tim(t+1)-tim(t),max(newrate)-min(newrate)+3],'facecolor',[cols(t) cols(t) cols(t)],'edgecolor',[cols(t) cols(t) cols(t)])
            hold on
        end
        plot(300:length(newrate)+299,newrate,'-m','linewidth',2)
        title('baseline subtracted')
        axis([0 tt min(newrate)-1 max(newrate)+2])
    end
    
    if info{currID,8}==1
        po='ON';
    elseif info{currID,8}==-1
        po='OFF';
    elseif info{currID,8}==2
        po='ON-OFF';
    end
    if info{currID,9}==1
        po2='ON';
    elseif info{currID,9}==-1
        po2='OFF';
    elseif info{currID,9}==0
        po2='none';
    end
    subplot(3,1,1)
    title([currd,'_',curru,' // manual pol: ',po, ' // LF: ',po2],'interpreter','none')
    set(gcf,'position',[1601 1 1600 824])
    
    
    
    saveas(gcf,fullfile(superpath,'FFFlash_new',[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
    close
    
end
save(fullfile(superpath,[Dtype,'_flash_new']),'flash_rates','all_spikes','togo')
%% 12d) plot mean full-field responses
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
typ='wph';
for n=1:3
    Dtype=typ(n);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_flashrates']))
    if ~isempty(regexp(Dtype,'w'))
        col='b';
    elseif  ~isempty(regexp(Dtype,'p'))
        col='g';
    else
        col='r';
    end
    
    lf=9; flash=8;
    
    lftype=[1 -1 0];
    lftext={'LF ON','LF OFF','no LF'};
    flashtype=[1 -1 2 0]; %2=onoff
    flashtext={' / Flash ON',' / Flash OFF',' / Flash ONOFF',' / no Flash'};
    
    consider=find(goodones==1);
    flashdata=cell2mat(info(consider,flash));
    lfdata=cell2mat(info(consider,lf));
    
    figure
    for l=1:length(lftype)
        currtype=lftype(l);
        
        if currtype~=0
            data=consider(find(lfdata==currtype));
            subplot(3,2,l)
            allrate=[];
            for dd=1:length(data)
                if dd==1 || length(flash_rates{data(dd)})==size(allrate,2)
                    currrate=flash_rates{data(dd)};
                    backgr=mean(currrate(300:550));
                    currrate=currrate-backgr;
                    allrate=[allrate;currrate];
                end
            end
            if ~isempty(data)
                plot(median(allrate(:,300:10500),1),'-','color',col,'linewidth',3)
                hold on
                %             if size(allrate,1)>1
                %                 plot(max(allrate(:,300:10500)),'-','color',col,'linewidth',1)
                %                 plot(min(allrate(:,300:10500)),'-','color',col,'linewidth',1)
                %             end
            end
            title([lftext{l},' / all flash polarities (n=',int2str(dd),')'])
        else
            for ty=1:length(flashtype)
                data=find(lfdata==currtype);
                tmp=find(flashdata(data)==flashtype(ty));
                data=consider(data(tmp));
                subplot(3,2,ty+2)
                allrate=[];
                for dd=1:length(data)
                    if dd==1 || length(flash_rates{data(dd)})==size(allrate,2)
                        currrate=flash_rates{data(dd)};
                        backgr=mean(currrate(300:550));
                        currrate=currrate-backgr;
                        allrate=[allrate;currrate];
                    end
                end
                if ~isempty(data)
                    plot(median(allrate(:,300:10500),1),'-','color',col,'linewidth',3)
                    hold on
                    %                 if size(allrate,1)>1
                    %                     plot(max(allrate(:,300:10500)),'-','color',col,'linewidth',1)
                    %                     plot(min(allrate(:,300:10500)),'-','color',col,'linewidth',1)
                    %                 end
                end
                title([lftext{l},flashtext{ty},' (n=',int2str(dd),')'])
            end
        end
        
    end
end
%% 12e) extract and plot latencies (only shortest per flash stim) no version for quick stimuli yet!!
clear
close all
savepath='C:\Users\Katja\Desktop\all_species_temp\overview_pics';

superpath='C:\Users\Katja\Desktop\all_species_temp';
typ='wph';
for n=1:3
    Dtype=typ(n);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_flashrates']))
    if ~isempty(regexp(Dtype,'w'))
        col='b';
    elseif  ~isempty(regexp(Dtype,'p'))
        col='g';
    else
        col='r';
    end
    
    lf=9; flash=8;
    
    lftype=[1 -1 0];
    lftext={'LF ON','LF OFF','no LF'};
    flashtype=[1 -1 2 0]; %2=onoff
    flashtext={' / Flash ON',' / Flash OFF',' / Flash ONOFF',' / no Flash'};
    
    consider=find(goodones==1);
    flashdata=cell2mat(info(consider,flash));
    lfdata=cell2mat(info(consider,lf));
    
    
    for l=1:length(lftype)
        currtype=lftype(l);
        
        if currtype~=0
            data=consider(find(lfdata==currtype));
            datesnow=info(data,1);
            
            maxs=[]; lats=[];
            if ~isempty(data)
                for dd=1:length(data)
                    currrate=flash_rates{data(dd)};
                    backgr=mean(currrate(300:550));
                    currrate=currrate-backgr;
                    if length(currrate)<13000
                        tim=[1984 3968 5952 7937];
                    else
                        tim=[2501 5002 7503 10004];
                    end
                    allmax=[]; alllat=[];
                    for ttt=1:length(tim)
                        [ma lat]=max(currrate(tim(ttt)+50:tim(ttt)+600));
                        if lat<70
                            if mean(currrate(tim(ttt)+50+lat-100:tim(ttt)+50+lat-10))+2*std(currrate(tim(ttt)+50+lat-100:tim(ttt)+50+lat-10))>=ma
                                lat=1000; ma=0;
                            end
                        end
                        allmax=[allmax ma];
                        alllat=[alllat lat+50];
                    end
                    [allmax id]=sort(allmax,'descend');
                    alllat=alllat(id);
                    [lat id]=min(alllat(1:2));
                    ma=allmax(id);
                    maxs=[maxs;ma];
                    lats=[lats;lat];
                end
                %             end
                
                numPerDate=[];
                prev=datesnow{1};
                for da=2:length(datesnow)
                    if isempty(regexp(datesnow{da},prev))
                        numPerDate=[numPerDate;da];
                        prev=datesnow{da};
                    end
                end
                numPerDate=[numPerDate;da];
                
                figure(1)
                subplot(3,2,l)
                for ll=1:length(lats)
                    plot(ll,lats(ll),'o','color',col)
                    hold on
                end
                old=1;
                for nn=1:length(numPerDate)
                    plot([numPerDate(nn)-0.5 numPerDate(nn)-0.5],[min(lats)-1 max(lats)+1],'-k')
                    currlats=median(lats(old:numPerDate(nn)-1));
                    plot([old numPerDate(nn)-1],[currlats currlats],'-','color',col)
                    old=numPerDate(nn);
                end
                [lats id]=sort(lats);
                maxs=maxs(id);
                figure(2)
                subplot(3,2,l)
                for ll=1:length(lats)
                    plot(ll,lats(ll),'o','color',col)
                    hold on
                end
                %                 old=0.5;
                %                 for nn=1:length(numPerDate)
                %                     plot([numPerDate(nn)+old numPerDate(nn)+old],[0 max(lats)+1],'-k')
                %                     old=numPerDate(nn)+old;
                %                 end
            end
            figure(1)
            title([lftext{l},' / all flash polarities (n=',int2str(length(data)),')'])
            figure(2)
            title([lftext{l},' / all flash polarities (n=',int2str(length(data)),')'])
            
            
            
        else
            for ty=1:length(flashtype)
                data=find(lfdata==currtype);
                tmp=find(flashdata(data)==flashtype(ty));
                data=consider(data(tmp));
                datesnow=info(data,1);
                
                maxs=[]; lats=[];
                if ~isempty(data)
                    for dd=1:length(data)
                        currrate=flash_rates{data(dd)};
                        backgr=mean(currrate(300:550));
                        currrate=currrate-backgr;
                        if length(currrate)<13000
                            tim=[1984 3968 5952 7937];
                        else
                            tim=[2501 5002 7503 10004];
                        end
                        allmax=[]; alllat=[];
                        for ttt=1:length(tim)
                            [ma lat]=max(currrate(tim(ttt)+50:tim(ttt)+600));
                            if lat<70
                                if mean(currrate(tim(ttt)+50+lat-100:tim(ttt)+50+lat-10))+2*std(currrate(tim(ttt)+50+lat-100:tim(ttt)+50+lat-10))>=ma
                                    lat=1000; ma=0;
                                end
                            end
                            allmax=[allmax ma];
                            alllat=[alllat lat+50];
                        end
                        [allmax id]=sort(allmax,'descend');
                        alllat=alllat(id);
                        [lat id]=min(alllat(1:2));
                        ma=allmax(id);
                        maxs=[maxs;ma];
                        lats=[lats;lat];
                    end
                    %             end
                    
                    numPerDate=[];
                    prev=datesnow{1};
                    for da=2:length(datesnow)
                        if isempty(regexp(datesnow{da},prev))
                            numPerDate=[numPerDate;da];
                            prev=datesnow{da};
                        end
                    end
                    numPerDate=[numPerDate;da];
                    
                    figure(1)
                    subplot(3,2,ty+2)
                    for ll=1:length(lats)
                        plot(ll,lats(ll),'o','color',col)
                        hold on
                    end
                    old=1;
                    for nn=1:length(numPerDate)
                        plot([numPerDate(nn)-0.5 numPerDate(nn)-0.5],[min(lats)-1 max(lats)+1],'-k')
                        currlats=median(lats(old:numPerDate(nn)-1));
                        plot([old numPerDate(nn)-1],[currlats currlats],'-','color',col)
                        old=numPerDate(nn);
                    end
                    [lats id]=sort(lats);
                    maxs=maxs(id);
                    figure(2)
                    subplot(3,2,ty+2)
                    for ll=1:length(lats)
                        plot(ll,lats(ll),'o','color',col)
                        hold on
                    end
                    %                       old=0.5;
                    %                 for nn=1:length(numPerDate)
                    %                     plot([numPerDate(nn)+old numPerDate(nn)+old],[0 max(lats)+1],'-k')
                    %                     old=numPerDate(nn)+old;
                    %                 end
                end
                figure(1)
                title([lftext{l},flashtext{ty},' (n=',int2str(dd),')'])
                figure(2)
                title([lftext{l},flashtext{ty},' (n=',int2str(dd),')'])
            end
            
        end
    end
    figure(1)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,[Dtype,'_latencyPerExp.bmp']))
    figure(2)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath,[Dtype,'_latencySorted.bmp']))
    close all
end
%% 13a) re-analysis of 6vel: without slowest speed

load(fullfile(superpath,[Dtype,'_info']))
currtype='Black';
idnow=16; %16 for black, 17 for white
load(fullfile(superpath,[Dtype,'_vel6_new_',currtype]))
cd(superpath)

names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[1 2 4 8 16];

peak6vel=zeros(size(info,1),6);
togo=[];
for i=1:length(info)
    if ~isempty(info{i,idnow})
        if info{i,idnow}>0
            togo=[togo;i];
        end
    end
end
for i=1:length(togo)
    i
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    
    
    if idnow==16
        avRate=vel6_newB{currID};
    else
        avRate=vel6_newW{currID};
    end
    if isempty(avRate)
        info{currID,idnow}=[];
    else
        
        procpath=fullfile(mainpath,currd,'processing');
        hpath=fullfile(mainpath,currd,'HEKA');
        load(fullfile(procpath,'stimuliList'))
        
        
        m=0; found=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(hekaNames{m},'bar_6veloc_downup_'))
                found=1;
            end
        end
        prot=read_header_field_heka(hpath,hekaNames{m},'Stimulus Protocol');
        cons=prot(1:end,2);
        changes=find(cons==-1);
        startpoints=round(prot(changes(1:7),1));
        
        
        avbackground=mean(avRate(200:startpoints(1)));
        clear data3
        for j=1:6
            [data3(j) data4(j)]=max(avRate(startpoints(j):startpoints(j+1)));
        end
        [maxval maxid]=max(data3);
        data3=data3-avbackground;
        peak6vel(currID,:)=data3;
        startpoints=startpoints(2:end);
        data3=data3(2:end);
        data3=(data3-min(data3))/(max(data3)-min(data3));
        
        
        
        total=sum(data3); sofar=0; nn=0;
        if total>0
            while sofar<total/2
                nn=nn+1;
                sofar=sofar+data3(nn);
            end
            if data3(nn)==sofar && sum(data3(nn+1:end))==0
                position=speeds(nn);
                value=total;
                cutoff=total;
            elseif data3(nn)==sofar && nn==1
                position=speeds(nn);
                value=total/2;
                cutoff=total/2;
            else
                cutoff=total/2;
                a=nn-1; b=nn;
                before=sum(data3(1:a));
                after=sum(data3(1:b));
                slope=(after-before)/(speeds(b)-speeds(a));
                position=(slope*speeds(b)-after+total/2)/slope;
                c=find(speeds<position);
                c=c(end);
                cc=c+1;
                slope2=(data3(cc)-data3(c))/(speeds(cc)-speeds(c));
                value=data3(cc)-slope2*speeds(cc)+position*slope2;
            end
            summed_coll(1)=data3(1);
            for rr=2:length(data3)
                summed_coll(rr)=summed_coll(rr-1)+data3(rr);
            end
        else
            position=0;
            value=0;
            summed_coll=repmat(0,length(data3),1);
            cutoff=total/2;
        end
        if idnow==16
            info{currID,24}=position;
        else
            info{currID,25}=position;
        end
    end
end

save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
if idnow==16
    save(fullfile(superpath,[Dtype,'_peaks6vel_black']),'peak6vel')
else
    save(fullfile(superpath,[Dtype,'_peaks6vel_white']),'peak6vel')
end
%% 13b) PLOT histograms for speed tuning
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');


names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[0.2; 1; 2; 4; 8; 16];
types='wph';
cols='bgr';
titles={'mouse','pig','human'};
figure(1)
figure(2)
for t=1:length(types)
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    
    goodB=[]; goodW=[];
    for i=1:length(info)
        if ~isempty(info{i,24})
            if info{i,24}>0
                goodB=[goodB;i];
            end
        end
        if ~isempty(info{i,25})
            if info{i,25}>0
                goodW=[goodW;i];
            end
        end
    end
    
    black=cell2mat(info(goodB,24));
    if t==1
        mblack=black;
    elseif t==2
        pblack=black;
    else
        hblack=black;
    end
    white=cell2mat(info(goodW,25));
    if t==1
        mwhite=white;
    elseif t==2
        pwhite=white;
    else
        hwhite=white;
    end
    black=[black;-1;17];
    white=[white;-1;17];
    
    figure(1)
    subplot(3,2,t*2-1)
    hist(black,38)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
    med=median(black(1:end-2));
    st=std(black(1:end-2));
    conf = 1.96*st/sqrt(length(black)-2);
    title([titles{t}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(black)-2),')'])
    hold on
    plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
    plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
    yy=get(gca,'ylim');
    plot([med-conf med-conf],[0 yy(2)+1],'--k')
    plot([med+conf med+conf],[0 yy(2)+1],'--k')
    %     set(gca,'xtick',speeds)
    % set(gca,'xticklabel',speeds)
    xlabel('speed (mm/s)')
    axis([0 16 0 yy(2)+1])
    
    subplot(3,2,[2 4])
    c=cdfplot(black(1:end-2))
    set(c,'color',cols(t),'linestyle','-','linewidth',2)
    hold on
    xlabel('speed (mm/s)')
    ylabel('probability')
    title('black bar')
    
    figure(2)
    subplot(3,2,t*2-1)
    hist(white,38)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
    med=median(white(1:end-2));
    st=std(white(1:end-2));
    conf = 1.96*st/sqrt(length(black)-2);
    title([titles{t}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(white)-2),')'])
    hold on
    plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
    plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
    xlabel('speed (mm/s)')
    yy=get(gca,'ylim');
    plot([med-conf med-conf],[0 yy(2)+1],'--k')
    plot([med+conf med+conf],[0 yy(2)+1],'--k')
    axis([0 16 0 yy(2)+1])
    %      set(gca,'xtick',speeds)
    %     set(gca,'xticklabel',speeds)
    
    subplot(3,2,[2 4])
    c=cdfplot(white(1:end-2))
    set(c,'color',cols(t),'linestyle','-','linewidth',2)
    hold on
    xlabel('speed (mm/s)')
    ylabel('probability')
    title('white bar')
    %      set(gca,'xtick',speeds)
    % set(gca,'xticklabel',speeds)
end
figure(1)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'speed_black_noSlow.bmp'))
saveas(gcf,fullfile(savepath,'speed_black_noSlow.fig'))
close
figure(2)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'speed_white_noSlow.bmp'))
saveas(gcf,fullfile(savepath,'speed_white_noSlow.fig'))
close
%% 13c) PLOT 6vel separated by ON/OFF

clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

types='wph';
cols='bgr';
titles={'mouse','pig','human'};

lfs={'ON LF','OFF LF','no LF'};
flashes={'ON manual','OFF manual','none','ON-OFF manual'};

flash=8; LF=9;
black=24; %noSlow: 24, all: 16
white=25; %noSlow: 25, all: 17
currtype='noSlow';
% currtype='all';


for t=1:length(types)
    Dtype=types(t);
    load(fullfile(superpath,[Dtype,'_info']))
    
    
    
    ids=find(goodones==1);
    lfdata=cell2mat(info(ids,LF));
    flashdata=cell2mat(info(ids,flash));
    blackdata=[]; whitedata=[];
    for i=1:length(ids)
        if ~isempty(info{ids(i),black})
            blackdata=[blackdata;info{ids(i),black}];
        else
            blackdata=[blackdata;0];
        end
        if ~isempty(info{ids(i),white})
            whitedata=[whitedata;info{ids(i),white}];
        else
            whitedata=[whitedata;0];
        end
    end
    
    counting=0;
    for currlf=[1 -1 0]
        figure(t)
        counting=counting+1;
        lfnow=find(lfdata==currlf);
        bl=blackdata(lfnow);
        wh=whitedata(lfnow);
        
        tmp=find(bl>0); bl=bl(tmp);
        tmp=find(wh>0); wh=wh(tmp);
        
        bl=[bl;-1;17];
        wh=[wh;-1;17];
        
        
        subplot(4,4,counting)
        hist(bl,38)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
        med=median(bl(1:end-2));
        st=std(bl(1:end-2));
        conf = 1.96*st/sqrt(length(bl)-2);
        title(['BLACK: ',lfs{counting}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(bl)-2),')'])
        hold on
        plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
        plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
        yy=get(gca,'ylim');
        plot([med-conf med-conf],[0 yy(2)+1],'--k')
        plot([med+conf med+conf],[0 yy(2)+1],'--k')
        xlabel('speed (mm/s)')
        axis([0 16 0 yy(2)+1])
        
        
        subplot(4,4,counting+8)
        hist(wh,38)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
        med=median(wh(1:end-2));
        st=std(wh(1:end-2));
        conf = 1.96*st/sqrt(length(wh)-2);
        title(['WHITE: ',lfs{counting}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(wh)-2),')'])
        hold on
        plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
        plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
        yy=get(gca,'ylim');
        plot([med-conf med-conf],[0 yy(2)+1],'--k')
        plot([med+conf med+conf],[0 yy(2)+1],'--k')
        xlabel('speed (mm/s)')
        axis([0 16 0 yy(2)+1])
        
        figure(101)
        subplot(4,4,counting)
        [c stats]=cdfplot(bl(1:end-2))
        set(c,'color',cols(t),'linestyle','-','linewidth',2)
        hold on
        title(['BLACK: ',lfs{counting}])
        
        
        figure(101)
        if length(wh)>2
            subplot(4,4,counting+8)
            [c stats]=cdfplot(wh(1:end-2))
            set(c,'color',cols(t),'linestyle','-','linewidth',2)
            hold on
            title(['WHITE: ',lfs{counting}])
        end
        
    end
    
    
    counting=0;
    for currflash=[1 -1 0 2]
        figure(t)
        counting=counting+1;
        flashnow=find(flashdata==currflash);
        bl=blackdata(flashnow);
        wh=whitedata(flashnow);
        
        tmp=find(bl>0); bl=bl(tmp);
        tmp=find(wh>0); wh=wh(tmp);
        
        bl=[bl;-1;17];
        wh=[wh;-1;17];
        
        
        subplot(4,4,counting+4)
        hist(bl,38)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
        med=median(bl(1:end-2));
        st=std(bl(1:end-2));
        conf = 1.96*st/sqrt(length(bl)-2);
        title(['BLACK: ',flashes{counting}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(bl)-2),')'])
        hold on
        plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
        plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
        yy=get(gca,'ylim');
        plot([med-conf med-conf],[0 yy(2)+1],'--k')
        plot([med+conf med+conf],[0 yy(2)+1],'--k')
        xlabel('speed (mm/s)')
        axis([0 16 0 yy(2)+1])
        
        
        subplot(4,4,counting+12)
        hist(wh,38)
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
        med=median(wh(1:end-2));
        st=std(wh(1:end-2));
        conf = 1.96*st/sqrt(length(wh)-2);
        title(['WHITE: ',flashes{counting}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(wh)-2),')'])
        hold on
        plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
        plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
        yy=get(gca,'ylim');
        plot([med-conf med-conf],[0 yy(2)+1],'--k')
        plot([med+conf med+conf],[0 yy(2)+1],'--k')
        xlabel('speed (mm/s)')
        axis([0 16 0 yy(2)+1])
        
        figure(101)
        subplot(4,4,counting+4)
        [c stats]=cdfplot(bl(1:end-2))
        set(c,'color',cols(t),'linestyle','-','linewidth',2)
        hold on
        title(['BLACK: ',flashes{counting}])
        
        
        figure(101)
        
        subplot(4,4,counting+12)
        if length(wh)>2
            [c stats]=cdfplot(wh(1:end-2))
            set(c,'color',cols(t),'linestyle','-','linewidth',2)
            hold on
            title(['WHITE: ',flashes{counting}])
        end
        
    end
    %     figure(1)
    %     set(gcf,'position',[1 31 1600 794])
    %     saveas(gcf,fullfile(savepath,['speed_polaritySep_',currtype,'_',Dtype,'.bmp']))
    %     saveas(gcf,fullfile(savepath,['speed_polaritySep_',currtype,'_',Dtype,'.fig']))
    %     close
end
%% 14a) chirp FFT


load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_chirp']))
cd(superpath)
savepath=fullfile(superpath,'chirp_FFT');

togo=[];
for i=1:length(info)
    if ~isempty(info{i,19})
        if info{i,19}>0
            togo=[togo;i];
        end
    end
end

chirpFFTabs=cell(length(togo),5);
chirpFFTsmooth=cell(length(togo),5);
chirpRates=cell(length(togo),2);
sicnt=0;

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
    
    prot=read_header_field_heka(hpath,hekalist{goodfiles(1)},'Stimulus Protocol');
    if prot(end,1)>24000
        goodfiles=[];
    end
    
    
    if ~isempty(goodfiles)
        sicnt=0;
        for sigma=[5 10 20 40 60]
            sicnt=sicnt+1;
            allrate=[];
            for f=goodfiles
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,sigma,25000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            chirpRates{i,1}=meanrate;
            
            
            int=prot(3:end-1,4);
            
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
            %         int2=[int2;repmat(128,missing,1)];
            
            %         %         smooint=[]; last=1;
            %         %         for iii=1:length(idds)-1
            %         %             curr1=int2(last);
            %         %             curr2=int2(last+idds(iii));
            %         %             newval2=mean([curr1 curr2]);
            %         %             newval1=mean([curr1 newval2]);
            %         %             newval3=mean([curr2 newval2]);
            %         %             each=round(idds(iii)/4);
            %         %             final=idds(iii)-3*each;
            %         %             smooint=[smooint;int2(last:last+each);repmat(newval1,each,1);repmat(newval2,each,1);repmat(newval3,final,1)];
            %         %
            %         %             %                newval2=mean([curr1 curr2]);
            %         %             %             smooint=[smooint;int2(last:last+round(idds(iii)/2));repmat(newval2,idds(iii)-round(idds(iii)/2),1)];
            %         %             last=last+idds(iii);
            %         %         end
            %         %         smooint=[smooint;int2(last:last+idds(iii)-1)];
            NFFT=9000;
            %         NFFT=2^nextpow2(numel(int2));
            % int2=int2(1:7900);
            Y=fft(int2,NFFT)/numel(int2);
            fS=1000/2*linspace(0,1,NFFT/2+1);
            cS=2*abs(Y(1:NFFT/2+1));
            
            %         %         NFFT=2^nextpow2(numel(smooint)); %stimulus smoothing
            %         %         Y=fft(smooint,NFFT)/numel(smooint);
            %         %         fS=1000/2*linspace(0,1,NFFT/2+1);
            %         %         cS=2*abs(Y(1:NFFT/2+1));
            
            
            rate=meanrate(round(time(1)):round(time(1))+7999);
            %         rate=meanrate(round(time(1)):round(time(firstpart)));
            NFFT=9000;
            %         NFFT=2^nextpow2(numel(rate));
            Y=fft(rate,NFFT)/numel(rate);
            f=1000/2*linspace(0,1,NFFT/2+1);
            c=2*abs(Y(1:NFFT/2+1));
            c=c';
            
            %         figure
            %         subplot(3,2,1)
            %         plot(rate)
            %         title('firing rate')
            %         subplot(3,2,3)
            %         plot(f,c)
            tmp1=find(f>0.8);
            tmp1=tmp1(1);
            tmp2=find(f<32);
            tmp2=tmp2(end);
            tmp3=find(f<=7.5);
            tmp3=tmp3(end);
            %         axis([0 32 0 max(c(tmp1:tmp2))])
            %         title('FFT response')
            %         subplot(3,2,5)
            %         plot(fS,cS,'-k')
            %         axis([0 32 0 max(cS(tmp1:tmp2))])
            %         title('FFT stimulus')
            %         %         subplot(3,2,2)
            %         %         plot(f,c)
            %         %         hold on
            %         %         plot(f,cS,'-k')
            %         %         axis([0 32 0 max(max(cS(tmp1:tmp2)),max(c(tmp1:tmp2)))])
            
            %         subplot(3,2,2)
            %         plot(f,c./cS,'-r')
            %         axis([0 32 0 max(c(tmp1:tmp2)./cS(tmp1:tmp2))])
            %         title('absolute')
            chirpFFTabs{i,1}=f;
            chirpFFTabs{i,sicnt+1}=c./cS;
            
            %         subplot(3,2,4)
            c=c(tmp1:tmp3); cS=cS(tmp1:tmp3);
            cc=c/max(c); ccS=cS/max(cS); ratio=cc./ccS;
            %         plot(f(tmp1:tmp3),rat,'-r')
            %         axis([0 8 0 max(rat)])
            %         title('normalized')
            
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
                fr=mean(f(to+tmp1-1:endval+tmp1-1));
                smooF=[smooF;fr];
            end
            
            %         subplot(3,2,6)
            %         plot(smooF,smoo,'-g')
            %         hold on
            %         title('smooth')
            
            %         axis tight
            
            chirpFFTsmooth{i,1}=smooF;
            chirpFFTsmooth{i,sicnt+1}=smoo;
            
            %         set(gcf,'position',[1 31 1600 794])
            %         saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
            %         close
            
            
        end
    end
end
save(fullfile(superpath,[Dtype,'_chirpFFT']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates')
%% 14b) PLOT

clear

siglist=[5 10 20 40 60];
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

types='wph';
cols='bgr';
titles={'mouse','pig','human'};

for sig=1:5;
    figure
    for tt=1:length(types)
        Dtype=types(tt);
        load(fullfile(superpath,[Dtype,'_info']))
        load(fullfile(superpath,[Dtype,'_chirpFFT']))
        
        
        %     figure
        coll=[]; collsmooth=[];
        currFr=chirpFFTabs{1,1};
        tmp1=find(currFr>0.8);
        tmp2=find(currFr<=7.5);
        currFr=currFr(tmp1(1):tmp2(end));
        for i=1:length(chirpFFTsmooth)
            currdata=chirpFFTabs{i,sig+1};
            currdata=currdata(tmp1(1):tmp2(end));
            coll=[coll currdata];
            subplot(4,3,tt)
            plot(currFr,currdata,'-','color',[0.6 0.6 0.6])
            hold on
            
            currdata=chirpFFTsmooth{i,sig+1};
            collsmooth=[collsmooth currdata];
            subplot(4,3,tt+6)
            plot(chirpFFTsmooth{i,1},currdata,'-','color',[0.6 0.6 0.6])
            hold on
        end
        subplot(4,3,tt+3)
        plot(currFr,median(coll,2),'-','color',cols(tt),'linewidth',3)
        axis tight
        
        subplot(4,3,tt+9)
        plot(chirpFFTsmooth{i,1},median(collsmooth,2),'-','color',cols(tt),'linewidth',3)
        axis tight
        
        xlabel('Hz')
        title('freq. modulation')
    end
    subplot(4,3,1)
    title(['sigma=',int2str(siglist(sig))])
end
% set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(savepath,'chirp_modulation.bmp'))
% saveas(gcf,fullfile(savepath,'chirp_modulation.fig'))
% close
%% 15a) DG values: peaks of first 4 harmonics
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)


DGallpeaks=cell(length(info),1);
DGallfreqs=cell(length(info),1);
togo=find(goodDG==1);
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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    thiscell=[];
    fouriers=[];
    rates=[];
    clear f2take
    for g=1:length(DGlist)
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            
            stimFreq=freqs(g);
            
            if max(meanrate)>0
                meanF=mean(meanrate);
                meanFR=round(meanF*10)/10;
                % fourier
                NFFT=2^nextpow2(numel(meanrate));
                Y=fft(meanrate,NFFT)/numel(meanrate);
                f=1000/2*linspace(0,1,NFFT/2+1);
                if ~exist('f2take','var') && ~isempty(find(f>0))
                    f2take=f;
                end
                c=2*abs(Y(1:NFFT/2+1));
                
                
                ampPeak=max(c(4:end));
                ampPeakR=round(ampPeak);
                posF=find(c==ampPeak);
                posF=posF(1);
                freq=f(posF);
                
                allpeaks=[]; allfreqs=[];
                testfreq=stimFreq:stimFreq:40;
                testfreq=testfreq(1:4);
                for te=1:length(testfreq)
                    stimFreqN=testfreq(te);
                    findfreq1=find(f<stimFreqN);
                    findfreq1=findfreq1(end)-1;
                    if findfreq1<=0
                        findfreq1=4;
                    end
                    
                    findfreq2=max(c(findfreq1:findfreq1+5)); %recorded amp at stimfreq
                    if te==1
                        pos1=find(c==findfreq2);
                        pos1=f(pos1);
                    end
                    allpeaks=[allpeaks;findfreq2];
                    allfreqs=[allfreqs;pos1];
                end
            else
                allpeaks=[]; allfreqs=[];
                
            end
            DGallpeaks{currID,g}=allpeaks;
            DGallfreqs{currID,g}=allfreqs;
        end
    end
    
end
save(fullfile(superpath,[Dtype,'_DGallpeaks']),'DGallpeaks','DGallfreqs')
%% 15b) rearrange data
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'DG_harmonics');
savepath2=fullfile(superpath,'overview_pics');


types='wph';
cols='kbgr';
leg={'F1','F2','F3','F4'};
mas=[5 5 3 3]; liw=[2 2 1 1];
titles={'mouse','pig','human'};

for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_DGallpeaks']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(DGallpeaks{m,1})
            found=1;
        end
    end
    DGstim=[stimInfo{m,1} stimInfo{m,2}];
    spats=unique(DGstim(:,2));
    temps=unique(DGstim(:,1));
    spes=[];
    for ss=1:length(DGstim)
        noww=DGstim(ss,1)*DGstim(ss,2);
        spes=[spes;noww];
    end
    spes=spes/1000;
    spe=unique(spes);
    
    
    togo=find(goodDG==1);
    clear spatmax spatsum tempmax tempsum
    for i=1:length(togo)
        figure
        currID=togo(i);
        currs=DGallpeaks(currID,:);
        
        for harm=1:4
            curr=[];
            for cu=1:24
                if ~isempty(currs{cu})
                    curr=[curr;currs{cu}(harm)];
                else
                    curr=[curr;0];
                end
            end
            
            
            allsum=[]; allmax=[];
            for s=1:length(spats)
                tmp=find(DGstim(:,2)==spats(s));
                su=sum((curr(tmp)));
                ma=max(curr(tmp));
                allsum=[allsum su];
                allmax=[allmax ma];
            end
            spatsum{i,harm}=allsum; %cell,harmonic
            spatmax{i,harm}=allmax;
            
            allsum=[];allmax=[];
            for t=1:length(temps)
                tmp=find(DGstim(:,1)==temps(t));
                allDataTemp{t,1}=curr(tmp);
                su=sum(curr(tmp));
                ma=max(curr(tmp));
                allsum=[allsum su];
                allmax=[allmax ma];
            end
            tempsum{i,harm}=allsum;
            tempmax{i,harm}=allmax;
            
            
            subplot(2,4,1)
            plot(spats,spatsum{i,harm},'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
            hold on
            title('spatial (sum)')
            subplot(2,4,2)
            plot(spats,spatmax{i,harm},'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
            hold on
            title('spatial (max)')
            subplot(2,4,5)
            plot(temps,tempsum{i,harm},'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
            hold on
            title('temporal (sum)')
            subplot(2,4,6)
            plot(temps,tempmax{i,harm},'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
            hold on
            title('temporal (max)')
            
            curr=curr/max(curr);
            
            allsum=[]; allmax=[];
            for s=1:length(spats)
                tmp=find(DGstim(:,2)==spats(s));
                su=sum((curr(tmp)));
                ma=max(curr(tmp));
                allsum=[allsum su];
                allmax=[allmax ma];
            end
            
            subplot(2,4,3)
            plot(spats,allsum,'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
            hold on
            title('norm spatial (sum)')
            subplot(2,4,4)
            plot(spats,allmax,'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
            hold on
            title('norm spatial (max)')
            
            allsum=[];allmax=[];
            for t=1:length(temps)
                tmp=find(DGstim(:,1)==temps(t));
                allDataTemp{t,1}=curr(tmp);
                su=sum(curr(tmp));
                ma=max(curr(tmp));
                allsum=[allsum su];
                allmax=[allmax ma];
            end
            
            
            
            
            subplot(2,4,7)
            plot(temps,allsum,'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
            hold on
            title('norm temporal (sum)')
            subplot(2,4,8)
            plot(temps,allmax,'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
            hold on
            title('norm temporal (max)')
        end
        subplot(2,4,4)
        legend(leg,'location','southeast')
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'.bmp']))
        close
    end
    
    
    
    save(fullfile(superpath,[Dtype,'_DG_FFTdata']),'tempmax','tempsum','spatmax','spatsum')
end
%% 15c) plotting of F1/F2 comparisons
clear
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
load(fullfile(superpath,'DGinfo'))
savepath2=fullfile(superpath,'overview_pics');


types='wph';
cols='kbgr';
leg={'F1','F2','F3','F4'};
mas=[5 5 3 3]; liw=[2 2 1 1];
titles={'mouse','pig','human'};

for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_DG_FFTdata']))
    
    
    figure(1)
    subplot(2,2,1)
    for harm=1:4
        spatm=cell2mat(spatmax(:,harm));
        plot(spats,median(spatm),'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
        hold on
    end
    legend(leg,'location','southeast')
    title('population: spatial max')
    subplot(2,2,2)
    for harm=1:4
        spatm=cell2mat(spatsum(:,harm));
        plot(spats,median(spatm),'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
        hold on
    end
    title('population: spatial sum')
    subplot(2,2,3)
    for harm=1:4
        tempm=cell2mat(tempmax(:,harm));
        plot(temps,median(tempm),'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
        hold on
    end
    title('population: temporal max')
    subplot(2,2,4)
    for harm=1:4
        tempm=cell2mat(tempsum(:,harm));
        plot(temps,median(tempm),'o-','color',cols(harm),'markerfacecolor',cols(harm),'markersize',mas(harm),'linewidth',liw(harm))
        hold on
    end
    title('population: temporal sum')
    
    
    for ss=1:length(spats)
        figure(2)
        subplot(4,6,ss)
        f1=cell2mat(spatmax(:,1));
        f1=f1(:,ss);
        f2=cell2mat(spatmax(:,2));
        f2=f2(:,ss);
        scatter(f1,f2,'o')
        xlabel('F1')
        ylabel('F2')
        hold on
        plot([0 max(max(f1),max(f2))+1],[0 max(max(f1),max(f2))+1],'-k')
        axis([0 max(f1)+max(f1)*0.1 0 max(f2)+max(f2)*0.1])
        title(['spatial (max) ',int2str(spats(ss)),'um'])
        figure(3)
        subplot(4,6,ss)
        rat=f1./f2;
        clear hh
        tmp=find(rat>=1);
        h=hist(rat(tmp),40)
        hist(rat(tmp),40)
        hold on
        tmp=find(rat<1);
        h=hist(rat(tmp),40)
        hist(rat(tmp),40)
        hh = findobj(gca,'Type','patch');
        set(hh(1),'FaceColor','r','EdgeColor','r')
        perc=100/length(rat)*length(tmp);
        title(['max ',int2str(spats(ss)),'um (',int2str(perc),'% F2>F1)'])
        tmp=find(~isnan(rat));
        me=median(rat(tmp));
        legend(num2str(me))
        axis tight
        
        figure(2)
        subplot(4,6,ss+6)
        f1=cell2mat(spatsum(:,1));
        f1=f1(:,ss);
        f2=cell2mat(spatsum(:,2));
        f2=f2(:,ss);
        scatter(f1,f2,'o')
        xlabel('F1')
        ylabel('F2')
        hold on
        plot([0 max(max(f1),max(f2))+1],[0 max(max(f1),max(f2))+1],'-k')
        axis([0 max(f1)+max(f1)*0.1 0 max(f2)+max(f2)*0.1])
        title(['spatial (sum) ',int2str(spats(ss)),'um'])
        figure(3)
        subplot(4,6,ss+6)
        ratio=f1./f2;
        clear hh
        tmp=find(rat>=1);
        h=hist(rat(tmp),40)
        hist(rat(tmp),40)
        hold on
        tmp=find(rat<1);
        h=hist(rat(tmp),40)
        hist(rat(tmp),40)
        hh = findobj(gca,'Type','patch');
        set(hh(1),'FaceColor','r','EdgeColor','r')
        perc=100/length(rat)*length(tmp);
        title(['sum ',int2str(spats(ss)),'um (',int2str(perc),'% F2>F1)'])
        tmp=find(~isnan(rat));
        me=median(rat(tmp));
        legend(num2str(me))
        axis tight
        
    end
    
    for ss=1:length(temps)
        figure(2)
        subplot(4,6,ss+12)
        f1=cell2mat(tempmax(:,1));
        f1=f1(:,ss);
        f2=cell2mat(tempmax(:,2));
        f2=f2(:,ss);
        scatter(f1,f2,'o')
        xlabel('F1')
        ylabel('F2')
        hold on
        plot([0 max(max(f1),max(f2))+1],[0 max(max(f1),max(f2))+1],'-k')
        axis([0 max(f1)+max(f1)*0.1 0 max(f2)+max(f2)*0.1])
        title(['temporal (max) ',int2str(temps(ss)),'Hz'])
        figure(3)
        subplot(4,6,ss+12)
        ratio=f1./f2;
        clear hh
        tmp=find(rat>=1);
        h=hist(rat(tmp),40)
        hist(rat(tmp),40)
        hold on
        tmp=find(rat<1);
        h=hist(rat(tmp),40)
        hist(rat(tmp),40)
        hh = findobj(gca,'Type','patch');
        set(hh(1),'FaceColor','r','EdgeColor','r')
        perc=100/length(rat)*length(tmp);
        title(['max ',int2str(temps(ss)),'Hz (',int2str(perc),'% F2>F1)'])
        tmp=find(~isnan(rat));
        me=median(rat(tmp));
        legend(num2str(me))
        axis tight
        
        figure(2)
        subplot(4,6,ss+18)
        f1=cell2mat(tempsum(:,1));
        f1=f1(:,ss);
        f2=cell2mat(tempsum(:,2));
        f2=f2(:,ss);
        scatter(f1,f2,'o')
        xlabel('F1')
        ylabel('F2')
        hold on
        plot([0 max(max(f1),max(f2))+1],[0 max(max(f1),max(f2))+1],'-k')
        axis([0 max(f1)+max(f1)*0.1 0 max(f2)+max(f2)*0.1])
        title(['temporal (sum) ',int2str(temps(ss)),'Hz'])
        figure(3)
        subplot(4,6,ss+18)
        ratio=f1./f2;
        clear hh
        tmp=find(rat>=1);
        h=hist(rat(tmp),40)
        hist(rat(tmp),40)
        hold on
        tmp=find(rat<1);
        h=hist(rat(tmp),40)
        hist(rat(tmp),40)
        hh = findobj(gca,'Type','patch');
        set(hh(1),'FaceColor','r','EdgeColor','r')
        perc=100/length(rat)*length(tmp);
        title(['sum ',int2str(temps(ss)),'Hz (',int2str(perc),'% F2>F1)'])
        tmp=find(~isnan(rat));
        me=median(rat(tmp));
        legend(num2str(me))
        axis tight
    end
    
    
    figure(1)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath2,[Dtype,'_tuningCuves.bmp']))
    close
    
    figure(2)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath2,[Dtype,'_F1vsF2.bmp']))
    close
    
    figure(3)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(savepath2,[Dtype,'_F1F2histo.bmp']))
    close
end
%% 15d) find maximal F2/F1
% see "identification and characterization of a y-like primate retinal
% ganglion cell type" Petrusca et al. 2007

clear
close all
superpath='S:\data\Katja\allstuff';
load(fullfile(superpath,'DGinfo'))
savepath2=fullfile(superpath,'overview_pics');


types='mph';
cols='bgr';
leg={'F1','F2','F3','F4'};
mas=[5 5 3 3]; liw=[2 2 1 1];
titles={'mouse','pig','human'};

figure
for tt=1:length(types)
    
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_DGallpeaks']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    togo=find(goodDG==1);
    
    
    maxf2f1=[];
    for ss=1:length(togo)
        curr=cell2mat(DGallpeaks(togo(ss),:));
        ratios=curr(2,:)./curr(1,:);
        maxf2f1=[maxf2f1;max(rats)];
    end
    maxf2f1=[maxf2f1;-1;10];
    subplot(2,2,tt)
    [h x]=hist(maxf2f1,22);
    bar(x(2:end-1),h(2:end-1),cols(tt))
    
    
end
%% 15d) look at ON-OFF cells


clear
close all
superpath='S:\data\Katja\allstuff';


lf=9; flash=8;
types='mph';
cols='bgr';
titles={'mouse','pig','human'};


for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_DG_FFTdata']))
    load(fullfile(superpath,'DGinfo'))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    
    consider=find(goodones==1);
    flashdata=cell2mat(info(consider,flash));
    onoff=find(flashdata==2);
    onoff=consider(onoff);
    id=find(goodDG==1);
    on=find(flashdata==1);
    on=consider(on);
    
    figure
    subplot(2,2,1)
    data1=[]; data2=[];
    for o=1:length(onoff)
        idnow=find(id==onoff(o));
        f1=cell2mat(spatmax(idnow,1));
        f2=cell2mat(spatmax(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            plot(spats,f1,'--k')
            hold on
            data2=[data2;f2];
            plot(spats,f2,'--b')
            hold on
        end
    end
    plot(spats,median(data1),'-k','linewidth',4)
    plot(spats,median(data2),'-b','linewidth',4)
    
    subplot(2,2,2)
    data1=[]; data2=[];
    for o=1:length(on)
        idnow=find(id==on(o));
        f1=cell2mat(spatsum(idnow,1));
        f2=cell2mat(spatsum(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            plot(spats,f1,'--k')
            hold on
            data2=[data2;f2];
            plot(spats,f2,'--b')
            hold on
        end
    end
    plot(spats,median(data1),'-k','linewidth',4)
    plot(spats,median(data2),'-b','linewidth',4)
    
    subplot(2,2,3)
    data1=[]; data2=[];
    for o=1:length(onoff)
        idnow=find(id==onoff(o));
        f1=cell2mat(tempsum(idnow,1));
        f2=cell2mat(tempsum(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            plot(temps,f1,'--k')
            hold on
            data2=[data2;f2];
            plot(temps,f2,'--b')
            hold on
        end
    end
    plot(temps,median(data1),'-k','linewidth',4)
    plot(temps,median(data2),'-b','linewidth',4)
    
    
    subplot(2,2,4)
    data1=[]; data2=[];
    for o=1:length(on)
        idnow=find(id==on(o));
        f1=cell2mat(tempmax(idnow,1));
        f2=cell2mat(tempmax(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            plot(temps,f1,'--k')
            hold on
            data2=[data2;f2];
            plot(temps,f2,'--b')
            hold on
        end
    end
    plot(temps,median(data1),'-k','linewidth',4)
    plot(temps,median(data2),'-b','linewidth',4)
    
end
%% 15e) look at ON-OFF cells: sum


clear
close all
superpath='S:\data\Katja\allstuff';


lf=9; flash=8;
types='mph';
cols='bgr';
titles={'mouse','pig','human'};


for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_DG_FFTdata']))
    load(fullfile(superpath,'DGinfo'))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    
    consider=find(goodones==1);
    flashdata=cell2mat(info(consider,flash));
    onoff=find(flashdata==2);
    onoff=consider(onoff);
    id=find(goodDG==1);
    on=find(flashdata==1);
    on=consider(on);
    off=find(flashdata==-1);
    off=consider(off);
    
    figure(tt)
    subplot(2,5,1)
    data1=[]; data2=[];
    for o=1:length(onoff)
        idnow=find(id==onoff(o));
        f1=cell2mat(spatsum(idnow,1));
        f2=cell2mat(spatsum(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            data2=[data2;f2];
            plot(spats,f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(spats,median(rat(tmp,:)),'-b','linewidth',4)
    title('on-off')
    subplot(2,5,5)
    plot(spats,median(rat(tmp,:)),'-b','linewidth',4)
    hold on
    figure(100)
    subplot(2,3,tt)
    plot(spats,median(rat(tmp,:)),'-b','linewidth',4)
    hold on
    
    figure(tt)
    subplot(2,5,2)
    data1=[]; data2=[];
    for o=1:length(on)
        idnow=find(id==on(o));
        f1=cell2mat(spatsum(idnow,1));
        f2=cell2mat(spatsum(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            data2=[data2;f2];
            plot(spats,f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(spats,median(rat(tmp,:)),'-g','linewidth',4)
    title('on')
    subplot(2,5,5)
    plot(spats,median(rat(tmp,:)),'-g','linewidth',4)
    figure(100)
    subplot(2,3,tt)
    plot(spats,median(rat(tmp,:)),'-g','linewidth',4)
    
    figure(tt)
    subplot(2,5,3)
    data1=[]; data2=[];
    for o=1:length(off)
        idnow=find(id==off(o));
        f1=cell2mat(spatsum(idnow,1));
        f2=cell2mat(spatsum(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            data2=[data2;f2];
            plot(spats,f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(spats,median(rat(tmp,:)),'-c','linewidth',4)
    title('off')
    subplot(2,5,5)
    plot(spats,median(rat(tmp,:)),'-c','linewidth',4)
    figure(100)
    subplot(2,3,tt)
    plot(spats,median(rat(tmp,:)),'-c','linewidth',4)
    
    figure(tt)
    subplot(2,5,4)
    data1=[]; data2=[];
    for o=1:length(consider)
        idnow=find(id==consider(o));
        f1=cell2mat(spatsum(idnow,1));
        f2=cell2mat(spatsum(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            data2=[data2;f2];
            plot(spats,f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(spats,median(rat(tmp,:)),'-r','linewidth',4)
    title('all')
    subplot(2,5,5)
    plot(spats,median(rat(tmp,:)),'-r','linewidth',4)
    figure(100)
    subplot(2,3,tt)
    plot(spats,median(rat(tmp,:)),'-r','linewidth',4)
    
    figure(tt)
    subplot(2,5,6)
    data1=[]; data2=[];
    for o=1:length(onoff)
        idnow=find(id==onoff(o));
        f1=cell2mat(tempsum(idnow,1));
        if ~isempty(f1)
            f1=f1(1:3);
            f2=cell2mat(tempsum(idnow,2));
            f2=f2(1:3);
            
            data1=[data1;f1];
            data2=[data2;f2];
            plot(temps(1:3),f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(temps(1:3),median(rat(tmp,:)),'-b','linewidth',4)
    title('on-off')
    subplot(2,5,10)
    plot(temps(1:3),median(rat(tmp,:)),'-b','linewidth',4)
    hold on
    figure(100)
    subplot(2,3,tt+3)
    plot(temps(1:3),median(rat(tmp,:)),'-b','linewidth',4)
    hold on
    
    figure(tt)
    subplot(2,5,7)
    data1=[]; data2=[];
    for o=1:length(on)
        idnow=find(id==on(o));
        f1=cell2mat(tempsum(idnow,1));
        if ~isempty(f1)
            f1=f1(1:3);
            f2=cell2mat(tempsum(idnow,2));
            f2=f2(1:3);
            
            data1=[data1;f1];
            data2=[data2;f2];
            plot(temps(1:3),f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(temps(1:3),median(rat(tmp,:)),'-g','linewidth',4)
    title('on')
    subplot(2,5,10)
    plot(temps(1:3),median(rat(tmp,:)),'-g','linewidth',4)
    figure(100)
    subplot(2,3,tt+3)
    plot(temps(1:3),median(rat(tmp,:)),'-g','linewidth',4)
    
    figure(tt)
    subplot(2,5,8)
    data1=[]; data2=[];
    for o=1:length(off)
        idnow=find(id==off(o));
        f1=cell2mat(tempsum(idnow,1));
        if ~isempty(f1)
            f1=f1(1:3);
            f2=cell2mat(tempsum(idnow,2));
            f2=f2(1:3);
            
            data1=[data1;f1];
            data2=[data2;f2];
            plot(temps(1:3),f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(temps(1:3),median(rat(tmp,:)),'-c','linewidth',4)
    title('off')
    subplot(2,5,10)
    plot(temps(1:3),median(rat(tmp,:)),'-c','linewidth',4)
    figure(100)
    subplot(2,3,tt+3)
    plot(temps(1:3),median(rat(tmp,:)),'-c','linewidth',4)
    
    figure(tt)
    subplot(2,5,9)
    data1=[]; data2=[];
    for o=1:length(consider)
        idnow=find(id==consider(o));
        f1=cell2mat(tempsum(idnow,1));
        if ~isempty(f1)
            f1=f1(1:3);
            f2=cell2mat(tempsum(idnow,2));
            f2=f2(1:3);
            
            data1=[data1;f1];
            data2=[data2;f2];
            plot(temps(1:3),f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(temps(1:3),median(rat(tmp,:)),'-r','linewidth',4)
    title('all')
    subplot(2,5,10)
    plot(temps(1:3),median(rat(tmp,:)),'-r','linewidth',4)
    figure(100)
    subplot(2,3,tt+3)
    plot(temps(1:3),median(rat(tmp,:)),'-r','linewidth',4)
    
end
legend('on-off','on','off','all')
%% 15e) look at ON-OFF cells: max


clear
close all
superpath='S:\data\Katja\allstuff';


lf=9; flash=8;
types='mph';
cols='bgr';
titles={'mouse','pig','human'};


for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_DG_FFTdata']))
    load(fullfile(superpath,'DGinfo'))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    
    consider=find(goodones==1);
    flashdata=cell2mat(info(consider,flash));
    onoff=find(flashdata==2);
    onoff=consider(onoff);
    id=find(goodDG==1);
    on=find(flashdata==1);
    on=consider(on);
    off=find(flashdata==-1);
    off=consider(off);
    
    figure(tt)
    subplot(2,5,1)
    data1=[]; data2=[];
    for o=1:length(onoff)
        idnow=find(id==onoff(o));
        f1=cell2mat(spatmax(idnow,1));
        f2=cell2mat(spatmax(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            data2=[data2;f2];
            plot(spats,f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(spats,median(rat(tmp,:)),'-b','linewidth',4)
    title('on-off')
    subplot(2,5,5)
    plot(spats,median(rat(tmp,:)),'-b','linewidth',4)
    hold on
    figure(100)
    subplot(2,3,tt)
    plot(spats,median(rat(tmp,:)),'-b','linewidth',4)
    hold on
    
    figure(tt)
    subplot(2,5,2)
    data1=[]; data2=[];
    for o=1:length(on)
        idnow=find(id==on(o));
        f1=cell2mat(spatmax(idnow,1));
        f2=cell2mat(spatmax(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            data2=[data2;f2];
            plot(spats,f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(spats,median(rat(tmp,:)),'-g','linewidth',4)
    title('on')
    subplot(2,5,5)
    plot(spats,median(rat(tmp,:)),'-g','linewidth',4)
    figure(100)
    subplot(2,3,tt)
    plot(spats,median(rat(tmp,:)),'-g','linewidth',4)
    
    figure(tt)
    subplot(2,5,3)
    data1=[]; data2=[];
    for o=1:length(off)
        idnow=find(id==off(o));
        f1=cell2mat(spatmax(idnow,1));
        f2=cell2mat(spatmax(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            data2=[data2;f2];
            plot(spats,f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(spats,median(rat(tmp,:)),'-c','linewidth',4)
    title('off')
    subplot(2,5,5)
    plot(spats,median(rat(tmp,:)),'-c','linewidth',4)
    figure(100)
    subplot(2,3,tt)
    plot(spats,median(rat(tmp,:)),'-c','linewidth',4)
    
    figure(tt)
    subplot(2,5,4)
    data1=[]; data2=[];
    for o=1:length(consider)
        idnow=find(id==consider(o));
        f1=cell2mat(spatmax(idnow,1));
        f2=cell2mat(spatmax(idnow,2));
        if ~isempty(f1)
            data1=[data1;f1];
            data2=[data2;f2];
            plot(spats,f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(spats,median(rat(tmp,:)),'-r','linewidth',4)
    title('all')
    subplot(2,5,5)
    plot(spats,median(rat(tmp,:)),'-r','linewidth',4)
    figure(100)
    subplot(2,3,tt)
    plot(spats,median(rat(tmp,:)),'-r','linewidth',4)
    
    figure(tt)
    subplot(2,5,6)
    data1=[]; data2=[];
    for o=1:length(onoff)
        idnow=find(id==onoff(o));
        f1=cell2mat(tempmax(idnow,1));
        if ~isempty(f1)
            f1=f1(1:3);
            f2=cell2mat(tempmax(idnow,2));
            f2=f2(1:3);
            
            data1=[data1;f1];
            data2=[data2;f2];
            plot(temps(1:3),f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(temps(1:3),median(rat(tmp,:)),'-b','linewidth',4)
    title('on-off')
    subplot(2,5,10)
    plot(temps(1:3),median(rat(tmp,:)),'-b','linewidth',4)
    hold on
    figure(100)
    subplot(2,3,tt+3)
    plot(temps(1:3),median(rat(tmp,:)),'-b','linewidth',4)
    hold on
    
    figure(tt)
    subplot(2,5,7)
    data1=[]; data2=[];
    for o=1:length(on)
        idnow=find(id==on(o));
        f1=cell2mat(tempmax(idnow,1));
        if ~isempty(f1)
            f1=f1(1:3);
            f2=cell2mat(tempmax(idnow,2));
            f2=f2(1:3);
            
            data1=[data1;f1];
            data2=[data2;f2];
            plot(temps(1:3),f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(temps(1:3),median(rat(tmp,:)),'-g','linewidth',4)
    title('on')
    subplot(2,5,10)
    plot(temps(1:3),median(rat(tmp,:)),'-g','linewidth',4)
    figure(100)
    subplot(2,3,tt+3)
    plot(temps(1:3),median(rat(tmp,:)),'-g','linewidth',4)
    
    figure(tt)
    subplot(2,5,8)
    data1=[]; data2=[];
    for o=1:length(off)
        idnow=find(id==off(o));
        f1=cell2mat(tempmax(idnow,1));
        if ~isempty(f1)
            f1=f1(1:3);
            f2=cell2mat(tempmax(idnow,2));
            f2=f2(1:3);
            
            data1=[data1;f1];
            data2=[data2;f2];
            plot(temps(1:3),f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(temps(1:3),median(rat(tmp,:)),'-c','linewidth',4)
    title('off')
    subplot(2,5,10)
    plot(temps(1:3),median(rat(tmp,:)),'-c','linewidth',4)
    figure(100)
    subplot(2,3,tt+3)
    plot(temps(1:3),median(rat(tmp,:)),'-c','linewidth',4)
    
    figure(tt)
    subplot(2,5,9)
    data1=[]; data2=[];
    for o=1:length(consider)
        idnow=find(id==consider(o));
        f1=cell2mat(tempmax(idnow,1));
        if ~isempty(f1)
            f1=f1(1:3);
            f2=cell2mat(tempmax(idnow,2));
            f2=f2(1:3);
            
            data1=[data1;f1];
            data2=[data2;f2];
            plot(temps(1:3),f1./f2,'--k')
            hold on
        end
    end
    ratio=data1./data2;
    tmp=[];
    for r=1:size(rat,1)
        if sum(~isnan(rat(r,:)))==size(rat,2);
            tmp=[tmp;r];
        end
    end
    plot(temps(1:3),median(rat(tmp,:)),'-r','linewidth',4)
    title('all')
    subplot(2,5,10)
    plot(temps(1:3),median(rat(tmp,:)),'-r','linewidth',4)
    figure(100)
    subplot(2,3,tt+3)
    plot(temps(1:3),median(rat(tmp,:)),'-r','linewidth',4)
    
end
legend('on-off','on','off','all')
%% 15f) look at full-field response of cells with low F1/F2 for narrow stimuli
clear
close all
superpath='S:\data\Katja\allstuff';
load(fullfile(superpath,'DGinfo'))
savepath2=fullfile(superpath,'overview_pics');


types='mph';
cols='bgr';
titles={'mouse','pig','human'};

figure
for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_DG_FFTdata']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_flash_new']))
    load(fullfile(superpath,'LFs',[Dtype,'_LFs']))
    
    id=find(goodDG==1);
    
    f1=cell2mat(spatmax(:,1));
    f2=cell2mat(spatmax(:,2));
    ratio=f1./f2;
    
    cand=[];
    for i=1:length(rat)
        if ratio(i,1)<1 || ratio(i,2)<1
            %         if ratio(i,3)<1
            cand=[cand id(i)];
        end
    end
    
    
    tmp2=[];
    for c=1:length(cand)
        if ~isempty(find(togo==cand(c)))
            tmp2=[tmp2;cand(c)];
        end
    end
    tmp3=find(goodones(tmp2)==1);
    su=sqrt(length(tmp3)); su=ceil(su);
    
    
    coll=[]; cnt=0;
    for c=1:length(cand)
        if goodones(cand(c))==1 && ~isempty(find(togo==cand(c)))
            cnt=cnt+1;
            id2=find(togo==cand(c));
            
            data=flash_rates{id2};
            data=data/max(data);
            coll=[coll;data];
            figure(1)
            subplot(2,2,tt)
            plot(data,'-','color',cols(tt))
            hold on
            figure(100+tt)
            subplot(su,su,cnt)
            plot(data,'-','color',cols(tt))
            hold on
            
        end
    end
    figure(1)
    subplot(2,2,tt)
    plot(median(coll),'-k','linewidth',2)
    
    
    tmp2=[];
    for c=1:length(cand)
        if ~isempty(find(LFs(:,cand(c))>0))
            tmp2=[tmp2;cand(c)];
        end
    end
    tmp3=find(goodones(tmp2)==1);
    tmp3=tmp2(tmp3);
    su=sqrt(length(tmp3)); su=ceil(su);
    
    cnt=0;
    for c=1:length(cand)
        if goodones(cand(c))==1 && ~isempty(find(LFs(:,cand(c))>0))
            cnt=cnt+1;
            data=LFs(:,cand(c));
            figure(200+tt)
            subplot(su,su,cnt)
            plot(data,'-','color',cols(tt))
            hold on
        end
    end
    
end
%% 15g) look at chirp response of cells with low F1/F2 for narrow stimuli
clear
close all
superpath='S:\data\Katja\allstuff';
load(fullfile(superpath,'DGinfo'))
savepath2=fullfile(superpath,'overview_pics');


types='ph';
cols='gr';
titles={'mouse','pig','human'};



for tt=1:length(types)
    Dtype=types(tt);
    load(fullfile(superpath,[Dtype,'_DG_FFTdata']))
    load(fullfile(superpath,[Dtype,'_DGallpeaks']))
    load(fullfile(superpath,[Dtype,'_goodDG']))
    load(fullfile(superpath,[Dtype,'_info']))
    load(fullfile(superpath,[Dtype,'_chirpFFT']))
    load(fullfile(superpath,'DGinfo'))
    spa=cell2mat(stimInfo(1,2));
    spa=find(spa<500);
    
    id=find(goodDG==1);
    
    cand=[];
    for c=1:length(id)
        curr=cell2mat(DGallpeaks(id(c),:));
        ratio=curr(1,:)./curr(2,:);
        ratio=ratio(spa);
        if ~isempty(find(ratio<1))
            cand=[cand;id(c)];
        end
    end
    
    tmp2=[];
    for c=1:length(cand)
        if ~isempty(find(togo==cand(c)))
            tmp2=[tmp2;cand(c)];
        end
    end
    su=sqrt(length(tmp2)); su=ceil(su);
    
    noncand=setdiff(1:length(info),cand);
    coll=[]; cnt=0; coll2=[];
    for c=1:length(cand)
        if ~isempty(find(togo==cand(c)))
            cnt=cnt+1;
            id2=find(togo==cand(c));
            data=chirpFFTsmooth{id2,1};
            
            coll=[coll data];
            figure(1)
            subplot(2,3,(tt*3)-2)
            plot(data,'-','color',cols(tt))
            hold on
            figure(100+tt)
            subplot(su,su,cnt)
            plot(data,'-','color',cols(tt))
            hold on
        end
    end
    figure(1)
    subplot(2,3,(tt*3)-2)
    plot(median(coll,2),'-k','linewidth',2)
    subplot(2,3,(tt*3))
    plot(median(coll,2),'-k','linewidth',2)
    hold on
    
    
    for c=1:length(noncand)
        if ~isempty(find(togo==noncand(c)))
            id2=find(togo==noncand(c));
            data=chirpFFTsmooth{id2,1};
            
            coll2=[coll2 data];
            figure(1)
            subplot(2,3,(tt*3)-1)
            plot(data,'-','color',cols(tt))
            hold on
            
        end
    end
    figure(1)
    subplot(2,3,(tt*3)-1)
    plot(median(coll2,2),'--k','linewidth',2)
    subplot(2,3,(tt*3))
    plot(median(coll2,2),'--k','linewidth',2)
    
end
%% --16a) DG FFT further analysis
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
savepath=fullfile(superpath,'DG_FFT');
savepath2=fullfile(superpath,'DG_allFFT');
plotting1=1; plotting2=1;

forplot=[9 13 14 17 18 21];
togo=find(goodDG==1);
DGfourier=zeros(length(togo),24,8900);%fourier result
DGfourierF=zeros(length(togo),24,8900); %fourier frequency
DGfourierCorr=zeros(length(togo),24,8900);%corrected fourier result
DGpeaks=zeros(length(togo),24,3);%peaks of first 3 harmonics
DGpeaksCorr=zeros(length(togo),24,3);%peaks of first 3 harmonics of corrected fourier
DGphase=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics
DGphasecorr=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics of corrected fourier

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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    thiscell=[];
    fouriers=[];
    rates=[];
    clear f2take
    dgcnt=0;
    for g=1:24
        if ~isempty(find(g==forplot))
            dgcnt=dgcnt+1;
        end
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            
            stimFreq=freqs(g);
            
            if max(meanrate)>0
                currfr=stimInfo{currID,1}(g);
                currspat=stimInfo{currID,2}(g);
                
                meanF=mean(meanrate);
                meanFR=round(meanF*10)/10;
                % fourier
                x=meanrate(2000:14000);
                if max(x)>0
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt)
                        plot(x)
                        axis tight
                        subplot(5,6,dgcnt+6)
                        plot(x(2000:5000))
                        axis tight
                        title([int2str(currspat),'/',int2str(currfr)])
                    end
                    
                    fs=1000;
                    n=2^nextpow2(numel(x));
                    y=fft(x,n)/numel(x);
                    f=fs/2*linspace(0,1,n/2+1);
                    y=y(1:n/2+1);
                    c=2*abs(y);
                    DGfourier(i,g,1:length(c))=c;
                    DGfourierF(i,g,1:length(f))=f;
                    
                    tmp=find(f>0.6); tmp=tmp(1);
                    tmp2=find(f<17); tmp2=tmp2(end);
                    tmp3=find(f>0.2);tmp3=tmp3(1);
                    tmp4=find(f>0.6);tmp4=tmp4(1);
                    ma=max(c(tmp:tmp2));
                    
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt+12)
                        plot(f,c,'-g')
                        axis([0 17 0 ma])
                    end
                    if plotting2==1
                        figure(100)
                        subplot(4,6,g)
                        plot(f,c,'-g')
                        axis([0 17 0 ma])
                        title([int2str(currspat),'/',int2str(currfr)])
                    end
                    
                    fi=fit(f(tmp3:end)',c(tmp3:end)','exp1');
                    aa=fi.a; bb=fi.b;
                    fi=aa*exp(bb*f(tmp3:end));
                    new=c(tmp3:end)-fi;
                    no=find(new<0);
                    new(no)=0;
                    st=std(new);
                    keepc=c;
                    %                     c=new-st;
                    c=new;
                    DGfourierCorr(i,g,1:length(c))=c;
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt+12)
                        hold on
                        plot(f(tmp3:end),fi)
                        plot(f(tmp3:end),c,'-m','linewidth',1)
                    end
                    if plotting2==1
                        figure(100)
                        subplot(4,6,g)
                        hold on
                        plot(f(tmp3:end),fi)
                        plot(f(tmp3:end),c,'-m','linewidth',1)
                    end
                    
                    
                    y0=fftshift(y);
                    f0=linspace(0,2,n/2+1)*(fs/4);
                    f0=f0-max(f0)/2;
                    phase = unwrap(angle(y0));
                    phases=phase*180/pi;
                    
                    tmp=find(f>currfr-0.3); tmp=tmp(1);
                    tmp2=find(f<33); tmp2=tmp2(end);
                    pks=[];locs=[]; clocs=[];
                    cc=keepc; realcurr=currfr;
                    for l=1:4
                        tmp=find(f>realcurr*l);
                        [test1 test2]=max(cc(tmp(1)-4:tmp(1)+3));
                        cc(tmp(1)-4+test2-3:tmp(1)-4+test2+3)=0;
                        %                         if l==1
                        %                             realcurr=f(test2+tmp(1)-7);
                        %                         end
                        locs=[locs test2+tmp(1)-7];
                        clocs=[clocs tmp(1)-4+test2];
                        pks=[pks test1];
                        
                    end
                    if locs(1)-10<1
                        beg=1;
                    else
                        beg=locs(1)-10;
                    end
                    restc=cc(locs(1)-10:locs(end)+10);
                    clocs=clocs-beg+1;
                    
                    restc=cc(locs(1)-10:locs(end)+10);
                    %                     cc=keepc;
                    %                     wrongpks=[];
                    %                     for l=1.5:4.5
                    %                       tmp=find(f>realcurr*l);
                    %                         [test1 test2]=max(cc(tmp(1)-4:tmp(1)+2));
                    %                         wrongpks=[wrongpks test1];
                    %                     end
                    for pp=1:3
                        if  pks(pp)<max(restc(clocs(pp:pp+1)-1))
                            pks(pp)=0;
                        end
                    end
                    if max(pks)<max(restc)+std(restc)
                        pks(1:4)=0;
                    end
                    pks=pks-std(restc);
                    DGpeaks(i,g,1:3)=pks(1:3);
                    frqs=f(locs);
                    t=0:0.001:12;
                    sines=zeros(12001,4);
                    phs=[];
                    cnt=0;
                    for fr=frqs
                        cnt=cnt+1;
                        tmp=find(f0==fr);
                        phs=[phs;phases(tmp)];
                        s=pks(cnt)*(4/pi)*sin(2*pi*fr*t+phs(cnt));
                        sines(:,cnt)=s;
                    end
                    di=frqs-currfr;
                    [mi id]=min(abs(di));
                    ph=phs(id);
                    di=frqs-2*currfr;
                    [mi id]=min(abs(di));
                    ph=[ph phs(id)];
                    di=frqs-3*currfr;
                    [mi id]=min(abs(di));
                    ph=[ph phs(id)];
                    currval=round((ph(1)-ph(2))/360*100)/100;
                    currval2=round((ph(1)-ph(3))/360*100)/100;
                    DGphase(i,g,1:2)=[currval currval2];
                    
                    
                    
                    pks=[];locs=[];clocs=[];
                    cc=c; realcurr=currfr;
                    for l=1:4
                        tmp=find(f>realcurr*l);
                        [test1 test2]=max(cc(tmp(1)-tmp3-3:tmp(1)+4-tmp3));%-4, +3 positions
                        cc(tmp(1)-tmp3-3+test2-1-3:tmp(1)-tmp3-3+test2-1+3)=0;
                        %                         if l==1
                        %                             realcurr=f(test2+tmp(1)-3);
                        %                         end
                        pks=[pks test1];
                        locs=[locs test2+tmp(1)-5];
                        clocs=[clocs tmp(1)-tmp3-3+test2-1];
                        
                    end
                    if locs(1)-tmp3-9<1
                        beg=1;
                    else
                        beg=locs(1)-tmp3-9;
                    end
                    restc=cc(beg:locs(end)-tmp3+11);
                    clocs=clocs-beg+1;
                    %                     cc=c;
                    %                     wrongpks=[];
                    %                    for l=1.5:4.5
                    %                       tmp=find(f>realcurr*l);
                    %                         [test1 test2]=max(cc(tmp(1)-tmp3-1:tmp(1)+5-tmp3));
                    %                         wrongpks=[wrongpks test1];
                    %                     end
                    for pp=1:3
                        if pks(pp)<max(restc(clocs(pp:pp+1)-1))
                            pks(pp)=0;
                        end
                    end
                    testing=find(pks<mean(restc)+2*std(restc));
                    pks(testing)=0;
                    if max(pks)<max(restc)+std(restc)
                        pks(1:4)=0;
                    end
                    
                    pks=pks-std(restc);
                    DGpeaksCorr(i,g,1:3)=pks(1:3);
                    
                    
                    frqs=f(locs);
                    
                    
                    t=0:0.001:12;
                    sines=zeros(12001,4);
                    phs=[];
                    cnt=0;
                    for fr=frqs
                        cnt=cnt+1;
                        tmp=find(f0==fr);
                        phs=[phs;phases(tmp)];
                        s=pks(cnt)*(4/pi)*sin(2*pi*fr*t+phs(cnt));
                        sines(:,cnt)=s;
                    end
                    
                    di=frqs-currfr;
                    [mi id]=min(abs(di));
                    ph=phs(id);
                    di=frqs-2*currfr;
                    [mi id]=min(abs(di));
                    ph=[ph phs(id)];
                    di=frqs-3*currfr;
                    [mi id]=min(abs(di));
                    ph=[ph phs(id)];
                    
                    currval=round((ph(1)-ph(2))/360*100)/100;
                    currval2=round((ph(1)-ph(3))/360*100)/100;
                    DGphaseCorr(i,g,1:2)=[currval currval2];
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt+12)
                        title(['FFT ',num2str(currval),'/',num2str(currval2),' periods'])
                        
                        
                        for s=1:3
                            subplot(5,6,dgcnt+18)
                            plot(sines(:,s),'-r')
                            hold on
                            subplot(5,6,dgcnt+24)
                            plot(sines(1:3000,s),'-r')
                            hold on
                        end
                        subplot(5,6,dgcnt+18)
                        plot(sum(sines'),'-k','linewidth',2)
                        axis tight
                        axis off
                        if dgcnt==1
                            title('reconstruction')
                        elseif dgcnt==2
                            title('with 3 strongest sines')
                        elseif dgcnt==4
                            title('si=peak*(4/pi)*sin(2pi*f*t+phase(f))')
                        end
                        
                        subplot(5,6,dgcnt+24)
                        plot(sum(sines(1:3000,:)'),'-k','linewidth',2)
                        axis tight
                    end
                    if plotting2==1
                        figure(100)
                        subplot(4,6,g)
                        for pp=1:3
                            if pks(pp)>0
                                plot(frqs(pp),pks(pp),'*k')
                            end
                        end
                    end
                    
                    %
                    %                                 sums=sum(sines(:,:)');
                    % fs=1000;
                    %                 n=2^nextpow2(numel(sums));
                    %                 y=fft(sums,n);
                    %                 f=(0:n-1)*(fs/n);
                    %                 c=abs(y);
                    
                    
                end
            end
        end
        
    end
    if plotting1==1
        figure(1)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_FFT.bmp']))
        close
    end
    if plotting2==1
        figure(100)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath2,[Dtype,'_',int2str(currID),'_allFFT.bmp']))
        close
    end
end
save(fullfile(superpath,[Dtype,'_DG_extendedFFT']),'DGfourier','DGfourierCorr','DGpeaks','DGpeaksCorr','DGphase','DGphaseCorr','DGfourierF')
%% 17) FFT of spont activity

load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)
spontFFT=zeros(size(info,1),100);

for i=1:length(info)
    
    currd=info{i,1};
    curru=info{i,2};
    currstart=info{i,4};
    currstop=info{i,5};
    
    %     hpath=fullfile(mainpath,currd,'HEKA');
    %     hlist=dir(fullfile(hpath,'*phys'));
    upath=fullfile(mainpath,currd,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,curru))
            found=1;
        end
    end
    load(fullfile(upath,ulist(m).name))
    spont=[];
    for f=currstart:currstop
        if ~isempty(regexp(unit{1,2}{f,1},'spont'))
            spont=[spont;f];
        end
    end
    
    if ~isempty(spont)
        allrate=[];
        for f=1:length(spont)
            spikes=unit{1,2}{spont(f),2};
            tmp=find(spikes<12000);
            spikes=spikes(tmp);
            rate=MEA_spikerates(spikes,40,12000);
            allrate=[allrate;rate];
        end
        if size(allrate,1)>1
            meanrate=mean(allrate);
        else
            meanrate=allrate;
        end
        fs=1000;
        n=6*fs;
        y=fft(meanrate,n)/numel(meanrate);
        f=fs/2*linspace(0,1,n/2+1);
        y=y(1:n/2+1);
        c=2*abs(y);
        spontFFT(i,1:length(c))=c;
        
        
    end
end
save(fullfile(superpath,[Dtype,'_spontFFT']),'spontFFT')
%% (-) 16a) DG FFT further analysis
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
load(fullfile(superpath,[Dtype,'_DGrates']))
cd(superpath)
savepath=fullfile(superpath,'DG_newFFTcalc');
savepath2=fullfile(superpath,'DG_FFTreconstr');
plotting1=1; plotting2=1;

forplot=[9 13 14 17 18 21];
togo=find(goodDG==1);
DGfourier=zeros(length(togo),24,8900);%fourier result
DGfourierF=zeros(length(togo),24,8900); %fourier frequency
DGfourierCorr=zeros(length(togo),24,8900);%corrected fourier result
DGpeaks=zeros(length(togo),24,3);%peaks of first 3 harmonics
DGpeaksCorr=zeros(length(togo),24,3);%peaks of first 3 harmonics of corrected fourier
DGphase=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics
DGphasecorr=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics of corrected fourier

for i=1:length(togo)
    i
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    
    
    
    thiscell=[];
    fouriers=[];
    clear f2take
    dgcnt=0;
    for g=1:24
        if ~isempty(find(g==forplot))
            dgcnt=dgcnt+1;
        end
        meanrate=DGrates{currID,g};
        x=meanrate(2000:14000);
        if ~isempty(x)
            
            stimFreq=stimInfo{currID,1}(g);
            
            if max(x)>0
                currfr=stimInfo{currID,1}(g);
                currspat=stimInfo{currID,2}(g);
                
                meanF=mean(x);
                meanFR=round(meanF*10)/10;
                % fourier
                
                if ~isempty(find(g==forplot)) && plotting1==1
                    figure(1)
                    subplot(5,6,dgcnt)
                    plot(x)
                    axis tight
                    subplot(5,6,dgcnt+6)
                    plot(x(2000:5000))
                    axis tight
                    title([int2str(currspat),'/',int2str(currfr)])
                end
                
                fs=1000; %sampling rate
                n=4*fs;
                y=fft(x,n)/numel(x);
                f=fs/2*linspace(0,1,n/2+1); %frequencies
                y=y(1:n/2+1);
                c=2*abs(y); %amplitudes
                
                
                
                
                
                
                DGfourier(i,g,1:length(c))=c;
                DGfourierF(i,g,1:length(f))=f;
                
                tmp=find(f>0.6); tmp=tmp(1);
                tmp2=find(f<17); tmp2=tmp2(end);
                tmp3=find(f>0.2);tmp3=tmp3(1);
                tmp4=find(f>0.6);tmp4=tmp4(1);
                ma=max(c(tmp:tmp2));
                
                if ~isempty(find(g==forplot)) && plotting1==1
                    figure(1)
                    subplot(5,6,dgcnt+12)
                    plot(f,c,'-g')
                    axis([0 17 0 ma])
                end
                if plotting2==1
                    figure(100)
                    subplot(4,6,g)
                    plot(f,c,'-g')
                    axis([0 17 0 ma])
                    title([int2str(currspat),'/',int2str(currfr)])
                end
                
                fi=fit(f(tmp3:end)',c(tmp3:end)','exp1');
                aa=fi.a; bb=fi.b;
                fi=aa*exp(bb*f(tmp3:end));
                new=c(tmp3:end)-fi;
                no=find(new<0);
                new(no)=0;
                st=std(new);
                keepc=c;
                c=[repmat(0,1,tmp3-1) new];
                
                DGfourierCorr(i,g,1:length(c))=c;
                if ~isempty(find(g==forplot)) && plotting1==1
                    figure(1)
                    subplot(5,6,dgcnt+12)
                    hold on
                    plot(f(tmp3:end),fi)
                    plot(f,c,'-m','linewidth',1)
                end
                if plotting2==1
                    figure(100)
                    subplot(4,6,g)
                    hold on
                    plot(f(tmp3:end),fi)
                    plot(f,c,'-m','linewidth',1)
                end
                
                
                %---extract phases---
                y0=fftshift(y);
                f0=linspace(0,2,n/2+1)*(fs/4);
                f0=f0-max(f0)/2;
                phases = angle(y0);
                
                tmp=find(f>currfr-0.3); tmp=tmp(1);
                tmp2=find(f<33); tmp2=tmp2(end);
                pks=[];locs=[];corrpks=[];
                restc=c;
                for l=1:4
                    tmp=find(f==currfr*l);
                    locs=[locs tmp];
                    pks=[pks keepc(tmp)];
                    corrpks=[corrpks c(tmp)];
                    restc(tmp-2:tmp+2)=0;
                end
                
                for pp=1:3
                    if corrpks(pp)<c(locs(pp)-1) || corrpks(pp)<c(locs(pp)+1)
                        corrpks(pp)=0;
                        pks(pp)=0;
                    end
                    if  corrpks(pp)<max(restc(tmp3:tmp2))*1.5
                        corrpks(pp)=0;
                        pks(pp)=0;
                    end
                end
                
                DGpeaks(i,g,1:3)=pks(1:3);
                DGpeaksCorr(i,g,1:3)=corrpks(1:3);
                
                t=0:0.001:12;
                cosines=zeros(12001,4);
                phs=[];
                cnt=0;
                for fr=currfr*(1:3)
                    cnt=cnt+1;
                    tmp=find(f0==fr);
                    phs=[phs;phases(tmp)];
                    s=pks(cnt)*(4/pi)*cos(2*pi*fr*t+phs(cnt));
                    cosines(:,cnt)=s;
                end
                
                currval=round((ph(1)-ph(2))/360*100)/100;
                currval2=round((ph(1)-ph(3))/360*100)/100;
                DGphase(i,g,1:2)=phs(1:2);
                
                if ~isempty(find(g==forplot)) && plotting1==1
                    figure(1)
                    subplot(5,6,dgcnt+12)
                    title(['FFT | phase ratio: ',num2str(phs(2)/phs(1))])
                    
                    
                    for s=1:3
                        subplot(5,6,dgcnt+18)
                        plot(cosines(:,s),'-r')
                        hold on
                        subplot(5,6,dgcnt+24)
                        plot(cosines(1:3000,s),'-r')
                        hold on
                    end
                    subplot(5,6,dgcnt+18)
                    plot(sum(cosines'),'-k','linewidth',2)
                    axis tight
                    axis off
                    if dgcnt==1
                        title('reconstruction')
                    elseif dgcnt==2
                        title('with 3 first harmonics')
                    elseif dgcnt==4
                        title('si=peak*(4/pi)*cos(2pi*f*t+phase(f))')
                    end
                    
                    subplot(5,6,dgcnt+24)
                    plot(sum(cosines(1:3000,:)'),'-k','linewidth',2)
                    axis tight
                end
                if plotting2==1
                    figure(100)
                    subplot(4,6,g)
                    for pp=1:3
                        if corrpks(pp)>0
                            plot(frqs(pp),corrpks(pp),'*k')
                        end
                    end
                end
                
                
                
            end
        end
    end
    
    
    if plotting1==1
        figure(1)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath2,[Dtype,'_',int2str(currID),'_FFT.bmp']))
        close
    end
    if plotting2==1
        figure(100)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_allFFT.bmp']))
        close
    end
end
save(fullfile(superpath,[Dtype,'_DG_newFFTanalysis']),'DGfourier','DGfourierCorr','DGpeaks','DGpeaksCorr','DGphase','DGphaseCorr','DGfourierF')
save(fullfile(superpath,[Dtype,'_DG_newFFTsimple']),'DGpeaks','DGpeaksCorr','DGphase','DGfourierF')

%% (-) 16b) DG FFT further analysis
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
savepath=fullfile(superpath,'DG_FFT');
savepath2=fullfile(superpath,'DG_allFFT');
savepath3=fullfile(superpath,'DG_newFFT');
plotting1=0; plotting2=1;

forplot=[9 13 14 17 18 21];
togo=find(goodDG==1);
DGfourier=zeros(length(togo),24,8900);%fourier result
DGfourierF=zeros(length(togo),24,8900); %fourier frequency
DGfourierCorr=zeros(length(togo),24,8900);%corrected fourier result
DGpeaks=zeros(length(togo),24,3);%peaks of first 3 harmonics
DGpeaksCorr=zeros(length(togo),24,3);%peaks of first 3 harmonics of corrected fourier
DGphase=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics
DGphaseCorr=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics of corrected fourier

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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    thiscell=[];
    fouriers=[];
    rates=[];
    clear f2take
    dgcnt=0;
    for g=1:24
        if ~isempty(find(g==forplot))
            dgcnt=dgcnt+1;
        end
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            
            stimFreq=freqs(g);
            
            if max(meanrate)>0
                currfr=stimInfo{currID,1}(g);
                currspat=stimInfo{currID,2}(g);
                
                meanF=mean(meanrate);
                meanFR=round(meanF*10)/10;
                % fourier
                x=meanrate(2000:14000);
                if max(x)>0
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt)
                        plot(x)
                        axis tight
                        subplot(5,6,dgcnt+6)
                        plot(x(2000:5000))
                        axis tight
                        title([int2str(currspat),'/',int2str(currfr)])
                    end
                    
                    fs=1000;
                    n=2^nextpow2(numel(x));
                    y=fft(x,n)/numel(x);
                    f=fs/2*linspace(0,1,n/2+1);
                    y=y(1:n/2+1);
                    c=2*abs(y); keepc=c;
                    DGfourier(i,g,1:length(c))=c;
                    DGfourierF(i,g,1:length(f))=f;
                    
                    tmp=find(f>0.6); tmp=tmp(1);
                    tmp2=find(f<17); tmp2=tmp2(end);
                    tmp3=find(f>0.2);tmp3=tmp3(1);
                    tmp4=find(f>0.6);tmp4=tmp4(1);
                    
                    ma=max(c(tmp:tmp2));
                    
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt+12)
                        plot(f,c,'-g')
                        axis([0 17 0 ma])
                    end
                    if plotting2==1
                        figure(100)
                        subplot(4,6,g)
                        plot(f,c,'-g')
                        axis([0 17 0 ma])
                        title([int2str(currspat),'/',int2str(currfr)])
                    end
                    
                    
                    y0=fftshift(y);
                    f0=linspace(0,2,n/2+1)*(fs/4);
                    f0=f0-max(f0)/2;
                    phase = unwrap(angle(y0));
                    phases=phase*180/pi;
                    
                    tmp=find(f>currfr-0.3); tmp=tmp(1);
                    tmp2=find(f<33); tmp2=tmp2(end);
                    tmp5=find(f<17); tmp5=tmp5(end);
                    pks=[];locs=[]; clocs=[];
                    cc=keepc; realcurr=currfr;
                    for l=1:4
                        tmp=find(f>realcurr*l);
                        nowID=tmp(1);
                        [test1 test2]=max(cc(nowID-4:nowID+3));% 1 ~ nowID-4
                        cc(nowID-4+test2-3: nowID-4+test2+3)=0;
                        locs=[locs test2+nowID-4-1];
                        clocs=[clocs test2+nowID-4-1];
                        pks=[pks test1];
                        
                    end
                    if clocs(1)-10<1
                        beg=1;
                    elseif currfr>1
                        beg=clocs(1)-20;
                    else
                        beg=clocs(1)-10;
                    end
                    restc=cc(beg:clocs(end)+10);
                    clocs=clocs-beg+1;
                    
                    for pp=1:3
                        if  pks(pp)<max(restc(clocs(pp)-10:clocs(pp+1)-1))*1.5
                            pks(pp)=0;
                        end
                    end
                    if max(pks)<max(restc)+std(restc)
                        pks(1:4)=0;
                    end
                    pks=pks-std(restc);
                    
                    DGpeaks(i,g,1:3)=pks(1:3);
                    frqs=f(locs);
                    t=0:0.001:12;
                    sines=zeros(12001,4);
                    phs=[];
                    cnt=0;
                    for fr=frqs
                        cnt=cnt+1;
                        tmp=find(f0==fr);
                        phs=[phs;phases(tmp)];
                        s=pks(cnt)*(4/pi)*sin(2*pi*fr*t+phs(cnt));
                        sines(:,cnt)=s;
                    end
                    di=frqs-currfr;
                    [mi id]=min(abs(di));
                    ph=phs(id);
                    di=frqs-2*currfr;
                    [mi id]=min(abs(di));
                    ph=[ph phs(id)];
                    di=frqs-3*currfr;
                    [mi id]=min(abs(di));
                    ph=[ph phs(id)];
                    currval=round((ph(1)-ph(2))/360*100)/100;
                    currval2=round((ph(1)-ph(3))/360*100)/100;
                    DGphase(i,g,1:2)=[currval currval2];
                    
                    
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt+12)
                        title(['FFT ',num2str(currval),'/',num2str(currval2),' periods'])
                        
                        
                        for s=1:3
                            subplot(5,6,dgcnt+18)
                            plot(sines(:,s),'-r')
                            hold on
                            subplot(5,6,dgcnt+24)
                            plot(sines(1:3000,s),'-r')
                            hold on
                        end
                        subplot(5,6,dgcnt+18)
                        plot(sum(sines'),'-k','linewidth',2)
                        axis tight
                        axis off
                        if dgcnt==1
                            title('reconstruction')
                        elseif dgcnt==2
                            title('with sines of first 3 harmonics')
                        elseif dgcnt==4
                            title('si=peak*(4/pi)*sin(2pi*f*t+phase(f))')
                        end
                        
                        subplot(5,6,dgcnt+24)
                        plot(sum(sines(1:3000,:)'),'-k','linewidth',2)
                        axis tight
                    end
                    
                    
                    
                    
                    
                    
                    pks=[];locs=[];clocs=[];
                    cc=c; realcurr=currfr;
                    ccc=c;
                    numOfSteps=32/realcurr;
                    for l=1:numOfSteps
                        tmp=find(f>realcurr*l);
                        nowID=tmp(1);
                        [test1 test2]=max(cc(nowID-4:nowID+3));%-4, +3 positions
                        cc(nowID-4+test2-2: nowID-4+test2+2)=0;
                        ccc(nowID-4+test2-2: nowID-4+test2+2)=median(ccc([nowID-4+test2-10:nowID-4+test2-5 nowID-4+test2+5:nowID-4+test2+10]));
                        pks=[pks test1];
                        locs=[locs test2+nowID-4-1];
                        clocs=[clocs test2+nowID-4-1];
                    end
                    
                    %                     wind=10;step=1;
                    %                     newc=repmat(0,1,tmp4-1);
                    %                     for q=tmp4:step:tmp2
                    %                         currc=mean(ccc(q:q+wind));
                    %                         newc=[newc currc];
                    %                     end
                    
                    %----------------------------------
                    %                     su=mean(ccc(tmp4:tmp5));
                    %                     testpart=c-su;
                    %                     if plotting2==1
                    %                         figure(200)
                    %                         subplot(4,6,g)
                    %                         plot(f,testpart,'-g')
                    %                         hold on
                    %                     end
                    %                     %                     test=find(cc==0);
                    %                     %                     newpart=cc;
                    %                     %                      newpart(test)=0;
                    %                     newpart=ccc-su;
                    %
                    %                     fi=fit(f(tmp3:end)',newpart(tmp3:end)','exp1');
                    %                     aa=fi.a; bb=fi.b;
                    %                     fi=aa*exp(bb*f(tmp3:end));
                    %                     testpart(tmp3:end)=testpart(tmp3:end)-fi;
                    %   testpart=c-su;
                    
                    %                     if plotting2==1
                    %                         figure(200)
                    %                         subplot(4,6,g)
                    %                         plot(f,testpart,'-g')
                    %                         hold on
                    %                     end
                    %-----------------------------------------
                    
                    %
                    if plotting2==1
                        figure(200)
                        subplot(4,6,g)
                        plot(f,c,'-g')
                        hold on
                    end
                    %                     test=find(cc==0);
                    %                     newpart=cc;
                    %                      newpart(test)=0;
                    newpart=ccc;
                    testpart=ccc;
                    fi=fit(f(tmp3:end)',newpart(tmp3:end)','exp1');
                    aa=fi.a; bb=fi.b;
                    fi=aa*exp(bb*f(tmp3:end));
                    ccc(tmp3:end)=ccc(tmp3:end)-fi;
                    su=mean(ccc(tmp4:tmp5));
                    testpart(tmp3:end)=c(tmp3:end)-fi-su;
                    
                    
                    
                    if plotting2==1
                        figure(200)
                        subplot(4,6,g)
                        %                         plot(f,newpart,'-c')
                        plot(f(tmp3:end),fi,'-k')
                        plot(f,testpart,'-m')
                        axis([0 17 0 max(max(fi),max(testpart(tmp3:end)))+0.01])
                        title([int2str(currspat),'/',int2str(currfr)])
                    end
                    
                    
                    pks=[];locs=[];clocs=[];
                    cc=testpart; realcurr=currfr;
                    
                    for l=1:3
                        tmp=find(f>realcurr*l);
                        nowID=tmp(1);
                        [test1 test2]=max(cc(nowID-4:nowID+3));%-4, +3 positions
                        cc(nowID-4+test2-3: nowID-4+test2+3)=0;
                        pks=[pks test1];
                        locs=[locs test2+nowID-4-1];
                    end
                    tmp=find(pks<0);
                    pks(tmp)=0;
                    frqs=f(locs);
                    DGpeaksCorr(i,g,1:3)=pks(1:3);
                    pks
                    
                    if plotting2==1
                        figure(200)
                        subplot(4,6,g)
                        for pp=1:3
                            if pks(pp)>0
                                plot(frqs(pp),pks(pp),'*k')
                            end
                        end
                    end
                    
                    
                end
            end
        end
        
    end
    if plotting1==1
        figure(1)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_FFT.bmp']))
        close
    end
    if plotting2==1
        figure(100)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath2,[Dtype,'_',int2str(currID),'_allFFT.bmp']))
        close
        figure(200)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath3,[Dtype,'_',int2str(currID),'_newFFT.bmp']))
        close
    end
end
save(fullfile(superpath,[Dtype,'_DG_extendedFFT']),'DGfourier','DGfourierCorr','DGpeaks','DGpeaksCorr','DGphase','DGphaseCorr','DGfourierF')
%% 16c) DG FFT further analysis
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
% load(fullfile(superpath,[Dtype,'_spontFFT']))
cd(superpath)
savepath=fullfile(superpath,'DG_FFT');
savepath2=fullfile(superpath,'DG_allFFT');
savepath3=fullfile(superpath,'DG_newFFT');
plotting1=0; plotting2=0;

forplot=[9 13 14 17 18 21];
togo=find(goodDG==1);
DGfourier=zeros(length(togo),24,8900);%fourier result
DGfourierF=zeros(length(togo),24,8900); %fourier frequency
DGfourierCorr=zeros(length(togo),24,8900);%corrected fourier result
DGpeaks=zeros(length(togo),24,3);%peaks of first 3 harmonics
DGpeaksCorr=zeros(length(togo),24,3);%peaks of first 3 harmonics of corrected fourier
DGphase=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics
DGphaseCorr=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics of corrected fourier
DGmeanrate=zeros(length(togo),24,12001);
for i=1:length(togo)
    i
    
    if plotting1==1
        figure(1)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_FFT.bmp']))
        close
    end
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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    %     spo=spontFFT(currID,:);
    
    
    thiscell=[];
    fouriers=[];
    rates=[];
    clear f2take
    dgcnt=0;
    for g=1:24
        if ~isempty(find(g==forplot))
            dgcnt=dgcnt+1;
        end
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            
            stimFreq=freqs(g);
            
            if max(meanrate)>0
                currfr=stimInfo{currID,1}(g);
                currspat=stimInfo{currID,2}(g);
                
                meanF=mean(meanrate);
                meanFR=round(meanF*10)/10;
                % fourier
                x=meanrate(2000:14000);
                DGmeanrate(i,g,:)=x;
                if max(x)>0
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt)
                        plot(x)
                        axis tight
                        subplot(5,6,dgcnt+6)
                        plot(x(2000:5000))
                        axis tight
                        title([int2str(currspat),'/',int2str(currfr)])
                    end
                    
                    fs=1000;
                    n=2^nextpow2(numel(x));
                    y=fft(x,n)/numel(x);
                    f=fs/2*linspace(0,1,n/2+1);
                    y=y(1:n/2+1);
                    c=2*abs(y); keepc=c;
                    DGfourier(i,g,1:length(c))=c;
                    DGfourierF(i,g,1:length(f))=f;
                    
                    tmp=find(f>0.6); tmp=tmp(1);
                    tmp2=find(f<17); tmp2=tmp2(end);
                    tmp3=find(f>0.2);tmp3=tmp3(1);
                    tmp4=find(f>0.6);tmp4=tmp4(1);
                    
                    ma=max(c(tmp:tmp2));
                    
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt+12)
                        plot(f,c,'-g')
                        axis([0 17 0 ma])
                    end
                    if plotting2==1
                        figure(100)
                        subplot(4,6,g)
                        plot(f,c,'-g')
                        axis([0 17 0 ma])
                        title([int2str(currspat),'/',int2str(currfr)])
                    end
                    
                    
                    y0=fftshift(y);
                    f0=linspace(0,2,n/2+1)*(fs/4);
                    f0=f0-max(f0)/2;
                    phase = unwrap(angle(y0));
                    phases=phase*180/pi;
                    
                    tmp=find(f>currfr-0.3); tmp=tmp(1);
                    tmp2=find(f<33); tmp2=tmp2(end);
                    tmp5=find(f<17); tmp5=tmp5(end);
                    pks=[];locs=[]; clocs=[];
                    cc=keepc; realcurr=currfr;
                    for l=1:4
                        tmp=find(f>realcurr*l);
                        nowID=tmp(1);
                        [test1 test2]=max(cc(nowID-4:nowID+3));% 1 ~ nowID-4
                        cc(nowID-4+test2-3: nowID-4+test2+3)=0;
                        locs=[locs test2+nowID-4-1];
                        clocs=[clocs test2+nowID-4-1];
                        pks=[pks test1];
                        
                    end
                    if clocs(1)-10<1
                        beg=1;
                    elseif currfr>1
                        beg=clocs(1)-20;
                    else
                        beg=clocs(1)-10;
                    end
                    restc=cc(beg:clocs(end)+10);
                    clocs=clocs-beg+1;
                    
                    for pp=1:3
                        if  pks(pp)<max(restc(clocs(pp)-10:clocs(pp+1)-1))*1.5
                            pks(pp)=0;
                        end
                    end
                    if max(pks)<max(restc)+std(restc)
                        pks(1:4)=0;
                    end
                    pks=pks-std(restc);
                    
                    DGpeaks(i,g,1:3)=pks(1:3);
                    frqs=f(locs);
                    t=0:0.001:12;
                    sines=zeros(12001,4);
                    phs=[];
                    cnt=0;
                    for fr=frqs
                        cnt=cnt+1;
                        tmp=find(f0==fr);
                        phs=[phs;phases(tmp)];
                        s=pks(cnt)*(4/pi)*sin(2*pi*fr*t+phs(cnt));
                        sines(:,cnt)=s;
                    end
                    di=frqs-currfr;
                    [mi id]=min(abs(di));
                    ph=phs(id);
                    di=frqs-2*currfr;
                    [mi id]=min(abs(di));
                    ph=[ph phs(id)];
                    di=frqs-3*currfr;
                    [mi id]=min(abs(di));
                    ph=[ph phs(id)];
                    currval=round((ph(1)-ph(2))/360*100)/100;
                    currval2=round((ph(1)-ph(3))/360*100)/100;
                    DGphase(i,g,1:2)=[currval currval2];
                    
                    
                    if ~isempty(find(g==forplot)) && plotting1==1
                        figure(1)
                        subplot(5,6,dgcnt+12)
                        title(['FFT ',num2str(currval),'/',num2str(currval2),' periods'])
                        
                        
                        for s=1:3
                            subplot(5,6,dgcnt+18)
                            plot(sines(:,s),'-r')
                            hold on
                            subplot(5,6,dgcnt+24)
                            plot(sines(1:3000,s),'-r')
                            hold on
                        end
                        subplot(5,6,dgcnt+18)
                        plot(sum(sines'),'-k','linewidth',2)
                        axis tight
                        axis off
                        if dgcnt==1
                            title('reconstruction')
                        elseif dgcnt==2
                            title('with sines of first 3 harmonics')
                        elseif dgcnt==4
                            title('si=peak*(4/pi)*sin(2pi*f*t+phase(f))')
                        end
                        
                        subplot(5,6,dgcnt+24)
                        plot(sum(sines(1:3000,:)'),'-k','linewidth',2)
                        axis tight
                    end
                    
                    
                    
                    %                     minspo=c-spo;
                    
                    
                    pks=[];locs=[];clocs=[];
                    cc=c; realcurr=currfr;
                    ccc=c;
                    numOfSteps=32/realcurr;
                    for l=1:numOfSteps
                        tmp=find(f>realcurr*l);
                        nowID=tmp(1);
                        [test1 test2]=max(cc(nowID-4:nowID+3));%-4, +3 positions
                        cc(nowID-4+test2-2: nowID-4+test2+2)=0;
                        ccc(nowID-4+test2-2: nowID-4+test2+2)=median(ccc([nowID-4+test2-10:nowID-4+test2-5 nowID-4+test2+5:nowID-4+test2+10]));
                        pks=[pks test1];
                        locs=[locs test2+nowID-4-1];
                        clocs=[clocs test2+nowID-4-1];
                    end
                    
                    
                    
                    %
                    if plotting2==1
                        figure(200)
                        subplot(4,6,g)
                        plot(f,c,'-g')
                        hold on
                    end
                    
                    %---version -mean-std
                    %                                         newpart=ccc;
                    %                                         testpart=ccc;
                    %                                         fi=fit(f(tmp3:end)',newpart(tmp3:end)','exp1');
                    %                                         aa=fi.a; bb=fi.b;
                    %                                         fi=aa*exp(bb*f(tmp3:end));
                    %                                         ccc(tmp3:end)=ccc(tmp3:end)-fi;
                    %                                         su=mean(ccc(tmp4:tmp5));
                    %                                         st=std(ccc(tmp4:tmp5));
                    %                                         testpart(tmp3:end)=c(tmp3:end)-fi-su-st;
                    %---------------------------------------------
                    
                    
                    
                    %-----version ^2
                    newpart=ccc;
                    fi=fit(f(tmp3:end)',newpart(tmp3:end)','exp1');
                    aa=fi.a; bb=fi.b;
                    fi=aa*exp(bb*f(tmp3:end));
                    testpart=c;
                    testpart(tmp3:end)=testpart(tmp3:end)-fi;
                    testpart=testpart.^2;
                    %-------------------------------------
                    
                    
                    
                    
                    
                    
                    %----------version /fit---------------------
                    %                     newpart=ccc;
                    %                     fi=fit(f(tmp3:end)',newpart(tmp3:end)','exp1');
                    %                     aa=fi.a; bb=fi.b;
                    %                     fi=aa*exp(bb*f(tmp3:end));
                    %                     testpart=c;
                    %                     testpart(tmp3:end)=testpart(tmp3:end)./fi;
                    %--------------------------------------------
                    
                    
                    
                    %-----version ^2 but with all set to at least 1
                    newpart=ccc;
                    %                                             fi=fit(f(tmp3:end)',newpart(tmp3:end)','exp1');
                    %                                             aa=fi.a; bb=fi.b;
                    %                                             fi=aa*exp(bb*f(tmp3:end));
                    %                                             testpart=c;
                    %                                             testpart(tmp3:end)=testpart(tmp3:end)-fi;
                    % tmp=find(testpart<0); testpart(tmp)=0; zer=find(testpart==0);
                    % testpart=testpart+1; testpart(zer)=0;
                    %                                             testpart=testpart.^2;
                    %-------------------------------------
                    
                    
                    
                    maxval=max(max(fi),max(testpart(tmp3:tmp5)));
                    maxval=max(maxval,max(c(tmp3:end)));
                    if plotting2==1
                        figure(200)
                        subplot(4,6,g)
                        %                         plot(f,newpart,'-c')
                        plot(f(tmp3:end),fi,'-k')
                        plot(f,testpart,'-m')
                        axis([0 17 0 maxval+0.01])
                        title([int2str(currspat),'/',int2str(currfr)])
                    end
                    
                    
                    pks=[];locs=[];clocs=[];
                    cc=testpart; realcurr=currfr;
                    
                    for l=1:3
                        tmp=find(f>realcurr*l);
                        nowID=tmp(1);
                        [test1 test2]=max(cc(nowID-4:nowID+3));%-4, +3 positions
                        cc(nowID-4+test2-3: nowID-4+test2+3)=0;
                        pks=[pks test1];
                        locs=[locs test2+nowID-4-1];
                    end
                    tmp=find(pks<0);
                    pks(tmp)=0;
                    frqs=f(locs);
                    DGpeaksCorr(i,g,1:3)=pks(1:3);
                    pks;
                    
                    if plotting2==1
                        figure(200)
                        subplot(4,6,g)
                        for pp=1:3
                            if pks(pp)>0
                                plot(frqs(pp),pks(pp),'*k')
                            end
                        end
                    end
                    
                    
                end
            end
        end
        
    end
    if plotting2==1
        figure(100)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath2,[Dtype,'_',int2str(currID),'_allFFT.bmp']))
        close
        figure(200)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath3,[Dtype,'_',int2str(currID),'_newFFT.bmp']))
        close
    end
end
save(fullfile(superpath,[Dtype,'_DG_meanrate']),'DGmeanrate')

save(fullfile(superpath,[Dtype,'_DG_extendedFFT']),'DGfourier','DGfourierCorr','DGpeaks','DGpeaksCorr','DGphase','DGphaseCorr','DGfourierF')
%% ...16d) reformat data from 16c
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
load(fullfile(superpath,[Dtype,'_DG_extendedFFT']))
save(fullfile(superpath,[Dtype,'_DG_extendedFFT_onlyPeaks']),'DGpeaks','DGpeaksCorr','DGphase','DGphaseCorr')
%% ...16e) plot some of the F1/F2 ratios
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
savepath=fullfile(superpath,'DG_FFT');
load(fullfile(superpath,[Dtype,'_DG_extendedFFT_onlyPeaks']))
load(fullfile(superpath,'DGinfo'))
spa=cell2mat(stimInfo(1,2));
tem=cell2mat(stimInfo(1,1));
for i=10:20
    %        for i=1:length(DGpeaks)
    currpeaks=DGpeaksCorr(i,:,:);
    f1=reshape(currpeaks(:,:,1),24,1);
    f2=reshape(currpeaks(:,:,2),24,1);
    tmp=find(f2<0); f2(tmp)=0;
    f1norm=f1/max(f1);
    f2norm=f2/max(f1);
    
    figure
    maxval=0;
    for s=1:length(spats)
        tmp=find(spa==spats(s));
        subplot(2,2,1)
        plot(repmat(spats(s),1),f1(tmp),'ob','markerfacecolor','b')
        hold on
        plot(repmat(spats(s),1),f2(tmp),'og','markerfacecolor','g')
        title('absolute F1 and F2')
        plot([spats(s) spats(s)],[0 max(max(f1),max(f2))],'-k')
        
        xlabel('spatial period (um)')
        ylabel('absolute F1 and F2')
        
        
        subplot(2,2,3)
        plot(repmat(spats(s),1),f2norm(tmp)./f1norm(tmp),'om','markerfacecolor','m')
        hold on
        title('all normalized F2/F1')
        xlabel('spatial period (um)')
        ylabel('normalized F2/F1')
        plot([spats(s) spats(s)],[0 max(f2norm./f1norm)],'-k')
        
        [so ids]=sort(f1norm(tmp),'descend');
        [so2 ids2]=sort(f2norm(tmp),'descend');
        if so2(1)>so(1)
            so=so2;
            ids=ids2;
        end
        %         [ma id]=max(f1norm(tmp));
        %         [ma id]=max(f2norm(tmp));
        %         if ma2>ma
        %             id=id2;
        %         end
        
        subplot(2,2,2)
        plot(spats(s),f1norm(tmp(ids(1))),'ob','markerfacecolor','b')
        hold on
        plot(spats(s),f1norm(tmp(ids(2))),'ob')
        plot(spats(s),f2norm(tmp(ids(1))),'og','markerfacecolor','g')
        plot(spats(s),f2norm(tmp(ids(2))),'og')
        title('max. F1 and corresponding F2')
        xlabel('spatial period (um)')
        ylabel('normalized F1 and F2')
        legend('max F1','second max F1','corresp. F2','corresp. F2')
        maxval=max(maxval,max(f2norm./f1norm));
        %         maxval=max(maxval,f2norm(tmp(ids(2)))./f1norm(tmp(ids(2))));
        
        
        subplot(2,2,4)
        plot(spats(s),f2norm(tmp(ids(1)))./f1norm(tmp(ids(1))),'or','markerfacecolor','r')
        hold on
        plot(spats(s),f2norm(tmp(ids(2)))./f1norm(tmp(ids(2))),'or')
        title('F2/F1 for max F1')
        xlabel('spatial period (um)')
        ylabel('normalized F2/F1')
        
    end
    
    for su=1:2
        subplot(2,2,su)
        axis([0 4100 0 inf])
    end
    for su=3:4
        subplot(2,2,su)
        plot([0 4100],[0 0],'-k')
        for s=1:length(spats)
            plot([spats(s) spats(s)],[0 max(1,maxval)],'-k')
        end
        axis([0 4100 0 max(1,maxval)])
    end
    subplot(2,2,2)
    for s=1:length(spats)
        plot([spats(s) spats(s)],[0 max(1,maxval)],'-k')
    end
    
end
%% 16f) look only at F1/F2 ratios with amplitude >=1
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
savepath=fullfile(superpath,'DG_FFT');
load(fullfile(superpath,[Dtype,'_DG_extendedFFT_onlyPeaks']))
load(fullfile(superpath,'DGinfo'))
spa=cell2mat(stimInfo(1,2));
tem=cell2mat(stimInfo(1,1));
DGratios=zeros(size(DGpeaksCorr,1),size(DGpeaksCorr,2));
for i=1:size(DGpeaksCorr,1)
    currpeaks=DGpeaksCorr(i,:,:);
    f1=reshape(currpeaks(:,:,1),24,1);
    f2=reshape(currpeaks(:,:,2),24,1);
    tmp=find(f2<0); f2(tmp)=0;
    f1norm=f1/max(f1);
    f2norm=f2/max(f1);
    goodvals=find(f2>=1);
    tmp=find(f1(goodvals)>=1);
    goodvals=goodvals(tmp);
    badvals=setdiff(1:length(f1norm),goodvals);
    f1norm(badvals)=0;
    f2norm(badvals)=0;
    ratios=f2norm./f1norm;
    ratios(badvals)=0;
    DGratios(i,:)=ratios;
    
    
end
save(fullfile(superpath,[Dtype,'_DGratios']),'DGratios','stimInfo','DGpeaksCorr')
%% 16g) plot ratios from 16f
load(fullfile(superpath,[Dtype,'_DGratios']))
cd(superpath)
spa=cell2mat(stimInfo(1,2));
tem=cell2mat(stimInfo(1,1));
unspa=unique(spa);
untem=unique(tem);

collmax=zeros(size(DGratios,1),6);
collmax2=zeros(size(DGratios,1),6); cnt=0; figure
for i=1:size(DGratios,1)
    cnt=cnt+1;
    currorig1=DGpeaksCorr(i,:,1);
    currorig2=DGpeaksCorr(i,:,2);
    %     currorig1=reshape(DGpeaksCorr(i,:,1:2),48,1);
    % [~,id]=sort(currorig,'descend');
    
    currdata=DGratios(i,:);
    
    if cnt>36
        cnt=1;
        figure
    end
    subplot(6,6,cnt)
    for s=1:length(unspa)
        tmp=find(spa==unspa(s));
        tmp2=find(currdata(tmp)>0);
        if ~isempty(tmp2)
            plotdata=currdata(tmp(tmp2));
            plot(repmat(unspa(s),1,length(tmp2)),plotdata,'o')
            hold on
            [ma,id]=max([currorig1(tmp(tmp2))' currorig2(tmp(tmp2))']);
            if length(ma)==2
                if ma(1)>ma(2)
                    id=id(1);
                else
                    id=id(2);
                end
            else
                id=1;
            end
            
            plot(unspa(s),plotdata(id),'o','markerfacecolor','b')
            collmax(i,s)=plotdata(id);
            collmax2(i,s)=max(plotdata);
        end
    end
    axis([0 4500 0 max([collmax(i,:) 1])])
    plot([0 4500],[0.5 0.5],'--k')
end

figure
counted=[]; f1resp=[];
subplot(1,2,1)
for s=1:length(unspa)
    tmp=find(collmax(:,s)>0);
    counted=[counted;length(tmp)];
    tester=find(spa==unspa(s));
    [row, col]=find(DGpeaksCorr(:,tester,1)>=1);
    f1resp=[f1resp;length(unique(row))];
    plot(repmat(unspa(s),1,length(tmp)),collmax(tmp,s),'o')
    hold on
    plot(unspa(s),median(collmax(tmp,s)),'or','markerfacecolor','r','markersize',8)
    plot([unspa(s) unspa(s)],[median(collmax(tmp,s))-std(collmax(tmp,s))/length(tmp) median(collmax(tmp,s))+std(collmax(tmp,s))/length(tmp)],'-r','linewidth',4)
end
axis([0 4500 0 1.5])
plot([0 4500],[0.5 0.5],'--k')
plot([0 4500],[0.25 0.25],'--k')
plot([0 4500],[1 1],'-k')
% title([int2str(counted'),' out of ',int2str(size(DGratios,1)),' cells'])
title('median +/- SE (F2/F1 for maximal F1 per spatial period considered)')
xlabel('spatial period (um)')
ylabel('F2/F1')
subplot(1,2,2)
for s=1:length(unspa)
    tmp=find(collmax2(:,s)>0);
    plot(repmat(unspa(s),1,length(tmp)),collmax2(tmp,s),'o')
    hold on
    plot(unspa(s),median(collmax2(tmp,s)),'or','markerfacecolor','r','markersize',8)
    plot([unspa(s) unspa(s)],[median(collmax2(tmp,s))-std(collmax2(tmp,s))/length(tmp) median(collmax2(tmp,s))+std(collmax2(tmp,s))/length(tmp)],'-r','linewidth',4)
    
end
axis([0 4500 0 1.5])
plot([0 4500],[0.5 0.5],'--k')
plot([0 4500],[0.25 0.25],'--k')
plot([0 4500],[1 1],'-k')
% title([int2str(f1resp'),' with F1 responses'])
title('median +/- SE (maximal F2/F1 per spatial period considered)')
xlabel('spatial period (um)')
ylabel('F2/F1')

figure
cols='kbcg';
for s=1:length(unspa)
    for t=1:length(untem)
        tmp=find(spa==unspa(s));
        tmp2=find(tem(tmp)==untem(t));
        currdata=DGratios(:,tmp(tmp2));
        tmp=find(currdata>0);
        plot(repmat(unspa(s),1,length(tmp)),currdata(tmp),'^','markeredgecolor',cols(t),'markersize',3)
        hold on
        plot(unspa(s),median(currdata(tmp)),'^','markeredgecolor',cols(t),'markerfacecolor',cols(t),'markersize',10)
    end
    tmp=find(spa==unspa(s));
    currdata=DGratios(:,tmp);
    test=find(currdata>0);
    currdata=currdata(test);
    plot(unspa(s),median(currdata),'or','markerfacecolor','r','markersize',7)
end
axis([0 4500 0 1.5])
plot([0 4500],[0.5 0.5],'--k')
plot([0 4500],[0.25 0.25],'--k')
plot([0 4500],[1 1],'-k')
%% 16h) find cells with an F2/F1>0.5

load(fullfile(superpath,[Dtype,'_DGratios']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
spa=cell2mat(stimInfo(1,2));
tem=cell2mat(stimInfo(1,1));


[row,col]=find(DGratios>0.5);
absid=find(DGratios>0.5);
rawdata=reshape(DGpeaksCorr(:,:,1),size(DGpeaksCorr,1),24);
row2=[]; col2=[];
for r=1:length(row)
    if rawdata(row(r),col(r))>=1
        row2=[row2;row(r)];
        col2=[col2;col(r)];
    end
end
% [row2,col2]=find(rawdata(absid)>=1);
[C,ia,ic] = unique(row2);

goodones=find(goodDG==1);
bigratios=goodones(C);
% spat=spa(col2);
%% 16i) plot various things for cells with F2/F1>0.5
mythreshold=0.4;

load(fullfile(superpath,[Dtype,'_DGratios']))
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
found=0;m=0;
while found==0
    m=m+1;
    curr=cell2mat(stimInfo(m,2));
    if ~isempty(curr)
        found=1;
    end
end
spa=cell2mat(stimInfo(m,2));
uspa=unique(spa);
tem=cell2mat(stimInfo(m,1));
utem=unique(tem);


[row,col]=find(DGratios>=mythreshold);
absid=find(DGratios>=mythreshold);
rawdata=reshape(DGpeaksCorr(:,:,1),size(DGpeaksCorr,1),24);
row2=[]; col2=[];
for r=1:length(row)
    if rawdata(row(r),col(r))>=1
        row2=[row2;row(r)];
        col2=[col2;col(r)];
    end
end
% [row2,col2]=find(rawdata(absid)>=1);
[C,ia,ic] = unique(row2);
goodones=find(goodDG==1);
bigratios=goodones(C); %these are my cells of interest

load(fullfile(superpath,[Dtype,'_DG']))
if  ~exist('DGrates','var')
    load(fullfile(superpath,[Dtype,'_DGrates']))
end
load(fullfile(superpath,[Dtype,'_flash_new']))
flashtogo=togo;
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_chirpFFT']))
chirptogo=togo;
flashtime=[0 2000 4000 6000 8000 11000];
flashcols=[0.5 0.3 0.5 0.9 0.5];

close all

for i=1:33
    %     for i=1:length(bigratios)
    figure
    
    absID=bigratios(i);
    
    relID=C(i);
    currcons=find(row2==relID);
    currcons=col2(currcons); %IDs for DGs with high F2/F1
    [ma maID]=max(DGratios(relID,currcons));
    subplot(4,4,[1 2])
    dgrate=DGrates{absID,currcons(maID)};
    backgr=mean(dgrate(300:1000));
    subdg=dgrate-backgr;
    dgrate=dgrate(2200:14200);
    plot(dgrate,'linewidth',2)
    axis tight
    title('DG response with max F2/F1')
    f0=[];   coll1=[];
    for ii=1:24
        dgrate=DGrates{absID,ii};
        backgr=mean(dgrate(300:1000));
        subdg=dgrate-backgr;
        maxval=max(subdg(2200:14200));
        f0=[f0 maxval];
    end
    f1=DGpeaksCorr(relID,:,1);
    f2=DGpeaksCorr(relID,:,2);
    rat01=f1./f0;
    
    
    flash=info{absID,8};
    if flash~=0
        subplot(4,4,[5 6])
        flashID=find(flashtogo==absID);
        for r=1:5
            rectangle('position',[flashtime(r),0,flashtime(r+1)-flashtime(r),max(flash_rates{flashID})+1],'facecolor',[flashcols(r) flashcols(r) flashcols(r)],'edgecolor',[flashcols(r) flashcols(r) flashcols(r)])
            hold on
        end
        plot(flash_rates{flashID},'-r','linewidth',2)
        axis([0 11000 0 max(flash_rates{flashID})+1])
    end
    
    
    
    coll1=[]; coll2=[]; coll0=[];
    for s=1:length(uspa)
        subplot(4,4,9)
        tmp=find(spa==uspa(s));
        tmp2=find(f1(tmp)>=1);
        tmp=tmp(tmp2);
        if ~isempty(tmp)
            plot(repmat(uspa(s),1,length(tmp)),f1(tmp),'ok')
            hold on
            plot(repmat(uspa(s),1,length(tmp)),f0(tmp),'o','color',[0.7 0.7 0.7])
            plot(uspa(s),max(f1(tmp)),'ok','markerfacecolor','k')
            coll1=[coll1 max(f1(tmp))];
            plot(uspa(s),max(f0(tmp)),'o','markeredgecolor',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7])
            coll0=[coll0 max(f0(tmp))];
        else
            coll1=[coll1 0];
            coll0=[coll0 0];
        end
        
        subplot(4,4,10)
        tmp=find(spa==uspa(s));
        tmp2=find(f2(tmp)>=1);
        tmp=tmp(tmp2);
        if ~isempty(tmp)
            plot(repmat(uspa(s),1,length(tmp)),f2(tmp),'og')
            hold on
            plot(uspa(s),max(f2(tmp)),'og','markerfacecolor','g')
            coll2=[coll2 max(f2(tmp))];
        else
            coll2=[coll2 0];
        end
    end
    subplot(4,4,9)
    plot(uspa,coll1,'-k')
    plot(uspa,coll0,'-','color',[0.7 0.7 0.7])
    title('F1 and F0 (gray) spatial')
    axis([0 4100 0 inf])
    subplot(4,4,10)
    plot(uspa,coll2,'-g')
    title('F2 spatial')
    axis([0 4100 0 inf])
    
    
    
    coll1=[]; coll2=[]; coll0=[];
    for s=1:length(utem)
        subplot(4,4,11)
        tmp=find(tem==utem(s));
        tmp2=find(f1(tmp)>=1);
        tmp=tmp(tmp2);
        if ~isempty(tmp)
            plot(repmat(utem(s),1,length(tmp)),f1(tmp),'ok')
            hold on
            plot(repmat(utem(s),1,length(tmp)),f0(tmp),'o','color',[0.7 0.7 0.7])
            plot(utem(s),max(f1(tmp)),'ok','markerfacecolor','k')
            coll1=[coll1 max(f1(tmp))];
            plot(utem(s),max(f0(tmp)),'o','markeredgecolor',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7])
            coll0=[coll0 max(f0(tmp))];
        else
            coll1=[coll1 0];
            coll0=[coll0 0];
        end
        
        subplot(4,4,12)
        tmp=find(tem==utem(s));
        tmp2=find(f2(tmp)>=1);
        tmp=tmp(tmp2);
        if ~isempty(tmp)
            plot(repmat(utem(s),1,length(tmp)),f2(tmp),'og')
            hold on
            plot(utem(s),max(f2(tmp)),'og','markerfacecolor','g')
            coll2=[coll2 max(f2(tmp))];
        else
            coll2=[coll2 0];
        end
    end
    subplot(4,4,11)
    plot(utem,coll1,'-k')
    plot(utem,coll0,'-','color',[0.7 0.7 0.7])
    title('F1 and F0 (gray) temporal')
    axis([0 8.5 0 inf])
    subplot(4,4,12)
    plot(utem,coll2,'-g')
    title('F2 temporal')
    axis([0 8.5 0 inf])
    
    
    
    chirpID=find(chirptogo==absID);
    if ~isempty(chirpID)
        subplot(4,4,[7 8])
        plot(chirpFFTsmooth{chirpID,2},chirpFFTsmooth{chirpID,1},'-c')
        title('chirp smooth FFT')
        axis tight
    end
    
    
    coll1=[];
    for s=1:length(uspa)
        subplot(4,4,13)
        tmp=find(spa==uspa(s));
        tmp2=find(f1(tmp)>=1);
        tmp=tmp(tmp2);
        if ~isempty(tmp)
            plot(repmat(uspa(s),1,length(tmp)),rat01(tmp),'om')
            hold on
            plot(uspa(s),max(rat01(tmp)),'om','markerfacecolor','m')
            coll1=[coll1 max(rat01(tmp))];
        else
            coll1=[coll1 0];
        end
    end
    plot(uspa,coll1,'-m')
    axis([0 4100 0 inf])
    title('F1/F0')
    
    
    
    
    coll1=[];
    for s=1:length(uspa)
        subplot(4,4,14)
        tmp=find(spa==uspa(s));
        tmp2=find(DGratios(relID,tmp)>=0);
        tmp=tmp(tmp2);
        if ~isempty(tmp)
            plot(repmat(uspa(s),1,length(tmp)),DGratios(relID,tmp),'or')
            hold on
            plot(uspa(s),max(DGratios(relID,tmp)),'or','markerfacecolor','r')
            coll1=[coll1 max(DGratios(relID,tmp))];
        else
            coll1=[coll1 0];
        end
    end
    plot(uspa,coll1,'-r')
    axis([0 4100 0 inf])
    title('F2/F1')
    
    
    
    coll1=[];
    for s=1:length(utem)
        subplot(4,4,15)
        tmp=find(tem==utem(s));
        tmp2=find(f1(tmp)>=1);
        tmp=tmp(tmp2);
        if ~isempty(tmp)
            plot(repmat(utem(s),1,length(tmp)),rat01(tmp),'om')
            hold on
            plot(utem(s),max(rat01(tmp)),'om','markerfacecolor','m')
            coll1=[coll1 max(rat01(tmp))];
        else
            coll1=[coll1 0];
        end
    end
    plot(utem,coll1,'-m')
    axis([0 8.5 0 inf])
    title('F1/F0')
    
    
    
    
    coll1=[];
    for s=1:length(utem)
        subplot(4,4,16)
        tmp=find(tem==utem(s));
        tmp2=find(DGratios(relID,tmp)>=0);
        tmp=tmp(tmp2);
        if ~isempty(tmp)
            plot(repmat(utem(s),1,length(tmp)),DGratios(relID,tmp),'or')
            hold on
            plot(utem(s),max(DGratios(relID,tmp)),'or','markerfacecolor','r')
            coll1=[coll1 max(DGratios(relID,tmp))];
        else
            coll1=[coll1 0];
        end
    end
    plot(utem,coll1,'-r')
    axis([0 8.5 0 inf])
    title('F2/F1')
    
end
%% 18a) DG FFT further analysis: new FFT calc
% DGpeaks contains all peaks
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
load(fullfile(superpath,[Dtype,'_DGrates']))
cd(superpath)
savepath=fullfile(superpath,'DG_newFFTcalc');
savepath3=fullfile(superpath,'DG_FFTreconstr');
plotting1=1; plotting2=1;

forplot=[9 13 14 17 18 21];
togo=find(goodDG==1);
DGfourier=zeros(length(togo),24,8900);%fourier result
DGfourierF=zeros(length(togo),24,8900); %fourier frequency
DGfourierCorr=zeros(length(togo),24,8900);%corrected fourier result
DGpeaks=zeros(length(togo),24,3);%peaks of first 3 harmonics
DGpeaksCorr=zeros(length(togo),24,3);%peaks of first 3 harmonics of corrected fourier
DGphase=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics
DGphaseCorr=zeros(length(togo),24,2);%phase shifts for 2nd and 3rd harmonics of corrected fourier
for i=1:length(togo)
    i
    
    if plotting1==1
        figure(1)
        set(gcf,'position',[1 31 1600 794])
        %         saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_FFT.bmp']))
        %         close
    end
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    
    
    
    dgcnt=0;
    for g=1:24
        if ~isempty(find(g==forplot))
            dgcnt=dgcnt+1;
        end
        
        
        
        meanrate=DGrates{currID,g};
        if ~isempty(meanrate)
            x=meanrate(2200:14200);
            if max(x)>0
                currfr=stimInfo{currID,1}(g);
                currspat=stimInfo{currID,2}(g);
                
                meanF=mean(x);
                meanFR=round(meanF*10)/10;
                % fourier
                
                
                if ~isempty(find(g==forplot)) && plotting1==1
                    figure(1)
                    subplot(5,6,dgcnt)
                    plot(x)
                    axis tight
                    subplot(5,6,dgcnt+6)
                    plot(x(2000:5000))
                    axis tight
                    title([int2str(currspat),'/',int2str(currfr)])
                end
                
                fs=1000;
                n=8*fs;
                y=fft(x,n)/numel(x);
                
                f=fs/2*linspace(0,1,n/2+1);
                y=y(1:n/2+1);
                c=2*abs(y); keepc=c;
                DGfourier(i,g,1:length(c))=c;
                DGfourierF(i,g,1:length(f))=f;
                
                tmp=find(f>0.6); tmp=tmp(1);
                tmp2=find(f<17); tmp2=tmp2(end);
                tmp3=find(f>0.2);tmp3=tmp3(1);
                tmp4=find(f>0.6);tmp4=tmp4(1);
                
                ma=max(c(tmp:tmp2));
                
                if ~isempty(find(g==forplot)) && plotting1==1
                    figure(1)
                    subplot(5,6,dgcnt+12)
                    plot(f,c,'-g')
                    axis([0 17 0 max(ma,0.1)])
                end
                
                
                
                fi=fit(f(tmp3:end)',c(tmp3:end)','exp1');
                aa=fi.a; bb=fi.b;
                fi=aa*exp(bb*f(tmp3:end));
                new=c(tmp3:end)-fi;
                no=find(new<0);
                new(no)=0;
                st=std(new);
                keepc=c;
                c=[repmat(0,1,tmp3-1) new];
                
                DGfourierCorr(i,g,1:length(c))=c;
                
                
                f0 = (-n/4:n/4)*(fs/n);
                y0=fftshift(y);
                %                 f0=linspace(0,2,n/2+1)*(fs/4);
                %                 f0=f0-max(f0)/2;
                %                 f0=round(f0*1000)/1000;
                phases = angle(y0);
                
                
                if plotting2==1
                    %                        if plotting2==1:end
                    figure(200)
                    subplot(4,6,g)
                    plot(f(tmp3:end),keepc(tmp3:end),'-g')
                    hold on
                    plot(f(tmp3:end),fi,'-k')
                    plot(f(tmp3:end),c(tmp3:end),'-m')
                    axis([0 17 0 max(max(keepc(tmp3:end)),0.1)])
                    title([int2str(currspat),'/',int2str(currfr)])
                end
                
                tmp=find(f>currfr-0.3); tmp=tmp(1);
                tmp2=find(f<33); tmp2=tmp2(end);
                pks=[];locs=[];corrpks=[];
                restc=c;
                for l=1:4
                    tmp=find(f==currfr*l);
                    [ma id]=max(keepc(tmp-1:tmp+1));
                    locs=[locs tmp-1+id-1];
                    pks=[pks keepc(tmp-1+id-1)];
                    corrpks=[corrpks c(tmp-1+id-1)];
                    restc(tmp-1+id-1-2:tmp-1+id-1+2)=0;
                end
                
                for pp=1:3
                    if corrpks(pp)<c(locs(pp)-1) || corrpks(pp)<c(locs(pp)+1)
                        corrpks(pp)=0;
                        pks(pp)=0;
                    end
                    if  corrpks(pp)<max(restc(locs(pp)-6:tmp2))*1.5
                        corrpks(pp)=0;
                        pks(pp)=0;
                    end
                end
                
                DGpeaks(i,g,1:3)=pks(1:3);
                DGpeaksCorr(i,g,1:3)=corrpks(1:3);
                
                t=0:0.001:12;
                cosines=zeros(12001,4);
                phs=[];
                cnt=0;
                for fr=f(locs(1:3))
                    cnt=cnt+1;
                    tmp=find(f0==fr);
                    phs=[phs;phases(tmp)];
                    s=pks(cnt)*(4/pi)*cos(2*pi*fr*t+phs(cnt));
                    cosines(:,cnt)=s;
                end
                
                
                DGphase(i,g,1:2)=phs(1:2);
                
                
                if ~isempty(find(g==forplot)) && plotting1==1
                    figure(1)
                    
                    
                    for s=1:3
                        subplot(5,6,dgcnt+18)
                        plot(cosines(:,s),'-r')
                        hold on
                        subplot(5,6,dgcnt+24)
                        plot(cosines(1:3000,s),'-r')
                        hold on
                    end
                    subplot(5,6,dgcnt+18)
                    plot(sum(cosines'),'-k','linewidth',2)
                    axis tight
                    axis off
                    if dgcnt==1
                        title('reconstruction')
                    elseif dgcnt==2
                        title('with 3 first harmonics')
                    elseif dgcnt==4
                        title('si=peak*(4/pi)*cos(2pi*f*t+phase(f))')
                    end
                    
                    subplot(5,6,dgcnt+24)
                    plot(sum(cosines(1:3000,:)'),'-k','linewidth',2)
                    axis tight
                end
                
                
                
                
                
                if plotting2==1
                    figure(200)
                    subplot(4,6,g)
                    for pp=1:3
                        if pks(pp)>0
                            plot(currfr*pp,corrpks(pp),'*k')
                        end
                    end
                end
            end
        end
    end
    
    
    
    
    
    if plotting2==1
        figure(200)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_allFFT.fig']))
        saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_allFFT.bmp']))
        close
        figure(1)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath3,[Dtype,'_',int2str(currID),'_newFFT.bmp']))
        close
    end
end
save(fullfile(superpath,[Dtype,'_DG_newFFTcalc']),'DGfourier','DGfourierCorr','DGpeaks','DGpeaksCorr','DGphase','DGphaseCorr','DGfourierF')
save(fullfile(superpath,[Dtype,'_DG_newFFTsimple']),'DGpeaks','DGpeaksCorr','DGphase')
%% 18b) look only at F1/F2 ratios with amplitude >=thr
thr=0.6;
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
savepath=fullfile(superpath,'DG_FFT');
load(fullfile(superpath,[Dtype,'_DG_newFFTsimple']))
load(fullfile(superpath,'DGinfo'))
spa=cell2mat(stimInfo(1,2));
tem=cell2mat(stimInfo(1,1));
DGratios=zeros(size(DGpeaksCorr,1),size(DGpeaksCorr,2));
for i=1:size(DGpeaksCorr,1)
    currpeaks=DGpeaksCorr(i,:,:);
    f1=reshape(currpeaks(:,:,1),24,1);
    f2=reshape(currpeaks(:,:,2),24,1);
    tmp=find(f2<0); f2(tmp)=0;
    f1norm=f1/max(f1);
    f2norm=f2/max(f1);
    goodvals=find(f2>=thr);
    tmp=find(f1(goodvals)>=thr);
    goodvals=goodvals(tmp);
    badvals=setdiff(1:length(f1norm),goodvals);
    f1norm(badvals)=0;
    f2norm(badvals)=0;
    ratios=f2norm./f1norm;
    ratios(badvals)=0;
    DGratios(i,:)=ratios;
    
    
end
save(fullfile(superpath,[Dtype,'_DGratiosNew']),'DGratios','stimInfo','DGpeaksCorr')
%% 18c) plot ratios from above
load(fullfile(superpath,[Dtype,'_DGratiosNew']))
cd(superpath)
spa=cell2mat(stimInfo(1,2));
tem=cell2mat(stimInfo(1,1));
unspa=unique(spa);
untem=unique(tem);

collmax=zeros(size(DGratios,1),6);
collmax2=zeros(size(DGratios,1),6); cnt=0; figure
for i=1:size(DGratios,1)
    cnt=cnt+1;
    currorig1=DGpeaksCorr(i,:,1);
    currorig2=DGpeaksCorr(i,:,2);
    %     currorig1=reshape(DGpeaksCorr(i,:,1:2),48,1);
    % [~,id]=sort(currorig,'descend');
    
    currdata=DGratios(i,:);
    
    if cnt>36
        cnt=1;
        figure
    end
    subplot(6,6,cnt)
    for s=1:length(unspa)
        tmp=find(spa==unspa(s));
        tmp2=find(currdata(tmp)>0);
        if ~isempty(tmp2)
            plotdata=currdata(tmp(tmp2));
            plot(repmat(unspa(s),1,length(tmp2)),plotdata,'o')
            hold on
            [ma,id]=max([currorig1(tmp(tmp2))' currorig2(tmp(tmp2))']);
            if length(ma)==2
                if ma(1)>ma(2)
                    id=id(1);
                else
                    id=id(2);
                end
            else
                id=1;
            end
            
            plot(unspa(s),plotdata(id),'o','markerfacecolor','b')
            collmax(i,s)=plotdata(id);
            collmax2(i,s)=max(plotdata);
        end
    end
    axis([0 4500 0 max([collmax(i,:) 1])])
    plot([0 4500],[0.5 0.5],'--k')
end

figure
counted=[]; f1resp=[];
subplot(1,2,1)
for s=1:length(unspa)
    tmp=find(collmax(:,s)>0);
    counted=[counted;length(tmp)];
    tester=find(spa==unspa(s));
    [row, col]=find(DGpeaksCorr(:,tester,1)>=1);
    f1resp=[f1resp;length(unique(row))];
    plot(repmat(unspa(s),1,length(tmp)),collmax(tmp,s),'o')
    hold on
    plot(unspa(s),median(collmax(tmp,s)),'or','markerfacecolor','r','markersize',8)
    plot([unspa(s) unspa(s)],[median(collmax(tmp,s))-std(collmax(tmp,s))/length(tmp) median(collmax(tmp,s))+std(collmax(tmp,s))/length(tmp)],'-r','linewidth',4)
end
axis([0 4500 0 1.5])
plot([0 4500],[0.5 0.5],'--k')
plot([0 4500],[0.25 0.25],'--k')
plot([0 4500],[1 1],'-k')
% title([int2str(counted'),' out of ',int2str(size(DGratios,1)),' cells'])
title('median +/- SE (F2/F1 for maximal F1 per spatial period considered)')
xlabel('spatial period (um)')
ylabel('F2/F1')
subplot(1,2,2)
for s=1:length(unspa)
    tmp=find(collmax2(:,s)>0);
    plot(repmat(unspa(s),1,length(tmp)),collmax2(tmp,s),'o')
    hold on
    plot(unspa(s),median(collmax2(tmp,s)),'or','markerfacecolor','r','markersize',8)
    plot([unspa(s) unspa(s)],[median(collmax2(tmp,s))-std(collmax2(tmp,s))/length(tmp) median(collmax2(tmp,s))+std(collmax2(tmp,s))/length(tmp)],'-r','linewidth',4)
    
end
axis([0 4500 0 1.5])
plot([0 4500],[0.5 0.5],'--k')
plot([0 4500],[0.25 0.25],'--k')
plot([0 4500],[1 1],'-k')
% title([int2str(f1resp'),' with F1 responses'])
title('median +/- SE (maximal F2/F1 per spatial period considered)')
xlabel('spatial period (um)')
ylabel('F2/F1')

figure
cols='kbcg';
for s=1:length(unspa)
    for t=1:length(untem)
        tmp=find(spa==unspa(s));
        tmp2=find(tem(tmp)==untem(t));
        currdata=DGratios(:,tmp(tmp2));
        tmp=find(currdata>0);
        plot(repmat(unspa(s),1,length(tmp)),currdata(tmp),'^','markeredgecolor',cols(t),'markersize',3)
        hold on
        plot(unspa(s),median(currdata(tmp)),'^','markeredgecolor',cols(t),'markerfacecolor',cols(t),'markersize',10)
    end
    tmp=find(spa==unspa(s));
    currdata=DGratios(:,tmp);
    test=find(currdata>0);
    currdata=currdata(test);
    plot(unspa(s),median(currdata),'or','markerfacecolor','r','markersize',7)
end
axis([0 4500 0 1.5])
plot([0 4500],[0.5 0.5],'--k')
plot([0 4500],[0.25 0.25],'--k')
plot([0 4500],[1 1],'-k')
%% 18d) look at F1, F2, F1+F2+F3, F2/F1: spatial
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
load(fullfile(superpath,[Dtype,'_DG_newFFTsimple']))
load(fullfile(superpath,'DGinfo'))
spat=[100 200 500 1000 2000 4000];
thr=0.6;

tmp=find(DGpeaksCorr<thr);
ma=max(DGpeaksCorr(:,:,1)');
DGpeaksCorr(:,:,:)=DGpeaksCorr(:,:,:)./repmat(ma',1,24,3);
DGpeaksCorr(tmp)=0;



% figure
cnt=0;
for t=[1 2 4 8]
    cnt=cnt+1;
    temp=find(stimInfo{1,1}==t);
    subplot(2,3,cnt)
    for s=spat
        sp=find(stimInfo{1,2}==s);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,1);
        tmp=find(dat>0);
        figure(1)
        subplot(2,3,cnt)
        plot(repmat(s,length(tmp),1),dat(tmp),'o')
        hold on
        plot(s,mean(dat(tmp)),'or','markerfacecolor','r')
        figure(100)
        subplot(2,3,cnt)
        plot(s,mean(dat(tmp)),'ob','markerfacecolor','b')
        hold on
    end
    %                        axis([0 4100 0 max(median(dat(tmp))+std(dat(tmp)),1)])
    figure(1)
    title(['F1: ',int2str(t),' Hz'])
    figure(100)
    title([int2str(t),' Hz'])
end
subplot(2,3,6)
for s=spat
    sp=find(stimInfo{1,2}==s);
    dat=DGpeaksCorr(:,sp,1);
    ma=max(dat');
    su=sum(dat');
    tmp=find(ma>0);
    figure(1)
    subplot(2,3,6)
    plot(repmat(s,length(tmp),1),ma(tmp),'o')
    hold on
    plot(s,mean(ma(tmp)),'or','markerfacecolor','r')
    figure(100)
    subplot(2,3,6)
    semilogx(s,mean(ma(tmp)),'ob','markerfacecolor','b')
    hold on
    semilogx([s s],[mean(ma(tmp))-std(ma(tmp)) mean(ma(tmp))+std(ma(tmp))],'-b')
    axis([0 4100 0 1.5])
end
%                        axis([0 4100 0 max(median(ma(tmp))+std(ma(tmp)),1)])
figure(1)
title(['F1: max'])


% figure
cnt=0;
for t=[1 2 4 8]
    cnt=cnt+1;
    temp=find(stimInfo{1,1}==t);
    subplot(2,3,cnt)
    for s=spat
        sp=find(stimInfo{1,2}==s);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,2);
        tmp=find(dat>0);
        figure(2)
        subplot(2,3,cnt)
        plot(repmat(s,length(tmp),1),dat(tmp),'ok')
        hold on
        plot(s,mean(dat(tmp)),'or','markerfacecolor','r')
        figure(100)
        subplot(2,3,cnt)
        plot(s,mean(dat(tmp)),'ok','markerfacecolor','k')
        hold on
        %               axis([0 4100 0 max(median(dat(tmp))+std(dat(tmp)),1)])
        
        
    end
    figure(2)
    title(['F2: ',int2str(t),' Hz'])
    figure(100)
    axis([0 4100 0 1.5])
end
subplot(2,3,6)
for s=spat
    sp=find(stimInfo{1,2}==s);
    dat=DGpeaksCorr(:,sp,2);
    ma=max(dat');
    tmp=find(ma>0);
    figure(2)
    subplot(2,3,6)
    plot(repmat(s,length(tmp),1),ma(tmp),'ok')
    hold on
    plot(s,mean(ma(tmp)),'or','markerfacecolor','r')
    figure(100)
    subplot(2,3,6)
    semilogx(s,mean(ma(tmp)),'ok','markerfacecolor','k')
    semilogx([s s],[mean(ma(tmp))-std(ma(tmp)) mean(ma(tmp))+std(ma(tmp))],'-k')
    
    hold on
end
%                        axis([0 4100 0 max(median(ma(tmp))+std(ma(tmp)),1)])
figure(2)
title(['F2: max'])



cnt=0;
for t=[1 2 4 8]
    cnt=cnt+1;
    temp=find(stimInfo{1,1}==t);
    subplot(2,3,cnt)
    for s=spat
        sp=find(stimInfo{1,2}==s);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,3);
        tmp=find(dat>0);
        figure(3)
        subplot(2,3,cnt)
        plot(repmat(s,length(tmp),1),dat(tmp),'om')
        hold on
        plot(s,mean(dat(tmp)),'or','markerfacecolor','r')
        figure(100)
        subplot(2,3,cnt)
        plot(s,mean(dat(tmp)),'om','markerfacecolor','m')
        %               axis([0 4100 0 max(median(dat(tmp))+std(dat(tmp)),1)])
        
        
    end
    figure(3)
    title(['F3: ',int2str(t),' Hz'])
end
subplot(2,3,6)
for s=spat
    sp=find(stimInfo{1,2}==s);
    dat=DGpeaksCorr(:,sp,3);
    ma=max(dat');
    tmp=find(ma>0);
    figure(3)
    subplot(2,3,6)
    plot(repmat(s,length(tmp),1),ma(tmp),'om')
    hold on
    semilogx(s,mean(ma(tmp)),'or','markerfacecolor','r')
    figure(100)
    subplot(2,3,6)
    semilogx(s,mean(ma(tmp)),'om','markerfacecolor','m')
end
%                        axis([0 4100 0 max(median(ma(tmp))+std(ma(tmp)),1)])
figure(3)
title(['F3: max'])




% figure
cnt=0;
for t=[1 2 4 8]
    cnt=cnt+1;
    temp=find(stimInfo{1,1}==t);
    subplot(2,3,cnt)
    for s=spat
        sp=find(stimInfo{1,2}==s);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,1);
        tmp=find(dat>0);
        dat2=DGpeaksCorr(:,id,2);
        tmp2=find(dat2(tmp)>0);
        tmp=tmp(tmp2);
        figure(4)
        subplot(2,3,cnt)
        plot(repmat(s,length(tmp),1),dat2(tmp)./dat(tmp),'og')
        hold on
        plot(s,mean(dat2(tmp)./dat(tmp)),'or','markerfacecolor','r')
        %         figure(100)
        %          subplot(2,3,cnt)
    end
    %          axis([0 4100 0 1.5])
    title(['F2/F1: ',int2str(t),' Hz'])
    
end
subplot(2,3,6)
for s=spat
    sp=find(stimInfo{1,2}==s);
    dat=DGpeaksCorr(:,sp,1);
    dat2=DGpeaksCorr(:,sp,2);
    [ma id]=max(dat2');
    dat=dat(:,id);
    dat2=dat2(:,id);
    tmp=find(dat>0);
    tmp2=find(dat2(tmp)>0);
    [ma id]=max(dat2');
    tmp=tmp(tmp2);
    plot(repmat(s,length(tmp),1),dat2(tmp)./dat(tmp),'og')
    hold on
    plot(s,mean(dat2(tmp)./dat(tmp)),'or','markerfacecolor','r')
end
%     axis([0 4100 0 1.5])
title(['F2/F1: max'])

figure
cnt=0;
colldata=zeros(5,6);
for t=[1 2 4 8]
    cnt=cnt+1;
    temp=find(stimInfo{1,1}==t);
    counter=0;
    for s=spat
        counter=counter+1;
        sp=find(stimInfo{1,2}==s);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,1);
        tmp=find(dat>0);
        colldata(cnt,counter)=mean(dat(tmp));
        
    end
end
counter=0;
for s=spat
    counter=counter+1;
    sp=find(stimInfo{1,2}==s);
    dat=DGpeaksCorr(:,sp,1);
    ma=max(dat');
    tmp=find(ma>0);
    colldata(5,counter)=mean(ma(tmp));
end
cnt=0;
for t=[1 2 4 8]
    cnt=cnt+1;
    temp=find(stimInfo{1,1}==t);
    subplot(2,3,cnt)
    counter=0;
    for s=spat
        counter=counter+1;
        sp=find(stimInfo{1,2}==s);
        id=intersect(temp,sp);
        dat=sum(reshape(DGpeaksCorr(:,id,:),size(DGpeaksCorr,1),size(DGpeaksCorr,3))');
        tmp=find(DGpeaksCorr(:,id,1)>0);
        figure(5)
        subplot(2,3,cnt)
        plot(repmat(s,length(tmp),1),dat(tmp),'oc')
        hold on
        plot(s,colldata(cnt,counter),'ob','markerfacecolor','b')
        plot(s,mean(dat(tmp)),'or','markerfacecolor','r')
        figure(100)
        subplot(2,3,cnt)
        plot(s,mean(dat(tmp)),'oc','markerfacecolor','c')
        %               axis([0 4100 0 max(median(dat(tmp))+std(dat(tmp)),1)])
        
    end
    figure(5)
    title(['F1+F2+F3: ',int2str(t),' Hz'])
    
end
subplot(2,3,6)
counter=0;
for s=spat
    counter=counter+1;
    sp=find(stimInfo{1,2}==s);
    dat=DGpeaksCorr(:,sp,2);
    dat=sum(reshape(DGpeaksCorr(:,sp,:),size(DGpeaksCorr,3),length(sp),size(DGpeaksCorr,1)));
    ma=max(reshape(dat,size(dat,2),size(dat,3)));
    tmp=find(ma>0);
    figure(5)
    subplot(2,3,6)
    plot(repmat(s,length(tmp),1),ma(tmp),'oc')
    hold on
    plot(s,colldata(5,counter),'ob','markerfacecolor','b')
    hold on
    plot(s,mean(ma(tmp)),'or','markerfacecolor','r')
    figure(100)
    subplot(2,3,6)
    semilogx(s,mean(ma(tmp)),'oc','markerfacecolor','c')
    hold on
end
figure(5)
title(['F1+F2+F3: max'])
%                        axis([0 4100 0 max(median(ma(tmp))+std(ma(tmp)),1)])
%% 18e) closer analysis of FFT
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
load(fullfile(superpath,[Dtype,'_DG_newFFTsimple']))
dat=reshape(DGpeaksCorr(:,:,2),size(DGpeaksCorr,1),size(DGpeaksCorr,2));
tmp=find(dat<0.6);
dat(tmp)=0;
datF1=reshape(DGpeaksCorr(:,:,1),size(DGpeaksCorr,1),size(DGpeaksCorr,2));
tmp=find(datF1<0.6);
datF1(tmp)=0;
ma=max(datF1');
datF1=datF1./repmat(ma',1,24);
dat=dat./repmat(ma',1,24);
tmp=find(isnan(datF1)); datF1(tmp)=0;
tmp=find(isnan(dat)); dat(tmp)=0;

meds=[]; stds=[];
for g=1:24
    tmp=find(dat(:,g)>0);
    med=median(dat(tmp,g));
    meds=[meds;med];
    stds=[stds;std(dat(tmp,g))];
end
highF2=cell(24,1);
for g=1:24
    tmp=find(dat(:,g)>meds(g)+stds(g));
    highF2{g}=tmp;
end


medsF1=[]; stdsF1=[];
for g=1:24
    tmp=find(datF1(:,g)>0);
    med=median(datF1(tmp,g));
    medsF1=[medsF1;med];
    stdsF1=[stdsF1;std(datF1(tmp,g))];
end
f1vals=cell(24,1);
clear f1zeros
for g=1:24
    curr=datF1(highF2{g,1},g);
    tmp=find(curr>0);
    f1vals{g}=curr(tmp);
    f1zeros(g)=length(curr(curr==0));
end
figure
for g=1:24
    plot(repmat(g,length(f1vals{g}),1),f1vals{g},'o','markerfacecolor','b')
    hold on
    plot(g,medsF1(g),'or','markerfacecolor','r')
    plot([g g],[medsF1(g)-stdsF1(g) medsF1(g)+stdsF1(g)],'-r')
    highF1=sort(datF1(:,g)); highF1=highF1(end-9:end);
    plot(repmat(g,10,1),highF1,'og','markersize',9)
end
title('F1 (blue) for cells/stimuli with high F2 (red: median F1; green: highest F1s)')




meds=[]; stds=[];
for g=1:24
    tmp=find(datF1(:,g)>0);
    med=median(datF1(tmp,g));
    meds=[meds;med];
    stds=[stds;std(datF1(tmp,g))];
end
highF1=cell(24,1);
for g=1:24
    tmp=find(datF1(:,g)>meds(g)+stds(g));
    highF1{g}=tmp;
end

medsF2=[]; stdsF2=[];
for g=1:24
    tmp=find(dat(:,g)>0);
    med=median(dat(tmp,g));
    medsF2=[medsF2;med];
    stdsF2=[stdsF2;std(dat(tmp,g))];
end
f2vals=cell(24,1);
clear f2zeros
for g=1:24
    curr=dat(highF1{g,1},g);
    tmp=find(curr>0);
    f2vals{g}=curr(tmp);
    f2zeros(g)=length(curr(curr==0));
end
figure
for g=1:24
    plot(repmat(g,length(f2vals{g}),1),f2vals{g},'o','markerfacecolor','b')
    hold on
    plot(g,medsF2(g),'or','markerfacecolor','r')
    plot([g g],[medsF2(g)-stdsF2(g) medsF2(g)+stdsF2(g)],'-r')
    highF2=sort(dat(:,g)); highF2=highF2(end-9:end);
    plot(repmat(g,10,1),highF2,'og','markersize',9)
end
title('F2 (blue) for cells/stimuli with high F1 (red: median F2; green: highest F2s)')


figure
for g=1:24
    subplot(4,6,g)
    scatter(datF1(:,g),dat(:,g))
end
%% 18f)
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
load(fullfile(superpath,[Dtype,'_DG_newFFTsimple']))
dat=reshape(DGpeaksCorr(:,:,2),size(DGpeaksCorr,1),size(DGpeaksCorr,2));
tmp=find(dat<0.6);
dat(tmp)=0;
datF1=reshape(DGpeaksCorr(:,:,1),size(DGpeaksCorr,1),size(DGpeaksCorr,2));
tmp=find(datF1<0.6);
datF1(tmp)=0;
ma=max(datF1');
datF1=datF1./repmat(ma',1,24);
dat=dat./repmat(ma',1,24);
tmp=find(isnan(datF1)); datF1(tmp)=0;
tmp=find(isnan(dat)); dat(tmp)=0;
datF3=reshape(DGpeaksCorr(:,:,3),size(DGpeaksCorr,1),size(DGpeaksCorr,2));
tmp=find(datF3<0.6);
datF31=datF3./repmat(ma',1,24);
tmp=find(isnan(datF3)); datF3(tmp)=0;

goodlist=find(goodDG==1);

meds=[];
for g=1:24
    tmp=find(dat(:,g)>0);
    meds=[meds;median(dat(tmp,g))];
end

figure
for g=1:24
    subplot(4,6,g)
    tmp=find(dat(:,g)>meds(g));
    curr=dat(tmp,g);
    currF1=datF1(tmp,g);
    currF3=datF3(tmp,g);
    for c=1:length(curr)
        plot(c,curr(c),'ob','markerfacecolor','b','markersize',8)
        hold on
        plot(c,currF1(c),'ok','markerfacecolor','k')
        plot(c,currF3(c),'og')
        %           plot(c,curr(c)+currF1(c)+currF3(c),'oc')
    end
    if ~isempty(tmp)
        title([int2str(stimInfo{goodlist(tmp(c)),1}(g)),'Hz / ',int2str(stimInfo{goodlist(tmp(c)),2}(g)),'um'])
    end
end


meds=[];
for g=1:24
    tmp=find(datF1(:,g)>0);
    meds=[meds;median(datF1(tmp,g))];
end

figure
for g=1:24
    subplot(4,6,g)
    tmp=find(datF1(:,g)>meds(g));
    currF1=datF1(tmp,g);
    curr=dat(tmp,g);
    currF3=datF3(tmp,g);
    for c=1:length(curr)
        plot(c,curr(c),'ob','markerfacecolor','b')
        hold on
        plot(c,currF1(c),'ok','markerfacecolor','k','markersize',8)
        plot(c,currF3(c),'og')
        %           plot(c,curr(c)+currF1(c)+currF3(c),'oc')
    end
    if ~isempty(tmp)
        title([int2str(stimInfo{goodlist(tmp(c)),1}(g)),'Hz / ',int2str(stimInfo{goodlist(tmp(c)),2}(g)),'um'])
    end
end
%% 18g) look at F1, F2, F1+F2+F3, F2/F1: temporal
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
load(fullfile(superpath,[Dtype,'_DG_newFFTsimple']))
load(fullfile(superpath,'DGinfo'))
spat=[100 200 500 1000 2000 4000];
tempo=[1 2 4 8];
thr=0.6;

tmp=find(DGpeaksCorr<thr);
ma=max(DGpeaksCorr(:,:,1)');
DGpeaksCorr(:,:,:)=DGpeaksCorr(:,:,:)./repmat(ma',1,24,3);
DGpeaksCorr(tmp)=0;



figure(1)
cnt=0;
for s=spat
    % for t=[1 2 4 8]
    cnt=cnt+1;
    temp=find(stimInfo{1,2}==s);
    
    for t=tempo
        sp=find(stimInfo{1,1}==t);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,1);
        tmp=find(dat>0);
        figure(1)
        subplot(3,3,cnt)
        plot(repmat(t,length(tmp),1),dat(tmp),'o')
        hold on
        
        plot(t,mean(dat(tmp)),'or','markerfacecolor','r')
        figure(100)
        subplot(3,3,cnt)
        plot(t,mean(dat(tmp)),'ob','markerfacecolor','b')
        hold on
        plot([t t],[mean(dat(tmp))-std(dat(tmp)) mean(dat(tmp))+std(dat(tmp))],'-b')
        
    end
    %                        axis([0 4100 0 max(median(dat(tmp))+std(dat(tmp)),1)])
    figure(1)
    title(['F1: ',int2str(s),' um'])
    figure(100)
    title([int2str(s),' um'])
end
subplot(3,3,9)
for s=tempo
    sp=find(stimInfo{1,1}==s);
    dat=DGpeaksCorr(:,sp,1);
    ma=max(dat');
    su=sum(dat');
    tmp=find(ma>0);
    figure(1)
    subplot(3,3,9)
    plot(repmat(s,length(tmp),1),ma(tmp),'o')
    hold on
    plot(s,mean(ma(tmp)),'or','markerfacecolor','r')
    figure(100)
    subplot(3,3,9)
    plot(s,mean(ma(tmp)),'ob','markerfacecolor','b')
    hold on
    plot([s s],[mean(ma(tmp))-std(ma(tmp)) mean(ma(tmp))+std(ma(tmp))],'-b')
    hold on
    figure(1)
    subplot(3,3,8)
    plot(repmat(s,length(tmp),1),su(tmp),'o')
    hold on
    plot(s,mean(su(tmp)),'or','markerfacecolor','r')
end
%
% axis([0 4100 0 max(median(ma(tmp))+std(ma(tmp)),1)])
figure(1)
subplot(3,3,9)
title(['F1: max'])
subplot(3,3,8)
title(['F1: sum'])

% figure
cnt=0;
for s=spat
    cnt=cnt+1;
    temp=find(stimInfo{1,2}==s);
    subplot(3,3,cnt)
    for t=tempo
        sp=find(stimInfo{1,1}==t);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,2);
        tmp=find(dat>0);
        figure(2)
        subplot(3,3,cnt)
        plot(repmat(t,length(tmp),1),dat(tmp),'ok')
        hold on
        plot(t,mean(dat(tmp)),'or','markerfacecolor','r')
        %               axis([0 4100 0 max(median(dat(tmp))+std(dat(tmp)),1)])
        figure(100)
        subplot(3,3,cnt)
        plot(t,mean(dat(tmp)),'ok','markerfacecolor','k')
    end
    figure(2)
    title(['F2: ',int2str(s),' um'])
end

for s=tempo
    sp=find(stimInfo{1,1}==s);
    dat=DGpeaksCorr(:,sp,2);
    ma=max(dat');
    tmp=find(ma>0);
    figure(2)
    subplot(3,3,9)
    plot(repmat(s,length(tmp),1),ma(tmp),'ok')
    hold on
    plot(s,mean(ma(tmp)),'or','markerfacecolor','r')
    figure(100)
    subplot(3,3,9)
    plot(s,mean(ma(tmp)),'ok','markerfacecolor','k')
end
%                        axis([0 4100 0 max(median(ma(tmp))+std(ma(tmp)),1)])
figure(2)
title(['F2: max'])




% figure
cnt=0;
for s=spat
    cnt=cnt+1;
    temp=find(stimInfo{1,2}==s);
    subplot(3,3,cnt)
    for t=tempo
        sp=find(stimInfo{1,1}==t);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,3);
        tmp=find(dat>0);
        figure(3)
        subplot(3,3,cnt)
        plot(repmat(t,length(tmp),1),dat(tmp),'om')
        hold on
        plot(t,mean(dat(tmp)),'or','markerfacecolor','r')
        figure(100)
        subplot(3,3,cnt)
        plot(t,mean(dat(tmp)),'om','markerfacecolor','m')
        hold on
        %               axis([0 4100 0 max(median(dat(tmp))+std(dat(tmp)),1)])
        
        
    end
    figure(3)
    title(['F3: ',int2str(s),' um'])
end

for s=tempo
    sp=find(stimInfo{1,1}==s);
    dat=DGpeaksCorr(:,sp,3);
    ma=max(dat');
    tmp=find(ma>0);
    figure(3)
    subplot(3,3,9)
    plot(repmat(s,length(tmp),1),ma(tmp),'om')
    hold on
    plot(s,mean(ma(tmp)),'or','markerfacecolor','r')
    figure(100)
    subplot(3,3,9)
    plot(s,mean(ma(tmp)),'om','markerfacecolor','m')
end
%                        axis([0 4100 0 max(median(ma(tmp))+std(ma(tmp)),1)])
figure(3)
title(['F3: max'])


% figure
cnt=0;
for s=spat
    cnt=cnt+1;
    temp=find(stimInfo{1,2}==s);
    subplot(3,3,cnt)
    for t=tempo
        sp=find(stimInfo{1,1}==t);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,1);
        tmp=find(dat>0);
        dat2=DGpeaksCorr(:,id,2);
        tmp2=find(dat2(tmp)>0);
        tmp=tmp(tmp2);
        figure(4)
        subplot(3,3,cnt)
        plot(repmat(t,length(tmp),1),dat2(tmp)./dat(tmp),'og')
        hold on
        plot(t,mean(dat2(tmp)./dat(tmp)),'or','markerfacecolor','r')
        %         figure(100)
        %         subplot(3,3,cnt)
    end
    %          axis([0 4100 0 1.5])
    title(['F2/F1: ',int2str(s),' um'])
    
end
subplot(3,3,9)
for s=tempo
    sp=find(stimInfo{1,1}==s);
    dat=DGpeaksCorr(:,sp,1);
    dat2=DGpeaksCorr(:,sp,2);
    [ma id]=max(dat2');
    dat=dat(:,id);
    dat2=dat2(:,id);
    tmp=find(dat>0);
    tmp2=find(dat2(tmp)>0);
    [ma id]=max(dat2');
    tmp=tmp(tmp2);
    plot(repmat(s,length(tmp),1),dat2(tmp)./dat(tmp),'og')
    hold on
    plot(s,mean(dat2(tmp)./dat(tmp)),'or','markerfacecolor','r')
end
%     axis([0 4100 0 1.5])
title(['F2/F1: max'])

% figure
cnt=0;
colldata=zeros(6,5);
for s=spat
    cnt=cnt+1;
    temp=find(stimInfo{1,2}==s);
    counter=0;
    for t=tempo
        counter=counter+1;
        sp=find(stimInfo{1,1}==t);
        id=intersect(temp,sp);
        dat=DGpeaksCorr(:,id,1);
        tmp=find(dat>0);
        colldata(cnt,counter)=mean(dat(tmp));
        
    end
end
counter=0;
for s=tempo
    counter=counter+1;
    sp=find(stimInfo{1,1}==s);
    dat=DGpeaksCorr(:,sp,1);
    ma=max(dat');
    tmp=find(ma>0);
    colldata(6,counter)=mean(ma(tmp));
end
cnt=0;
for s=spat
    cnt=cnt+1;
    temp=find(stimInfo{1,2}==s);
    subplot(3,3,cnt)
    counter=0;
    for t=tempo
        counter=counter+1;
        sp=find(stimInfo{1,1}==t);
        id=intersect(temp,sp);
        dat=sum(reshape(DGpeaksCorr(:,id,:),size(DGpeaksCorr,1),size(DGpeaksCorr,3))');
        tmp=find(DGpeaksCorr(:,id,1)>0);
        figure(5)
        subplot(3,3,cnt)
        plot(repmat(t,length(tmp),1),dat(tmp),'oc')
        hold on
        plot(t,colldata(cnt,counter),'ob','markerfacecolor','b')
        plot(t,mean(dat(tmp)),'or','markerfacecolor','r')
        figure(100)
        subplot(3,3,cnt)
        plot(t,mean(dat(tmp)),'oc','markerfacecolor','c')
        axis([0 9 0 1.5])
        
    end
    figure(5)
    title(['F1+F2+F3: ',int2str(s),' um'])
    
end

counter=0;
for s=tempo
    counter=counter+1;
    sp=find(stimInfo{1,1}==s);
    dat=DGpeaksCorr(:,sp,2);
    dat=sum(reshape(DGpeaksCorr(:,sp,:),size(DGpeaksCorr,3),length(sp),size(DGpeaksCorr,1)));
    ma=max(reshape(dat,size(dat,2),size(dat,3)));
    tmp=find(ma>0);
    figure(5)
    subplot(3,3,9)
    plot(repmat(s,length(tmp),1),ma(tmp),'oc')
    hold on
    plot(s,colldata(5,counter),'ob','markerfacecolor','b')
    plot(s,mean(ma(tmp)),'or','markerfacecolor','r')
    figure(100)
    subplot(3,3,9)
    plot(s,mean(ma(tmp)),'oc','markerfacecolor','c')
    axis([0 9 0 1.5])
end
figure(5)
title(['F1+F2+F3: max'])
%                        axis([0 4100 0 max(median(ma(tmp))+std(ma(tmp)),1)])
%% 18h) plot raw data for cells with high F2 at 100-200um (human)
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
load(fullfile(superpath,[Dtype,'_DGrates']))
load(fullfile(superpath,[Dtype,'_DG_newFFTcalc']))
load(fullfile(superpath,'DGinfo'))
spat=[100 200 500 1000 2000 4000];
tempo=[1 2 4 8];
thr=0.6;
goodlist=find(goodDG==1);

tmp=find(DGpeaksCorr<thr);
ma=max(DGpeaksCorr(:,:,1)');
DGpeaksCorr(:,:,:)=DGpeaksCorr(:,:,:)./repmat(ma',1,24,3);
DGpeaksCorr(tmp)=0;

tmp=find(stimInfo{1,2}<500);
tempo=stimInfo{1,1}(tmp);
data=DGpeaksCorr(:,tmp,2);
for i=1:length(data)
    test=find(data(i,:)>0);
    if ~isempty(test)
        figure
        su=length(test);
        for t=1:length(test)
            subplot(2,su,t)
            plot(reshape(DGfourierF(i,tmp(test(t)),5:90),86,1),reshape(DGfourier(i,tmp(test(t)),5:90),86,1))
            title([int2str(tempo(test(t))),'Hz'])
            subplot(2,su,t+su)
            plot(DGrates{goodlist(i),tmp(test(t))},'-k')
        end
    end
end

%% 18i) plot raw data for cells with high F2 at 2000-4000um (mouse)
load(fullfile(superpath,[Dtype,'_goodDG']))
cd(superpath)
load(fullfile(superpath,[Dtype,'_DGrates']))
load(fullfile(superpath,[Dtype,'_DG_newFFTcalc']))
load(fullfile(superpath,'DGinfo'))
spat=[100 200 500 1000 2000 4000];
tempo=[1 2 4 8];
thr=0.6;
goodlist=find(goodDG==1);

tmp=find(DGpeaksCorr<thr);
ma=max(DGpeaksCorr(:,:,1)');
DGpeaksCorr(:,:,:)=DGpeaksCorr(:,:,:)./repmat(ma',1,24,3);
DGpeaksCorr(tmp)=0;

tmp=find(stimInfo{1,2}<500);
tempo=stimInfo{1,1}(tmp);
data=DGpeaksCorr(:,tmp,2);
for i=1:length(data)
    test=find(data(i,:)>0);
    if ~isempty(test)
        figure
        su=length(test);
        for t=1:length(test)
            subplot(2,su,t)
            plot(reshape(DGfourierF(i,tmp(test(t)),5:90),86,1),reshape(DGfourier(i,tmp(test(t)),5:90),86,1))
            title([int2str(tempo(test(t))),'Hz'])
            subplot(2,su,t+su)
            plot(DGrates{goodlist(i),tmp(test(t))},'-k')
        end
    end
end
%% 18j) !!! get spike times for DG into a variable
% col18) fourier peak
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_goodDG']))
togo=find(goodDG==1);
DGspikes=cell(length(togo),24,4);

for i=1:length(togo)
    i
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    currstart=info{currID,4};
    currstop=info{currID,5};
    
    %     hpath=fullfile(mainpath,currd,'HEKA');
    %     hlist=dir(fullfile(hpath,'*phys'));
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
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
    
    for g=1:length(DGlist)
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        
        
        
        fcnt=0;
        for f=goodnow
            fcnt=fcnt+1;
            spikes=unit{1,2}{f,2};
            DGspikes{i,g,fcnt}=spikes;
        end
        
        
    end
    
    
    
end
save(fullfile(superpath,[Dtype,'_DGspikes']),'DGspikes')
%% 19a) !!! chirp FFT: correct
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_chirp']))
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

chirpFFTabs=cell(length(togo),2);
chirpFFTsmooth=cell(length(togo),2);
chirpRates=cell(length(togo),2);
for i=1:length(togo)
    i
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    currstart=info{currID,4};
    currstop=info{currID,5};
    
    %----version if server works
    %     hpath=fullfile(mainpath,currd,'HEKA');
    %     hlist=dir(fullfile(hpath,'*phys'));
    %     upath=fullfile(mainpath,currd,'units');
    %     ulist=dir(fullfile(upath,'*mat'));
    %
    %     found=0;m=0;
    %     while found==0
    %         m=m+1;
    %         if ~isempty(regexp(ulist(m).name,['unit_',curru]))
    %             found=1;
    %         end
    %     end
    %     load(fullfile(upath,ulist(m).name))
    %
    %     hekalist=unit{1,2}(:,1);
    %
    %     chirp=[];
    %     for f=currstart:currstop
    %         if ~isempty(regexp(hekalist{f},'chirp%'))
    %             ftype=1;
    %             chirp=[chirp;1];
    %         else
    %             chirp=[chirp;0];
    %         end
    %     end
    %     goodfiles=currstart:currstop;
    %     goodfiles=goodfiles(chirp==1);
    %
    %     prot=read_header_field_heka(hpath,hekalist{goodfiles(1)},'Stimulus Protocol');
    %     if prot(end,1)>24000
    %         goodfiles=[];
    %     end
    %
    %
    %     if ~isempty(goodfiles)
    %         allrate=[];
    %         for f=goodfiles
    %             spikes=unit{1,2}{f,2};
    %             rate=MEA_spikerates(spikes,40,25000);
    %             allrate=[allrate;rate];
    %         end
    %         meanrate=mean(allrate);
    %         chirpRates{i,1}=meanrate;
    
    %--------------------------
    
    if i==1
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
        
        prot=read_header_field_heka(hpath,hekalist{goodfiles(1)},'Stimulus Protocol');
    end
    
    
    
    meanrate=chirpFreq{currID,1};
    if ~isempty(meanrate) && max(meanrate)>0
        %            sicnt=0;
        %         for sigma=[5 10 20 40 60]
        int=prot(3:end-1,4);
        
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
        
        
        
        
        NFFT=8000;
        Y=fft(int2,NFFT)/numel(int2);
        fS=1000/2*linspace(0,1,NFFT/2+1);
        cS=2*abs(Y(1:NFFT/2+1));
        
        
        rate=meanrate(round(time(1)):round(time(firstpart)));
        NFFT=8000;
        Y=fft(rate,NFFT)/numel(rate);
        f=1000/2*linspace(0,1,NFFT/2+1);
        c=2*abs(Y(1:NFFT/2+1));
        c=c';
        
        figure
        subplot(3,2,1)
        plot(rate)
        title('firing rate')
        subplot(3,2,3)
        plot(f,c)
        tmp1=find(f>0.4);
        tmp1=tmp1(1);
        tmp2=find(f<32);
        tmp2=tmp2(end);
        tmp3=find(f<=7.5);
        tmp3=tmp3(end);
        axis([0 32 0 max(c(tmp1:tmp2))])
        title('FFT response')
        subplot(3,2,5)
        plot(fS,cS,'-k')
        axis([0 32 0 max(cS(tmp1:tmp2))])
        title('FFT stimulus')
        %                         subplot(3,2,2)
        %                         plot(f,c)
        %                         hold on
        %                         plot(f,cS,'-k')
        %                         axis([0 32 0 max(max(cS(tmp1:tmp2)),max(c(tmp1:tmp2)))])
        
        %                 subplot(2,2,2)
        %                 plot(f,c./cS,'-r')
        %                 axis([0 32 0 max(c(tmp1:tmp2)./cS(tmp1:tmp2))])
        %                 title('absolute')
        chirpFFTabs{i,1}=c./cS;
        chirpFFTabs{i,2}=f;
        
        subplot(3,2,2)
        c=c(tmp1:tmp3); cS=cS(tmp1:tmp3);
        cc=c/max(c); ccS=cS/max(cS); ratio=cc./ccS;
        plot(f(tmp1:tmp3),ratio,'-r')
        axis([0 8 0 max(ratio)])
        title('normalized ratio')
        
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
        
        subplot(3,2,4)
        plot(smooF,smoo,'-g')
        hold on
        title('smooth')
        axis tight
        
        chirpFFTsmooth{i,sicnt}=smoo;
        chirpFFTsmooth{i,1}=smooF;
        
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
        close
        
        
    end
    
end
save(fullfile(superpath,[Dtype,'_chirpFFTnew']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates')
%% 19a) !!! BETTER chirp FFT: correct
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_chirp']))
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

chirpFFTabs=cell(length(togo),2);
chirpFFTsmooth=cell(length(togo),2);
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
    
    prot=read_header_field_heka(hpath,hekalist{goodfiles(1)},'Stimulus Protocol');
    if prot(end,1)>24000
        goodfiles=[];
    end
    
    
    if ~isempty(goodfiles)
        sicnt=0;
        for sigma=[5 10 20 40 60]
            sicnt=sicnt+1;
            
            allrate=[];
            for f=goodfiles
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,sigma,25000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
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
            
            
            
            
            NFFT=8000;
            Y=fft(int2,NFFT)/numel(int2);
            fS=1000/2*linspace(0,1,NFFT/2+1);
            cS=2*abs(Y(1:NFFT/2+1));
            
            
            rate=meanrate(round(time(1)):round(time(firstpart)));
            NFFT=8000;
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
            chirpFFTabs{i,1}=f;
            chirpFFTabs{i,sicnt+1}=c./cS;
            
            
            c=c(tmp1:tmp3); cS=cS(tmp1:tmp3);
            cc=c/max(c); ccS=cS/max(cS); ratio=cc./ccS;
            %         figure
            %         subplot(3,2,1)
            %         plot(rate)
            %         title('firing rate')
            %         subplot(3,2,3)
            %         plot(f,c)
            %         tmp1=find(f>0.4);
            %         tmp1=tmp1(1);
            %         tmp2=find(f<32);
            %         tmp2=tmp2(end);
            %         tmp3=find(f<=7.5);
            %         tmp3=tmp3(end);
            %         axis([0 32 0 max(c(tmp1:tmp2))])
            %         title('FFT response')
            %         subplot(3,2,5)
            %         plot(fS,cS,'-k')
            %         axis([0 32 0 max(cS(tmp1:tmp2))])
            %         title('FFT stimulus')
            %
            %         chirpFFTabs{i,1}=c./cS;
            %         chirpFFTabs{i,2}=f;
            %
            %         subplot(3,2,2)
            %         c=c(tmp1:tmp3); cS=cS(tmp1:tmp3);
            %         cc=c/max(c); ccS=cS/max(cS); ratio=cc./ccS;
            %         plot(f(tmp1:tmp3),ratio,'-r')
            %         axis([0 8 0 max(ratio)])
            %         title('normalized ratio')
            
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
            
            %         subplot(3,2,4)
            %         plot(smooF,smoo,'-g')
            %         hold on
            %         title('smooth')
            %         axis tight
            
            chirpFFTsmooth{i,sicnt+1}=smoo;
            chirpFFTsmooth{i,1}=smooF;
            %
            %         set(gcf,'position',[1 31 1600 794])
            %         saveas(gcf,fullfile(savepath,[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
            %         close
        end
        
    end
    
end
save(fullfile(superpath,[Dtype,'_chirpFFTnew']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates')
%% 19b) plotting chirp
load(fullfile(superpath,[Dtype,'_chirpFFTnew']))
figure
coll=[];
for i=1:length(chirpFFTsmooth)
    plot(chirpFFTsmooth{i,2},chirpFFTsmooth{i,1})
    hold on
    coll=[coll;chirpFFTsmooth{i,1}'];
end
plot(chirpFFTsmooth{i,2},mean(coll),'-r','linewidth',2)
plot(chirpFFTsmooth{i,2},mean(coll)-std(coll),'-r')
plot(chirpFFTsmooth{i,2},mean(coll)+std(coll),'-r')
plot(chirpFFTsmooth{i,2},mean(coll)-std(coll)/size(coll,1),'--r')
plot(chirpFFTsmooth{i,2},mean(coll)+std(coll)/size(coll,1),'--r')
% figure(1)
% plot(chirpFFTsmooth{i,2},mean(coll),'-g','linewidth',2)
% plot(chirpFFTsmooth{i,2},mean(coll)-std(coll),'-g')
% plot(chirpFFTsmooth{i,2},mean(coll)+std(coll),'-g')
% plot(chirpFFTsmooth{i,2},mean(coll)-std(coll)/size(coll,1),'--g')
% plot(chirpFFTsmooth{i,2},mean(coll)+std(coll)/size(coll,1),'--g')
%% 19c) plotting for all species
clear

siglist=[5 10 20 40 60];
close all
superpath='C:\Users\Katja\Desktop\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

types='wph';
cols='bgr';
titles={'mouse','pig','human'};

for sig=1:5;
    figure
    for tt=1:length(types)
        Dtype=types(tt);
        load(fullfile(superpath,[Dtype,'_info']))
        load(fullfile(superpath,[Dtype,'_chirpFFTnew']))
        
        
        %     figure
        coll=[]; collsmooth=[];
        currFr=chirpFFTabs{1,1};
        tmp1=find(currFr>0.8);
        tmp2=find(currFr<=7.5);
        currFr=currFr(tmp1(1):tmp2(end));
        for i=1:length(chirpFFTsmooth)
            currdata=chirpFFTabs{i,sig+1};
            currdata=currdata(tmp1(1):tmp2(end));
            coll=[coll currdata];
            subplot(4,3,tt)
            plot(currFr,currdata,'-','color',[0.6 0.6 0.6])
            hold on
            
            currdata=chirpFFTsmooth{i,sig+1};
            collsmooth=[collsmooth currdata];
            subplot(4,3,tt+6)
            plot(chirpFFTsmooth{i,1},currdata,'-','color',[0.6 0.6 0.6])
            hold on
        end
        subplot(4,3,tt+3)
        plot(currFr,median(coll,2),'-','color',cols(tt),'linewidth',3)
        axis tight
        
        subplot(4,3,tt+9)
        plot(chirpFFTsmooth{i,1},median(collsmooth,2),'-','color',cols(tt),'linewidth',3)
        axis tight
        
        xlabel('Hz')
        title('freq. modulation')
    end
    subplot(4,3,1)
    title(['sigma=',int2str(siglist(sig))])
end
%% 20) !!! FFT of spont activity

load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_spont']))
cd(superpath)
spontFFT=zeros(size(info,1),4001);
spontFFTf=zeros(size(info,1),4001);

for i=1:length(info)
    
    currd=info{i,1};
    curru=info{i,2};
    meanrate=spontRates{i};
    
    if ~isempty(meanrate)
        
        fs=1000;
        n=8*fs;
        y=fft(meanrate,n)/numel(meanrate);
        f=fs/2*linspace(0,1,n/2+1);
        y=y(1:n/2+1);
        c=2*abs(y);
        spontFFT(i,:)=c;
        spontFFTf(i,:)=f;
        
    end
end
save(fullfile(superpath,[Dtype,'_spontFFT']),'spontFFT','spontFFTf')
%% 21) ---just testing FFT: after discussion with Mallot: interpolate firing rates to have 2^x

load(fullfile(superpath,[Dtype,'_DGspikes']))
load(fullfile(superpath,'DGinfo'))
% factor=2^10;
% timeInSec=16;
% stepsize=factor/1000;
% lookup=1:stepsize:timeInSec*factor;

targetSR=2^10;
factor=targetSR/1000;
fs1=1000;
fs2=1024;
fact1=4; fact2=6;
plotit=1;

for ce=1:size(DGspikes,1)
    for dg=1:size(DGspikes,2)
        allrate1=[]; allrate2=[];
        for r=1:size(DGspikes,3)
            curr=DGspikes{ce,dg,r};
            rate1=MEA_spikerates(curr,40,14500,fs1);
            allrate1=[allrate1;rate1];
            %             rate2=MEA_spikerates(curr,40,14500,fs2);
            %             allrate2=[allrate2;rate2];
        end
        meanrate1=mean(allrate1);
        %         meanrate2=mean(allrate2);
        x1=meanrate1(2200:14200);
        %         x2=meanrate2(2200:14488);
        n1=fs1*fact1;
        y1=fft(x1,n1)/numel(x1);
        f1=fs1/2*linspace(0,1,n1/2+1);
        y1=y1(1:n1/2+1);
        c1=2*abs(y1);
        n2=fs1*fact2;
        y2=fft(x1,n2)/numel(x1);
        f2=fs1/2*linspace(0,1,n2/2+1);
        y2=y2(1:n2/2+1);
        c2=2*abs(y2);
        y01=fftshift(y1);
        f01 = (-n1/4:n1/4)*(fs1/n1);
        phases1 = angle(y01);
        y02=fftshift(y2);
        f02 = (-n2/4:n2/4)*(fs1/n2);
        phases2 = angle(y02);
        
        if plotit==1
            figure
            subplot(3,2,1)
            plot(x1)
            title([num2str(stimInfo{1}(dg)),' Hz / ',num2str(stimInfo{2}(dg)),' um'])
            subplot(3,2,2)
            plot(x1)
            subplot(3,2,3)
            plot(f1,c1)
            axis([0 17 -inf inf])
            title(['n = ',int2str(fact1),' * fs'])
            tmp1=find(f1==stimInfo{1}(dg));
            tmp2=find(f1==stimInfo{1}(dg)*2);
            hold on
            plot(f1(tmp1),c1(tmp1),'*r')
            plot(f1(tmp2),c1(tmp2),'*r')
            subplot(3,2,4)
            plot(f2,c2)
            axis([0 17 -inf inf])
            title(['n = ',int2str(fact2),' * fs'])
            tmp1=find(f2==stimInfo{1}(dg));
            tmp2=find(f2==stimInfo{1}(dg)*2);
            hold on
            plot(f2(tmp1),c2(tmp1),'*r')
            plot(f2(tmp2),c2(tmp2),'*r')
            
            subplot(3,2,5)
            plot(f01,phases1)
            axis([0 17 -inf inf])
            tmp1=find(f01==stimInfo{1}(dg));
            tmp2=find(f01==stimInfo{1}(dg)*2);
            ratios=(phases1(tmp2)-2*pi*(-2:2))/phases1(tmp1);
            ratios=round(ratios*100)/100;
            title(['phase1 = ',num2str(phases1(tmp1)),' / phase2 = ',num2str(phases1(tmp2)),' || ratio = ',num2str(ratios)])
            %             title(['phase1 = ',num2str(phases1(tmp1)),' / phase2 = ',num2str(phases1(tmp2)),' || ratio = ',num2str(phases1(tmp2)/phases1(tmp1))])
            subplot(3,2,6)
            plot(f02/factor,phases2)
            axis([0 17 -inf inf])
            tmp1=find(f02==stimInfo{1}(dg));
            tmp2=find(f02==stimInfo{1}(dg)*2);
            ratios=(phases2(tmp2)-2*pi*(-2:2))/phases2(tmp1);
            ratios=round(ratios*100)/100;
            title(['phase1 = ',num2str(phases2(tmp1)),' / phase2 = ',num2str(phases2(tmp2)),' || ratio = ',num2str(ratios)])
            %             title(['phase1 = ',num2str(phases2(tmp1)),' / phase2 = ',num2str(phases2(tmp2)),' || ratio = ',num2str(phases2(tmp2)/phases2(tmp1))])
        end
    end
end

%% ----- here starts a new attempt at chirp and DG -----
%% FFT chirp for various sigmas and various frequency resolutions (NFFT)

plotting=0;

sigmas=[5 10 20 40];
NFFTs=[2000 4000 6000 8000 10000 15000];
varInfo=cell(1,6,6);
for n=1:length(NFFTs)
    varInfo{1,1,n}=['frq for NFFT=',int2str(NFFTs(n))];
    for s=1:length(sigmas)-1
        varInfo{1,s+1,n}=['result for sigma=',int2str(sigmas(s)),' and NFFT=',int2str(NFFTs(n))];
    end
    varInfo{1,s+1,n}=['result for binary spike rate and NFFT=',int2str(NFFTs(n))];
end

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
                chirpFFTabs{i,1,ncnt}=f;
                chirpFFTabs{i,sicnt+1,ncnt}=c./cS;
                
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
save(fullfile(superpath,'newAttempt','chirp',[Dtype,'_chirpFFT']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates','varInfo')
%% plot also chirp rates
close all

sigmas=[5 10 20 40];
NFFTs=[2000 4000 6000 8000 10000 15000];

load(fullfile(superpath,'newAttempt','chirp',[Dtype,'_chirpFFT']))

for i=1:length(togo)
    figure
    for si=1:5
        for ni=1:6
            subplot(5,6,(si-1)*6+ni)
            data=chirpRates{i,2};
            data=data(3000:11000);
            plot(data)
            axis([0 8000 0 Inf])
            if si<5
                title(['sigma=',int2str(sigmas(si)),' / NFFT=',int2str(NFFTs(ni))])
            else
                title(['binary / NFFT=',int2str(NFFTs(ni))])
            end
        end
    end
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(superpath,'newAttempt','chirp',[Dtype,'_',int2str(currID),'_chirp.bmp']))
    close
end
%% FFT DG for various sigmas and various frequency resolutions (NFFT)
% run for human and pig!!
dglist=1:24;%200/1Hz, 1000/2Hz, 2000/1Hz, 4000/1Hz
sigmas=[5 10 20 40];
NFFTs=[2000 4000 6000 8000 10000 15000];
varInfo=cell(1,6,6);
for n=1:length(NFFTs)
    varInfo{1,1,n}=['frq for NFFT=',int2str(NFFTs(n))];
    for s=1:length(sigmas)-1
        varInfo{1,s+1,n}=['result for sigma=',int2str(sigmas(s)),' and NFFT=',int2str(NFFTs(n))];
    end
    varInfo{1,s+1,n}=['result for binary spike rate and NFFT=',int2str(NFFTs(n))];
end

co='kbcgrm';

load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_DGspikes']))
load(fullfile(superpath,'DGinfo'))
cd(superpath)


togo=[];
for i=1:length(info)
    if ~isempty(info{i,19})
        if info{i,19}>0
            togo=[togo;i];
        end
    end
end

DGFFTabs=cell(length(togo),2,2,24);
DGRates=cell(length(togo),2,24);
for i=1:length(togo)
    i
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    currstart=info{currID,4};
    currstop=info{currID,5};
    
    
    for dg=1:size(DGspikes,2)
        sicnt=0;
        for si=1:5
            sicnt=sicnt+1;
            if si<5
                sigma=sigmas(si);
                allrate1=[];
                for r=1:size(DGspikes,3)
                    curr=DGspikes{i,dg,r};
                    rate1=MEA_spikerates(curr,sigma,14500);
                    allrate1=[allrate1;rate1];
                end
                meanrate1=mean(allrate1);
                meanrate=meanrate1(2200:14200);
                
            else
                allrate=[];
                for r=1:size(DGspikes,3)
                    spikes=DGspikes{i,dg,r};
                    spikes=round(spikes);
                    rate=zeros(1,14500);
                    rate(spikes)=1;
                    rate(spikes(diff(spikes)==0))=2;
                    allrate=[allrate;rate(1:14500)];
                end
                meanrate1=mean(allrate);
                meanrate=meanrate1(2200:14200);
            end
            
            DGRates{i,sicnt,dg}=meanrate;
            
            
            
            ncnt=0;
            for NFFT=[2000 4000 6000 8000 10000 15000]
                ncnt=ncnt+1;
                
                Y=fft(meanrate,NFFT)/numel(meanrate);
                f=1000/2*linspace(0,1,NFFT/2+1);
                c=2*abs(Y(1:NFFT/2+1));
                c=c';
                
                tmp1=find(f>0.8);
                tmp1=tmp1(1);
                tmp3=find(f<=32);
                tmp3=tmp3(end);
                
                DGFFTabs{i,1,ncnt,dg}=f;
                DGFFTabs{i,sicnt+1,ncnt,dg}=c;
                
                if ~isempty(find(dg==dglist))
                    figure(1)
                    subplot(5,6,(sicnt-1)*6+ncnt)
                    plot(f(tmp1:tmp3),c(tmp1:tmp3)/max(c(tmp1:tmp3)),'color',co(sicnt))
                    hold on
                    currF=stimInfo{1,1}(dg);
                    test=find(f==currF);
                    plot(f(test),c(test)/max(c(tmp1:tmp3)),'*m','markersize',4)
                end
            end
            
            
        end
        if ~isempty(find(dg==dglist))
            figure(1)
            for si=1:5
                for ni=1:6
                    subplot(5,6,(si-1)*6+ni)
                    axis([0 32 0 Inf])
                    if si<5
                        title(['sigma=',int2str(sigmas(si)),' / NFFT=',int2str(NFFTs(ni))])
                    else
                        title(['binary / NFFT=',int2str(NFFTs(ni))])
                    end
                end
            end
            
            set(gcf,'position',[1 31 1600 794])
            saveas(gcf,fullfile(superpath,'newAttempt','DG',[Dtype,'_',int2str(currID),'_',int2str(stimInfo{1,1}(dg)),'_',int2str(stimInfo{1,2}(dg)),'_responseFFT.bmp']))
            close
            
        end
        
        
    end
    figure
    currdata=reshape(DGFFTabs(i,2:end,:,:),size(DGFFTabs,2)-1,size(DGFFTabs,3),24);
    freq=reshape(DGFFTabs(i,1,:,1),1,6)';
    for si=1:5
        for ni=1:6
            now=cell2mat(reshape(currdata(si,ni,:),1,24));
            currF=freq{ni};
            maxnow=[];
            for m=1:24
                fr=stimInfo{1,1}(m);
                test=find(currF==fr);
                maxnow=[maxnow;now(test,m)];
            end
            normthis=maxnow/max(maxnow);
            spats=unique(stimInfo{1,2});
            temps=unique(stimInfo{1,1});
            subplot(5,6,(si-1)*6+ni)
            zi=reshape(normthis,4,6);
            zi=zi';
            zi=flipud(zi);
            imagesc(1:4,1:6,zi)
            colormap(gray)
            set(gca,'YTick',1:6)
            set(gca,'XTick',1:4)
            set(gca,'XTickLabel',temps)
            xlabel('Hz')
            set(gca,'YTickLabel',flipud(spats))
            ylabel('period (um)')
            if si<5
                title(['sigma=',int2str(sigmas(si)),' / NFFT=',int2str(NFFTs(ni))])
            else
                title(['binary / NFFT=',int2str(NFFTs(ni))])
            end
        end
    end
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(superpath,'newAttempt','DG',[Dtype,'_',int2str(currID),'_heatmapsFFT.bmp']))
    close
    
end
for dg=1:24
    FFTabs=reshape(DGFFTabs(:,:,:,dg),size(DGFFTabs,1),size(DGFFTabs,2),size(DGFFTabs,3));
    save(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_',int2str(stimInfo{1,2}(dg)),'_',int2str(stimInfo{1,1}(dg))]),'FFTabs','togo','varInfo','stimInfo','-v7.3')
end
% save(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT']),'DGFFTabs','togo','varInfo','stimInfo','-v7.3')
save(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGrates']),'togo','DGRates','varInfo','stimInfo','-v7.3')