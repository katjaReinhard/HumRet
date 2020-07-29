%% GENREAL COMMENTS
%The script is more or less chronlogical. Some steps might not be necessary or are redundant.

%Cells starting with a number or !! are needed and were relevant for the data in the paper.
%Cells starting with -- or (-) are not needed.


%Careful, some cells might change the variables "h_info", "p_info", and "m_info", which are essential for
%all analysis. Do not execute random cells without safety copy of these variables.

%For many cells, the cells "adjustable variables" and "data2 have to be
%executed first.


%% adjustable variables
clear
close all


Dtype='h'; % h = human, p = pig, w = wildtype mouse
% superpath='S:\data\Katja\allspecies';
% superpath='/gpfs01/muench/data/Katja/allspecies';
% superpath='C:\Users\Katja\Desktop\all_species_temp';
% superpath='S:\data\Katja\all_species';
superpath='R:\Users\Katja\human_retina';
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
        dnum=3;
    case 'p'
        % Pig data
        dates={'20121127','20121211a','20130613a','20130613b','20130613c','20140723d','20120213','20120319a','20120319b'};
        mainpath='S:\data\Katja\Pig';
        starts=[1 1 1 1 1 156 383 1 361];
        stops=[190 190 164 164 164 285 502 270 630];
        NDs=[4 4 4 4 4 4 3 3 3];
        dnum=2;
        
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
        dnum=1;
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

apath='S:\data\Katja\human';
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
    if ~isempty(find(Dtype=='w'))
        file_list=file_list(3:end);
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
%% -> 3) plotting afterwards due to data loss 
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

inquestion=[];
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
            inquestion=[inquestion;i];
            %             winopen(fullfile(superpath,'LFs',[Dtype,'_',int2str(i),'_',currD,'_',currU,'.bmp']))
            %
            %             choice = questdlg('want to correct LF?', ...
            %                 'LF', ...
            %                 'ON','OFF','none','ON');
            %             switch choice
            %                 case 'ON'
            %                     info{i,9}=1;
            %                 case 'OFF'
            %                     info{i,9}=-1;
            %                 case 'none'
            %                     info{i,9}=0;
            %             end
        end
    end
end

% save(fullfile(superpath,[Dtype,'_info']),'info','goodones')
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
        %         break
    else
        info{curID,12}=0;
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
            winopen(fullfile(superpath,'6vel_new_corr',[Dtype,'_',int2str(i),'_',currd,'_',curru,'.bmp']))
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
    if Dtype=='w'
        currstart=currstart+91;
    end
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
        
        
        
        
%         saveas(gcf,fullfile(superpath,'6vel_new_corr',[Dtype,'_',int2str(currID),'_',currd,'_',curru,'.bmp']))
%         close
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
    if ~isempty(find(Dtype=='w'))
        goodfiles=goodfiles(6:end);
    end
    
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
            case 'no'
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
%% !!! => read out some values: new
clear
close all

Dtype='w';
superpath='f:\all_species_temp';
load(fullfile(superpath,[Dtype,'_info']))

good=find(goodones==1);
fullfield=cell2mat(info(good,8));
keep=fullfield; id1=find(fullfield~=0);
fullfield=length(find(fullfield~=0));
LF=cell2mat(info(good,9));
keep2=LF; id2=find(LF~=0);
LF=length(find(LF~=0));
lfOFF=find(keep2==-1);
lfON=find(keep2==1);
flashOFF=find(keep==-1);
flashON=find(keep==1);
flashONOFF=find(keep==2);
DS=cell2mat(info(good,12));
DSkeep=DS;
DS=length(find(DS==1));
DSpercentage=100/length(DSkeep)*DS;
tmp=find(DSkeep==1);
ONDS=find(keep(tmp)==1);
ONOFFDS=find(keep(tmp)==2);

vel61=cell2mat(info(good,16));
vel61=find(vel61>0);
vel62=cell2mat(info(good,17));
vel62=find(vel62>0);
vel6=length(unique([vel61;vel62]));
DG=cell2mat(info(good,18)); id3=find(DG~=0);
DG=length(find(DG>0));
chirp=cell2mat(info(good,19)); id4=find(chirp~=0);
chirp=length(find(chirp>0));

exp=unique(info(:,1));

all=unique([id1;id2;vel61;vel62;id3;id4]);
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
superpath='E:\all_species_temp';
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
                for a=1:size(allrate,1)
                    plot(allrate(a,300:10500),'-','color',[0.6 0.6 0.6],'linewidth',1)
                    hold on
                end
                plot(median(allrate(:,300:10500),1),'-','color',col,'linewidth',3)
                hold on
                if size(allrate,1)>1
                    med=median(allrate(:,300:10500));
                    dev=allrate(:,300:10500)-repmat(med,size(allrate,1),1);
                    dev2=dev.^2;
                    posSTD=[]; negSTD=[];
                    for xx=1:length(med)
                        tmp=find(dev(:,xx)>0);
                        posSTD=[posSTD ((1/(length(tmp)-1))*sum(dev2(tmp,xx))).^(1/2)];
                        tmp=find(dev(:,xx)<0);
                        negSTD=[negSTD ((1/(length(tmp)-1))*sum(dev2(tmp,xx))).^(1/2)];
                    end
                    %                                 plot(min(allrate(:,300:10500)),'-','color',col,'linewidth',1)
                    %                                 plot(max(allrate(:,300:10500)),'-','color',col,'linewidth',1)
                end
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
                    for a=1:size(allrate,1)
                        plot(allrate(a,300:10500),'-','color',[0.6 0.6 0.6],'linewidth',1)
                        hold on
                    end
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
%% 13a) !!! re-analysis of 6vel: without slowest speed

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
          data3=(data3)/(max(data3)-min(data3));
        
        
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
%% 13b) !!! PLOT histograms for speed tuning
clear
close all
superpath='f:\all_species_temp';
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
%% 13c) !!! PLOT 6vel separated by ON/OFF

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
superpath='F:\all_species_temp';
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
    if Dtype=='w'
        currstart=currstart+91;
    end
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
%% 19a) chirp FFT: correct
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

%% ----- HERE STARTS NEW ATTEMPT AT DG AND CHIRP CALCULATIONS -----
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
save(fullfile(superpath,'newAttempt2','chirp',[Dtype,'_chirpFFT']),'chirpFFTabs','chirpFFTsmooth','togo','chirpRates','varInfo')
%% plot also chirp rates
close all

sigmas=[5 10 20 40];
NFFTs=[2000 4000 6000 8000 10000 15000];

load(fullfile(superpath,'newAttempt','chirp',[Dtype,'_chirpFFT']))

for i=1:length(togo)
    figure
    for si=1:5
        
        subplot(2,3,si)
        data=chirpRates{i,si};
        data=data(3000:11000);
        plot(data)
        axis([0 8000 0 Inf])
        if si<5
            title(['sigma=',int2str(sigmas(si))])
        else
            title(['binary'])
        end
    end
    
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(superpath,'newAttempt','chirp',[Dtype,'_',int2str(togo(i)),'_chirp.bmp']))
    close
end
%% FFT DG for various sigmas and various frequency resolutions (NFFT)

% dglist=1:24;
dglist=[1:6 18 19];
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
load(fullfile(superpath,[Dtype,'_goodDG']))
load(fullfile(superpath,'DGinfo'))
cd(superpath)


togo=find(goodDG==1);
% for i=1:length(info)
%     if ~isempty(info{i,18})
%         if info{i,18}>0
%             togo=[togo;i];
%         end
%     end
% end

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
    saveas(gcf,fullfile(superpath,'newAttempt2','DG',[Dtype,'_',int2str(currID),'_heatmapsFFT.bmp']))
    close
    
end
for dg=1:24
    FFTabs=reshape(DGFFTabs(:,:,:,dg),size(DGFFTabs,1),size(DGFFTabs,2),size(DGFFTabs,3));
    save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGFFT_',int2str(stimInfo{1,2}(dg)),'_',int2str(stimInfo{1,1}(dg))]),'FFTabs','togo','varInfo','stimInfo','-v7.3')
end
% save(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT']),'DGFFTabs','togo','varInfo','stimInfo','-v7.3')
save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGrates']),'togo','DGRates','varInfo','stimInfo','-v7.3')

%% plot also DG rates
close all

sigmas=[5 10 20 40];
NFFTs=[2000 4000 6000 8000 10000 15000];

load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGrates']))

for i=1:length(togo)
    figure
    for dg=1:24
        for si=1:5
            
            subplot(3,2,si)
            data=DGRates{i,si,dg};
            plot(data)
            axis([0 12000 0 Inf])
            if si<5
                title(['sigma=',int2str(sigmas(si))])
            else
                title(['binary'])
            end
        end
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(superpath,'newAttempt','DG',[Dtype,'_',int2str(togo(i)),'_',int2str(stimInfo{1,2}(dg)),'_',int2str(stimInfo{1,1}(dg)),'.bmp']))
        close
    end
    
    
end
%% should save FFT for DG also method-wise...
listSpat=[100 200 500 1000 2000 4000];
listFreq=[1 2 4 8];

cnt=0;
for s=listSpat
    for f=listFreq
        cnt=cnt+1;
        disp(int2str(cnt))
        load(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGFFT_',int2str(s),'_',int2str(f)]));
        if s==100 && f==1
            FFTfreq=reshape(FFTabs(:,1,:),size(FFTabs,1),6);%freq
            FFT1=reshape(FFTabs(:,2,:),size(FFTabs,1),6);%5
            FFT2=reshape(FFTabs(:,3,:),size(FFTabs,1),6);%10
            FFT3=reshape(FFTabs(:,4,:),size(FFTabs,1),6);%20
            FFT4=reshape(FFTabs(:,5,:),size(FFTabs,1),6);%40
            FFT5=reshape(FFTabs(:,6,:),size(FFTabs,1),6);%binary
        else
            FFTfreq=cat(3,FFTfreq,reshape(FFTabs(:,1,:),size(FFTabs,1),6));
            FFT1=cat(3,FFT1,reshape(FFTabs(:,2,:),size(FFTabs,1),6));
            FFT2=cat(3,FFT2,reshape(FFTabs(:,3,:),size(FFTabs,1),6));
            FFT3=cat(3,FFT3,reshape(FFTabs(:,4,:),size(FFTabs,1),6));
            FFT4=cat(3,FFT4,reshape(FFTabs(:,5,:),size(FFTabs,1),6));
            FFT5=cat(3,FFT5,reshape(FFTabs(:,6,:),size(FFTabs,1),6));
        end
        clear FFTabs
    end
end

FFTabs=FFT1;
save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGFFT_sigma5']),'FFTabs','stimInfo','togo','-v7.3')
clear FFT1

FFTabs=FFT2;
save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGFFT_sigma10']),'FFTabs','stimInfo','togo','-v7.3')
clear FFT2

FFTabs=FFT3;
save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGFFT_sigma20']),'FFTabs','stimInfo','togo','-v7.3')
clear FFT3

FFTabs=FFT4;
save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGFFT_sigma40']),'FFTabs','stimInfo','togo','-v7.3')
clear FFT4

FFTabs=FFT5;
save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGFFT_binary']),'FFTabs','stimInfo','togo','-v7.3')
clear FFT5

FFTabs=FFTfreq;
save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_DGFFT_freq']),'FFTabs','stimInfo','togo','-v7.3')
clear FFTfreq
%% judge again mouse DG
clear
close all
path1='F:\all_species_temp\newAttempt2\DG';
load(fullfile(path1,'w_DGFFT_binary'))
keep=FFTabs;
FFT=FFTabs(:,6,:);
load(fullfile(path1,'w_DGrates'))
load(fullfile(path1,'w_DGFFT_freq'))
freq=FFTabs(:,6,:);
load(fullfile('F:\all_species_temp','W_info'))

goodies=[];
for t=1:length(togo)
    curr=reshape(FFT(t,1,:),1,24);
    maxs=[];
    for f=1:24
        ma=max(curr{f});
        maxs=[maxs;ma];
    end
    [ma id]=max(maxs);
    
    
    
    figure
    subplot(2,2,[1 2])
    plot(DGRates{t,3,id})
    title('firing rate, sigma=20')
    subplot(2,2,4)
    tmp=find(freq{t,1,id}<=15);
    plot(freq{t,1,id}(tmp),curr{id}(tmp))
    title('max FFT')
    
    
    normthis=maxs/max(maxs);
    
    
    
    zi=reshape(normthis,4,6);
    
    zi=flipud(zi');
    subplot(2,2,3)
    imagesc(1:4,1:6,zi)
    
    xlabel('Hz')
    
    ylabel('period (um)')
    
    set(gcf,'position',[1601 1 1600 824])
    
    
    
    
    choice = questdlg('does the cell respond?', ...
        'DG', ...
        'yes','no','yes');
    switch choice
        case 'yes'
            info{togo(t),18}=maxs(id);
            goodies=[goodies;t];
        case 'OFF'
            info{togo(t),18}=0;
    end
    
    saveas(gcf,fullfile('F:\all_species_temp\newAttempt2\judgeDG',['w_',int2str(togo(t)),'.bmp']))
close
end
   

togo=togo(goodies);
FFTabs=keep(goodies,:,:);
save(fullfile(path1,'w_DGFFT_binary_corr'),'FFTabs','togo')
FFTabs=freq(goodies,:,:);
save(fullfile(path1,'w_DGFFT_freq_corr'),'FFTabs','togo')
DGRates=DGRates(goodies,:,:);
save(fullfile(path1,'w_DGrates_corr'),'DGRates','togo')
save(fullfile('F:\all_species_temp','w_info'),'info','goodones')
%% TEST: extract peaks (f1,f2,f3,1/2 f1) for "difficult" conditions and check if conditions make sense
% binary, NFFT=15000
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_100_1']));
FFT=FFTabs(:,[1 5 6],6); % cell, (freq, sig40,binary) for NFFT=15000
clear FFTabs
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_100_2']));
FFT=cat(3,FFT,FFTabs(:,[1 5 6],6));% cell, (freq, sig40,binary) for NFFT=15000, DG condition
clear FFTabs
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_100_4']));
FFT=cat(3,FFT,FFTabs(:,[1 5 6],6));% cell, (freq, sig40,binary) for NFFT=15000, DG condition
clear FFTabs
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_100_8']));
FFT=cat(3,FFT,FFTabs(:,[1 5 6],6));% cell, (freq, sig40,binary) for NFFT=15000, DG condition
clear FFTabs
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_200_1']));
FFT=cat(3,FFT,FFTabs(:,[1 5 6],6));% cell, (freq, sig40,binary) for NFFT=15000, DG condition
clear FFTabs
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_2000_8']));
FFT=cat(3,FFT,FFTabs(:,[1 5 6],6));% cell, (freq, sig40,binary) for NFFT=15000, DG condition
clear FFTabs
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGrates']));

idlist=[1 2 3 10 11 12];
dglist=[1 2 3 4 5 20];
frs=[1 2 4 8 1 8];
spats={'100','100','100','100','200','2000'};
temps={'1','2','4','8','1','8'};
for i=1:length(togo)
    figure
    for j=1:6
        subplot(6,3,idlist(j))
        plot(DGRates{i,4,dglist(j)}) % plot firing rates with sig40
        axis tight
        
        title([spats{j},' / ',temps{j}],'fontweight','bold')
        
        subplot(6,3,idlist(j)+3)
        data=cell2mat(FFT(i,3,j));
        if ~isempty(data)
            freq=cell2mat(FFT(i,1,j));
            tmp1=find(freq>=0.35);
            tmp2=find(freq<=30);
            data=data(tmp1(1):tmp2(end));
            freq=freq(tmp1(1):tmp2(end));
            
            tempData=data;
            excluded=[];
            for f=[0.5 1 2 3];
                tmp=find(freq==frs(j)*f);
                excluded=[excluded tmp-2:tmp+2];
            end
            included=1:length(tempData);
            start=find(freq<=frs(j)*0.5);
            start=start(end)-5;
            if start<0
                start=1;
            end
            endpoint=find(freq==frs(j)*3);
            endpoint=endpoint+5;
            included=included(start:endpoint);
            % included=setdiff(included,excluded-start);
            included=setdiff(included,excluded);
            tempData=tempData(included);
            
            bckr=mean(tempData);
            bckrSTD=std(tempData);
            thr=bckr+3*bckrSTD;
            peaks=[]; ids=[]; maxs=[];
            testf=[0.5 1 2 3];
            for f=1:4;
                %             frs(j)*testf(f)
                tmp=find(freq==frs(j)*testf(f));
                if isempty(tmp)
                    tst=find(freq<frs(j)*testf(f));
                    tst=tst(end);
                    flist=tst:tst+1;
                    [~,totake]=max(data(flist));
                else
                    flist=tmp-1:tmp+1;
                    [~,totake]=max(data(tmp-1:tmp+1));
                end
                tmp=flist(totake);
                aa=find(freq<=frs(j)*testf(f)-frs(j)/2);
                if isempty(aa)
                    aa=1;
                end
                aa=aa(end);
                bb=find(freq>=frs(j)*testf(f)+frs(j)/2);
                bb=bb(1);
                
                maxs=[maxs max(data(aa:bb))];
                peaks=[peaks data(tmp)];
                ids=[ids tmp];
            end
            
            tmp=find(peaks<=thr);
            peaks(tmp)=0;
            pks=peaks(find(peaks>0));
            is=ids(find(peaks>0));
            
            mx=maxs(find(peaks>0));
            plot(freq,data,'-k','linewidth',2)
            hold on
            plot(freq,repmat(bckr,length(freq),1),'-g')
            plot(freq,repmat(thr,length(freq),1),'-r')
            for p=1:length(pks)
                plot(freq(is(p)),pks(p),'o','markerfacecolor',[0 0.6 0],'markersize',8)
                if pks(p)<mx(p)
                    plot(freq(is(p)),pks(p),'o','markerfacecolor','m','markersize',8)
                end
            end
            axis tight
            
            if peaks(3)
                title(['FFT with binary // F2=',num2str(100/(peaks(2)-thr)*(peaks(3)-thr)),'%'])
            else
                title('FFT with binary // no F2')
            end
            
            
            
            
            
            subplot(6,3,idlist(j)+6)
            data=cell2mat(FFT(i,2,j));
            freq=cell2mat(FFT(i,1,j));
            tmp1=find(freq>=0.35);
            tmp2=find(freq<=30);
            data=data(tmp1(1):tmp2(end));
            freq=freq(tmp1(1):tmp2(end));
            plot(freq,data)
            axis tight
            title('FFT with sig=40')
        end
    end
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(superpath,'newAttempt','DGtest',[Dtype,'_',int2str(togo(i)),'.bmp']))
    close
    
end
%% DG: extract F1, F2, F3 peaks
% binary, NFFT=15000
% threshold=mean(noise)+3STD
% only if highest peak within fn-f/2 to fn+f/2

load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_binary']));
FFT=FFTabs;
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_DGFFT_freq']));

frs=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
FFTpeaks=zeros(size(FFT,1),24,4);
for i=1:length(togo)
    
    for j=1:24
        
        freq=cell2mat(FFTabs(i,6,j));
        data=cell2mat(FFT(i,6,j));
        if ~isempty(data)
            
            tmp1=find(freq>=0.35);
            tmp2=find(freq<=30);
            data=data(tmp1(1):tmp2(end));
            freq=freq(tmp1(1):tmp2(end));
            
            tempData=data;
            excluded=[];
            for f=[0.5 1 2 3];
                tmp=find(freq==frs(j)*f);
                excluded=[excluded tmp-2:tmp+2];
            end
            included=1:length(tempData);
            start=find(freq<=frs(j)*0.5);
            start=start(end)-5;
            if start<0
                start=1;
            end
            endpoint=find(freq==frs(j)*3);
            endpoint=endpoint+5;
            included=included(start:endpoint);
            % included=setdiff(included,excluded-start);
            included=setdiff(included,excluded);
            tempData=tempData(included);
            
            bckr=mean(tempData);
            bckrSTD=std(tempData);
            thr=bckr+3*bckrSTD;
            peaks=[]; ids=[]; maxs=[];
            testf=[0.5 1 2 3];
            for f=1:4;
                %             frs(j)*testf(f)
                tmp=find(freq==frs(j)*testf(f));
                if isempty(tmp)
                    tst=find(freq<frs(j)*testf(f));
                    tst=tst(end);
                    flist=tst:tst+1;
                    [~,totake]=max(data(flist));
                else
                    flist=tmp-1:tmp+1;
                    [~,totake]=max(data(tmp-1:tmp+1));
                end
                tmp=flist(totake);
                aa=find(freq<=frs(j)*testf(f)-frs(j)/2);
                if isempty(aa)
                    aa=1;
                end
                aa=aa(end);
                bb=find(freq>=frs(j)*testf(f)+frs(j)/2);
                bb=bb(1);
                
                maxs=[maxs max(data(aa:bb))];
                peaks=[peaks data(tmp)];
                ids=[ids tmp];
                [t1,t2]=max(data(aa:bb))
                plot(freq(t2+aa-1),t1,'ok')
                plot([freq(aa) freq(bb)],[0.0015 0.0015],'-k')
            end
            
            tmp=find(peaks<=thr);
            peaks(tmp)=0;
            pks=peaks(find(peaks>0));
            is=ids(find(peaks>0));
            
            for p=1:4
                if maxs(p)>peaks(p)
                    peaks(p)=0;
                end
            end
            FFTpeaks(i,j,:)=peaks;
        end
    end
    
    
end
save(fullfile(superpath,'newAttempt2','DG',[Dtype,'_FFTpeaks']),'FFTpeaks','togo')
%% look at DG peaks
clear
close all
superpath='f:\all_species_temp';
singleplot=0;

Dtype='wph';
spatials=[100 200 500 1000 2000 4000];
freqs=[1 2 4 8];
harmlist=[0.5 1 2 3];
lstyle={'-','-','-.','--'};
lwidth=[1 2 1 1];

cols='bgr';
for d=1:3
    load(fullfile(superpath,'newAttempt2','DG',[Dtype(d),'_FFTpeaks']))
    normPeaks=zeros(size(FFTpeaks));
    
    for i=1:size(FFTpeaks,1)
        curr=reshape(FFTpeaks(i,:,:),24,4);
        ma=max(curr(:,2));
        curr=curr./ma;
        normPeaks(i,:,:)=curr;
        %
        if singleplot==1
            for fp=1:4
                scnt=0;
                for s=1:4:24
                    scnt=scnt+1;
                    figure(d*2-1)
                    subplot(6,6,(fp-1)*6+scnt)
                    plot(curr(s:s+3,fp),'o','markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6],'linestyle','none','markersize',3)
                    hold on
                    title([Dtype(d),' / ',int2str(spatials(scnt)),'um / f',num2str(harmlist(fp))])
                end
                
                for s=1:4
                    figure(d*2)
                    subplot(6,4,(fp-1)*4+s)
                    plot(curr(s:4:24,fp),'o','markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6],'linestyle','none','markersize',3)
                    hold on
                    
                    title([Dtype(d),' / ',int2str(freqs(s)),'Hz / f',num2str(harmlist(fp))])
                end
            end
        end
    end
    save(fullfile(superpath,'newAttempt2','DG',[Dtype(d),'_normPeaks']),'normPeaks','togo')
    for fp=1:4
        scnt=0;
        for s=1:4:24
            scnt=scnt+1;
            if singleplot==1
                figure(d*2-1)
                subplot(6,6,(fp-1)*6+scnt)
                plot(mean(normPeaks(:,s:s+3,fp)),'color',cols(d),'linewidth',2)
                hold on
                patch([1:4 4:-1:1],[mean(normPeaks(:,s:s+3,fp))-std(normPeaks(:,s:s+3,fp)) fliplr(mean(normPeaks(:,s:s+3,fp))+std(normPeaks(:,s:s+3,fp)))],cols(d),'facealpha',0.4)
                set(gca,'xtick',1:4)
                set(gca,'xticklabel',freqs)
                xlabel('frequency [Hz]')
                axis([0.9 4.1 0 1.2])
                subplot(6,6,24+scnt)
                plot(mean(normPeaks(:,s:s+3,fp)),'color',cols(d),'linestyle',lstyle{fp})
                hold on
                set(gca,'xtick',1:4)
                set(gca,'xticklabel',freqs)
                xlabel('frequency [Hz]')
                title('all harmonics')
                axis tight
            end
            figure(100)
            subplot(4,6,scnt)
            plot(mean(normPeaks(:,s:s+3,fp)),'color',cols(d),'linestyle',lstyle{fp},'linewidth',lwidth(fp))
            hold on
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',freqs)
            xlabel('frequency [Hz]')
            title(['all harmonics / ',int2str(spatials(scnt)),'um'])
            axis tight
            figure(200)
            subplot(2,6,scnt)
            plot(mean(FFTpeaks(:,s:s+3,fp)),'color',cols(d),'linestyle',lstyle{fp},'linewidth',lwidth(fp))
            hold on
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',freqs)
            xlabel('frequency [Hz]')
            title(['all harmonics / ',int2str(spatials(scnt)),'um'])
            axis tight
        end
        for s=1:4
            if singleplot==1
                figure(d*2)
                subplot(6,4,(fp-1)*4+s)
                plot(mean(normPeaks(:,s:4:24,fp)),'color',cols(d),'linewidth',2)
                hold on
                patch([1:6 6:-1:1],[mean(normPeaks(:,s:4:24,fp))-std(normPeaks(:,s:4:24,fp)) fliplr(mean(normPeaks(:,s:4:24,fp))+std(normPeaks(:,s:4:24,fp)))],cols(d),'facealpha',0.4)
                set(gca,'xtick',1:6)
                set(gca,'xticklabel',spatials)
                xlabel('spatial period [um]')
                axis([0.9 6.1 0 1.2])
                subplot(6,4,16+s)
                plot(mean(normPeaks(:,s:4:24,fp)),'color',cols(d),'linestyle',lstyle{fp})
                hold on
                set(gca,'xtick',1:6)
                set(gca,'xticklabel',spatials)
                xlabel('spatial period [um]')
                title('all harmonics')
                axis tight
            end
            figure(100)
            subplot(4,6,12+s)
            plot(mean(normPeaks(:,s:4:24,fp)),'color',cols(d),'linestyle',lstyle{fp},'linewidth',lwidth(fp))
            hold on
            set(gca,'xtick',1:6)
            set(gca,'xticklabel',spatials)
            xlabel('spatial period [um]')
            title(['all harmonics / ',int2str(freqs(s)),'Hz'])
            axis tight
            figure(200)
            subplot(2,6,6+s)
            plot(mean(FFTpeaks(:,s:4:24,fp)),'color',cols(d),'linestyle',lstyle{fp},'linewidth',lwidth(fp))
            hold on
            set(gca,'xtick',1:6)
            set(gca,'xticklabel',spatials)
            xlabel('spatial period [um]')
            title(['all harmonics / ',int2str(freqs(s)),'Hz'])
            axis tight
        end
    end
    
    scnt=0;
    for s=1:4:24
        scnt=scnt+1;
        if singleplot==1
            figure(d*2-1)
            subplot(6,6,30+scnt)
        end
        rat=normPeaks(:,s:s+3,3)./normPeaks(:,s:s+3,2);
        tmp=find(rat==Inf);
        rat(tmp)=0;
        
        ratmean=[]; ratst=[]; NoC=[];
        for x=1:4
            now=rat(:,x);
            tmp=find(now>0);
            ratmean=[ratmean mean(now(tmp))];
            ratst=[ratst std(now(tmp))];
            NoC=[NoC length(tmp)];
            if singleplot==1
                if ~isempty(tmp)
                    plot(x,now(tmp),'o','markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6],'linestyle','none','markersize',3)
                end
                hold on
            end
        end
        if singleplot==1
            plot(ratmean,'color',cols(d),'linewidth',2)
            hold on
            patch([1:4 4:-1:1],[ratmean-ratst fliplr(ratmean+ratst)],cols(d),'facealpha',0.4)
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',freqs)
            xlabel('frequency [Hz]')
            axis([0.9 4.1 0 1.2])
            NoC=int2str(NoC);
            title(['%F2 (',NoC,' of ', int2str(size(normPeaks,1)),' cells)'])
        end
        figure(100)
        subplot(4,6,6+scnt)
        plot(ratmean,'color',cols(d),'linewidth',2)
        hold on
        patch([1:4 4:-1:1],[ratmean-ratst fliplr(ratmean+ratst)],cols(d),'facealpha',0.2)
        set(gca,'xtick',1:4)
        set(gca,'xticklabel',freqs)
        xlabel('frequency [Hz]')
        axis([0.9 4.1 0 1.2])
        title('%F2')
    end
    for s=1:4
        if singleplot==1
            figure(d*2)
            subplot(6,4,20+s)
        end
        rat=normPeaks(:,s:4:24,3)./normPeaks(:,s:4:24,2);
        tmp=find(rat==Inf);
        rat(tmp)=0;
        ratmean=[]; ratst=[]; NoC=[];
        for x=1:6
            now=rat(:,x);
            tmp=find(now>0);
            ratmean=[ratmean mean(now(tmp))];
            ratst=[ratst std(now(tmp))];
            NoC=[NoC length(tmp)];
            if singleplot==1
                if ~isempty(tmp)
                    plot(x,now(tmp),'o','markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6],'linestyle','none','markersize',3)
                end
                hold on
            end
        end
        if singleplot==1
            plot(ratmean,'color',cols(d),'linewidth',2)
            hold on
            patch([1:6 6:-1:1],[ratmean-ratst fliplr(ratmean+ratst)],cols(d),'facealpha',0.4)
            set(gca,'xtick',1:6)
            set(gca,'xticklabel',spatials)
            xlabel('spatial period [um]')
            NoC=int2str(NoC);
            title(['%F2 (',NoC,' of ', int2str(size(normPeaks,1)),' cells)'])
            axis([0.9 6.1 0 1.2])
        end
        figure(100)
        subplot(4,6,18+s)
        plot(ratmean,'color',cols(d),'linewidth',2)
        hold on
        patch([1:6 6:-1:1],[ratmean-ratst fliplr(ratmean+ratst)],cols(d),'facealpha',0.4)
        set(gca,'xtick',1:6)
        set(gca,'xticklabel',spatials)
        xlabel('spatial period [um]')
        title('%F2')
        axis([0.9 6.1 0 1.2])
        
    end
    
end
if singleplot==1
    for s=1:3
        figure(s*2-1)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(superpath,'newAttempt','stat_figs',[Dtype(s),'_DGfrequency']))
        saveas(gcf,fullfile(superpath,'newAttempt','stat_figs',[Dtype(s),'_DGfrequency.bmp']))
        figure(s*2)
        set(gcf,'position',[1 31 1600 794])
        saveas(gcf,fullfile(superpath,'newAttempt','stat_figs',[Dtype(s),'_DGspatial']))
        saveas(gcf,fullfile(superpath,'newAttempt','stat_figs',[Dtype(s),'_DGspatial.bmp']))
    end
end
% figure(100)
% set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(superpath,'newAttempt','stat_figs','DGcomparison'))
% saveas(gcf,fullfile(superpath,'newAttempt','stat_figs','DGcomparison.bmp'))
% figure(200)
% set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(superpath,'newAttempt','stat_figs','DGcomparisonAbs'))
% saveas(gcf,fullfile(superpath,'newAttempt','stat_figs','DGcomparisonAbs.bmp'))
%% DG further plotting
%this allows to compare the shape of the curves; for "absolute" response
%strength comparison see pictures from cell above
clear
close all
superpath='E:\all_species_temp';

cols='bgr';
Dtype='wph';
spatials=[100 200 500 1000 2000 4000];
freqs=[1 2 4 8];
harmlist=[0.5 1 2 3];


for d=1:3
    load(fullfile(superpath,'newAttempt','DG',[Dtype(d),'_normPeaks']))
    for fp=1:4
        scnt=0;
        for s=1:4:24
            scnt=scnt+1;
            figure(1)
            subplot(4,6,(fp-1)*6+scnt)
            meanPeak=mean(normPeaks(:,s:s+3,fp));
            normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
            plot(normMean,'color',cols(d),'linewidth',2)
            hold on
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',freqs)
            xlabel('frequency [Hz]')
            axis tight
            title([int2str(spatials(scnt)),'um / f',num2str(harmlist(fp))])
        end
        for s=1:4
            figure(2)
            subplot(4,6,(fp-1)*6+s)
            meanPeak=mean(normPeaks(:,s:4:24,fp));
            normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
            plot(normMean,'color',cols(d),'linewidth',2)
            hold on
            
            set(gca,'xtick',1:6)
            set(gca,'xticklabel',spatials)
            xlabel('spatial period [um]')
            axis tight
            title([int2str(freqs(s)),'Hz / f',num2str(harmlist(fp))])
        end
    end
end
%% DG check spatiotemporla tuning for cells with and without flash responses
%red is with flash responses
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,'newAttempt','DG',[Dtype,'_normPeaks']))
dgtogo=togo;
load(fullfile(superpath,'newAttempt','chirp','chirpRatios_8000'))
good=find(goodones==1);
flash=[];
for i=1:length(good)
    if ~isempty(info{good(i),8}) && isempty(find(info{good(i),8}==0))
        flash=[flash;good(i)];
    end
end

doublegood=intersect(dgtogo,flash);
ids=[];
for d=1:length(doublegood)
    tmp=find(dgtogo==doublegood(d));
    ids=[ids;d];
end

dg=normPeaks(ids,:,:);
oids=setdiff(1:length(dgtogo),ids);
odg=normPeaks(oids,:,:);

check=21:24;
figure
id2=[];
for c=1:size(dg,1)
    if ~isempty(find(dg(c,check,2)>0))
        id2=[id2;c];
    end
end
subplot(2,2,1)
for ii=1:length(id2)
    plot(dg(id2(ii),check,2),'-r')
    hold on
end
subplot(2,2,3)
plot(median(dg(id2,check,2)),'-r')
% for cc=1:length(check)
%     tmp=find(dg(id2,check(cc),2)>0);
%     plot(cc,median(dg(id2(tmp),check(cc),2)),'or')
%     hold on
% end
hold on
oid2=[];
for c=1:size(odg,1)
    if ~isempty(find(odg(c,check,2)>0))
        oid2=[oid2;c];
    end
end
subplot(2,2,2)
for ii=1:length(oid2)
    plot(odg(oid2(ii),check,2))
    hold on
end
subplot(2,2,3)
% for cc=1:length(check)
%       tmp=find(odg(oid2,check(cc),2)>0);
%     plot(cc,median(odg(oid2(tmp),check(cc),2)),'d')
%     hold on
% end
plot(median(odg(oid2,check,2)))

% figure
% plot(mean(dg(:,[1 5 9 13 17 21],2)),'-r')
% hold on
% plot(mean(odg(:,[1 5 9 13 17 21],2)))

for j=21:24
    [a,b]=ranksum(dg(:,j,2),odg(:,j,2))
end
% for j=[1 5 9 13 17 21]
% [a,b]=ranksum(dg(:,j,2),odg(:,j,2))
% end

data=colldata2{dnum};
togo=togos{dnum};
doublegood=intersect(togo,flash);
ids=[];
for d=1:length(doublegood)
    tmp=find(togo==doublegood(d));
    ids=[ids;d];
end

chirp=data(ids,:);
oids=setdiff(1:length(togo),ids);
ochirp=data(oids,:);

figure
subplot(2,2,1)
for ii=1:size(chirp,1)
    plot(chirp(ii,:),'-r')
    hold on
end
subplot(2,2,3)
plot(median(chirp),'-r')
hold on
subplot(2,2,2)
for ii=1:size(ochirp,1)
    plot(ochirp(ii,:))
    hold on
end
subplot(2,2,3)
plot(median(ochirp))
pvals=[];
for c=1:size(chirp,2)
    p=ranksum(chirp(:,c),ochirp(:,c));
    pvals=[pvals p];
end
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
%% useless contrast tuning FFT
clear
close all

nfft=4; %2000,4000,6000,8000,10000,15000
nlist=[2000 4000 6000 8000 10000 15000];
superpath='E:\all_species_temp';

Dtype='wph';
cols='bgr';

colldata2=cell(3,1); colldata1a=cell(3,1); colldata1b=cell(3,1);colldata3=cell(3,1);
frs1=cell(3,1); frs2=cell(3,1); togos=cell(3,1);
figure(1)
for d=1:3
    load(fullfile(superpath,'newAttempt','chirp',[Dtype(d),'_chirpFFT']))
    
    resp=cell2mat(reshape(chirpFFTabs(:,6,nfft),size(chirpFFTabs,1),1)')'; %non-smooth binary
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
saveas(gcf,fullfile(superpath,'newAttempt',['chirp_comparison',int2str(nlist(nfft)),'.bmp']))
saveas(gcf,fullfile(superpath,'newAttempt',['chirp_comparison',int2str(nlist(nfft)),'.fig']))
save(fullfile(superpath,'newAttempt','chirp',['chirpRatios_',int2str(nlist(nfft))]),'colldata2','colldata3','colldata1a','colldata1b','frs1','frs2','togos')

figure(2)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(superpath,'newAttempt',['chirp_comparison',int2str(nlist(nfft)),'.eps']))

%% contrast tuning (chirp, 2nd half)

plotting=0;

load(fullfile(superpath,[Dtype,'_info']))
% load(fullfile(superpath,[Dtype,'_chirp']))
cd(superpath)


togo=[];
for i=1:length(info)
    if ~isempty(info{i,19})
        if info{i,19}>0
            togo=[togo;i];
        end
    end
end

contrastMod=cell(length(togo),2);
contrastMod2=cell(length(togo),1);
backMod=cell(length(togo),2);
backMod2=cell(length(togo),1);
NewcontrastMod=cell(length(togo),2);
NewcontrastMod2=cell(length(togo),1);
NewbackMod=cell(length(togo),2);
NewbackMod2=cell(length(togo),1);
contrastRates=cell(length(togo),3);
backRates=cell(length(togo),3);
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
        allrate=[]; allrate40=[];
        for f=goodfiles
            spikes=unit{1,2}{f,2};
            rate=MEA_spikerates(spikes,10,25000);
            allrate=[allrate;rate];
            rate=MEA_spikerates(spikes,40,25000);
            allrate40=[allrate40;rate];
        end
        meanrate=mean(allrate);
        mean40=mean(allrate40);
        backr=meanrate(300:3000);
        back40=mean40(300:3000);
        di=diff(prot(2:end-1,1));
        tmp=find(di>100);
        if length(tmp)==2
            
            start=tmp(1)+2;
            start=round(prot(start,1));
            endpoint=tmp(2)+1;
            endpoint=round(prot(endpoint,1));
            meanrate=meanrate(start:endpoint);
            mean40=mean40(start:endpoint);
            
            int=prot(:,4);
            time=prot(2:end-2,1);
            timediff=diff(time);
            su=sum(timediff);
            timediff=round(timediff);
            roundsum=sum(timediff);
            changing=ceil(roundsum-su);
            stepping=floor(length(timediff)/changing);
            for ti=1:stepping:length(timediff)
                timediff(ti)=16;
            end
            tst=find(timediff>17);
            timediff(tst)=17;
            
            cnt=0;
            int2=[];
            for tt=tmp(1)+2:tmp(2)+1
                cnt=cnt+1;
                int2=[int2;repmat(int(tt),timediff(cnt),1)];
            end
            while length(int2)<8100
                tt=tt+1; cnt=cnt+1;
                int2=[int2;repmat(int(tt),timediff(cnt),1)];
            end
            
            contrastRates{i,1}=meanrate;
            contrastRates{i,2}=int2;
            contrastRates{i,3}=mean40;
            
            backRates{i,1}=backr;
            backRates{i,3}=back40;
            
            window=500;
            keepmean=meanrate;
            meanrate=meanrate-mean(backr);
            absrate=abs(meanrate);
            cnt=1;
            clear param datapoint modu
            param2=zeros(length(1:50:length(meanrate)-window),length(absrate));
            modu2=zeros(length(1:50:length(meanrate)-window),length(absrate));
            for fc=1:50:length(meanrate)-window
                param(cnt)=sum(absrate(fc:fc+window-1))/window; %modulation strength
                datapoint(cnt)=mean(fc:fc+window-1);
                param2(cnt,fc:fc+window-1)=repmat(sum(absrate(fc:fc+window))/window,500,1);
                modu(cnt)=max(keepmean(fc:fc+window-1))-min(keepmean(fc:fc+window-1));
                modu2(cnt,fc:fc+window-1)=repmat(max(keepmean(fc:fc+window-1))-min(keepmean(fc:fc+window-1)),500,1);
                cnt=cnt+1;
                %                 param2=[param2;repmat(sum(absrate(fc:fc+window))/window,500,1)];
            end
            param3=[];
            for p=1:size(param2,2)
                tmp=find(param2(:,p)>0);
                if isempty(tmp)
                    param3=[param3 0];
                else
                    param3=[param3 mean(param2(tmp,p))];
                end
            end
            param3=param3(1:8100);
            param4=[];
            for p=1:size(param2,2)
                tmp=find(param2(:,p)>0);
                if isempty(tmp)
                    param4=[param4 0];
                else
                    param4=[param4 max(param2(tmp,p))];
                end
            end
            param4=param4(1:8100);
            modu3=[];
            for p=1:size(modu2,2)
                tmp=find(modu2(:,p)>0);
                if isempty(tmp)
                    modu3=[modu3 0];
                else
                    modu3=[modu3 mean(modu2(tmp,p))];
                end
            end
            modu3=modu3(1:8100);
            int3=int2/max(int2);
            int3=int3*max(param);
            int4=int2/max(int2);
            int4=int4*max(modu);
            
            keepback=backr;
            backr=backr-mean(backr);
            backabs=abs(backr);
            cnt=1;
            clear backP backDP backM
            backP2=zeros(length(1:50:length(backr)-window),length(backabs));
            backM2=zeros(length(1:50:length(backr)-window),length(backabs));
            for fc=1:50:length(backr)-window
                backP(cnt)=sum(backabs(fc:fc+window-1))/window; %modulation strength
                backDP(cnt)=mean(fc:fc+window-1);
                backP2(cnt,fc:fc+window-1)=repmat(sum(backabs(fc:fc+window))/window,500,1);
                backM(cnt)=max(keepback(fc:fc+window-1))-min(keepback(fc:fc+window-1));
                backM2(cnt,fc:fc+window-1)=repmat(max(keepback(fc:fc+window-1))-min(keepback(fc:fc+window-1)),500,1);
                cnt=cnt+1;
                %                 param2=[param2;repmat(sum(absrate(fc:fc+window))/window,500,1)];
            end
            backP3=[];
            for p=1:size(backP2,2)
                tmp=find(backP2(:,p)>0);
                if isempty(tmp)
                    backP3=[backP3 0];
                else
                    backP3=[backP3 mean(backP2(tmp,p))];
                end
            end
            backP3=backP3(1:2700);
            backP4=[];
            for p=1:size(backP2,2)
                tmp=find(backP2(:,p)>0);
                if isempty(tmp)
                    backP4=[backP4 0];
                else
                    backP4=[backP4 max(backP2(tmp,p))];
                end
            end
            backP4=backP4(1:2700);
            backM3=[];
            for p=1:size(backM2,2)
                tmp=find(backM2(:,p)>0);
                if isempty(tmp)
                    backM3=[backM3 0];
                else
                    backM3=[backM3 mean(backM2(tmp,p))];
                end
            end
            backM3=backM3(1:2700);
            
            if plotting==1
                figure
                subplot(2,3,1)
                plot(back40,'-c')
                hold on
                plot(4001:4000+length(mean40),mean40)
                title('sigma=40')
                subplot(2,3,2)
                plot(backP3,'-c','linewidth',2)
                hold on
                plot(backP4,'-g','linewidth',2)
                plot(4001:4000+length(param3),param3,'linewidth',2)
                plot(4001:4000+length(param4),param4,'k','linewidth',2)
                hold on
                plot(4001:4000+length(int3),int3,'-','color',[0.6 0.6 0.6])
                title('resp. mod. averaged for overlapping parts')
                subplot(2,3,3)
                plot(backM3,'-c','linewidth',2)
                hold on
                plot(4001:4000+length(modu3),modu3,'linewidth',2)
                hold on
                plot(4001:4000+length(int4),int4,'-','color',[0.6 0.6 0.6])
                title('resp. mod. averaged for overlapping parts (method2)')
                subplot(2,3,5)
                plot(backDP,backP,'-c','linewidth',2)
                hold on
                plot(datapoint+4000,param,'linewidth',2)
                title('resp. mod. non-smooth')
                subplot(2,3,6)
                plot(backDP,backM,'-c','linewidth',2)
                hold on
                plot(datapoint+4000,modu,'linewidth',2)
                title('resp. mod. non-smooth (method2)')
                subplot(2,3,4)
                plot(keepback,'-c')
                hold on
                plot(4001:4000+length(keepmean),keepmean)
                title('sigma=10')
                set(gcf,'position',[1 31 1600 794])
                saveas(gcf,fullfile(superpath,'newAttempt','contrast',[Dtype,'_',int2str(togo(i)),'_contrast.bmp']))
                close
            end
            contrastMod{i,1}=param;
            contrastMod{i,2}=datapoint;
            contrastMod2{i,1}=param3;
            contrastMod2{i,2}=param4;
            backMod{i,1}=backP;
            backMod{i,1}=backDP;
            backMod2{i,1}=backP3;
            backMod2{i,1}=backP4;
            NewcontrastMod{i,1}=modu;
            NewcontrastMod{i,2}=datapoint;
            NewcontrastMod2{i,1}=modu3;
            NewbackMod{i,1}=backM;
            NewbackMod{i,1}=backDP;
            NewbackMod2{i,1}=backM3;
        end
    end
end
save(fullfile(superpath,'newAttempt','contrast',[Dtype,'_contrast']),'backRates','NewcontrastMod','NewcontrastMod2','NewbackMod','NewbackMod2','contrastMod','contrastMod2','contrastRates','backMod','backMod2','togo')
%% -> check if real response
% a)2Hz peak
% b)no 2Hz peak in background
% c)only parts with higher modulation than background are real responses

load(fullfile(superpath,'newAttempt','contrast',[Dtype,'_contrast']))
%     contrastRates{i,1}=meanrate;
%     contrastRates{i,2}=int2;
%     contrastRates{i,3}=mean40;
%     contrastMod{i,1}=param; %original
%     contrastMod{i,2}=datapoint; %time axis
%     contrastMod2{i,1}=param3; %av for each overlapping window
%     contrastMod2{i,2}=param4; %max for each overlapping window
%     backMod{i,1}=backP;
%     backMod{i,1}=backDP;
%     backMod2{i,1}=backP3;
%     backMod2{i,1}=backP4;
%     NewcontrastMod{i,1}=modu;
%     NewcontrastMod{i,2}=datapoint;
%     NewcontrastMod2{i,1}=modu3;
%     NewbackMod{i,1}=backM;
%     NewbackMod{i,1}=backDP;
%     NewbackMod2{i,1}=backM3;

contrastResponse=cell(2,4);
contrastTogo=togo; abscnt=0; relIndex=[];
for i=1:length(togo)
    mod=contrastMod2{i,1};
    rate=contrastRates{i,1};
    brate=backRates{i,1};
    bmod=backMod2{i,1};
    int=contrastRates{i,2};
    
    if ~isempty(int)
        Y=fft(rate,8000)/numel(rate);
        f=1000/2*linspace(0,1,8000/2+1);
        c=2*abs(Y(1:8000/2+1));
        c=c';
        
        tmp1=find(f>1.2); tmp1=tmp1(1);
        tmp2=find(f<12); tmp2=tmp2(end);
        c=c(tmp1:tmp2);
        f=f(tmp1:tmp2);
        
        exclude=[]; pks=[];
        for fr=2:2:10
            tmp=find(f==fr);
            exclude=[exclude tmp-2:tmp+2];
            pks=[pks c(tmp)];
        end
        bamp=c;
        include=1:length(c);
        include=setdiff(include,exclude);
        bamp=bamp(include);
        bf=f(include);
        noise=mean(bamp);
        noiseSTD=std(bamp);
        tst=find(pks>noise+2*noiseSTD);
        
        if ~isempty(find(tst==1)) | ~isempty(find(tst==2))
            Y=fft(brate,8000)/numel(brate);
            f=1000/2*linspace(0,1,8000/2+1);
            c=2*abs(Y(1:8000/2+1));
            c=c';
            
            tmp1=find(f>1.2); tmp1=tmp1(1);
            tmp2=find(f<12); tmp2=tmp2(end);
            c=c(tmp1:tmp2);
            f=f(tmp1:tmp2);
            
            exclude=[]; pks=[];
            for fr=2:2:10
                tmp=find(f==fr);
                exclude=[exclude tmp-2:tmp+2];
                pks=[pks c(tmp)];
            end
            bamp=c;
            include=1:length(c);
            include=setdiff(include,exclude);
            bamp=bamp(include);
            bf=f(include);
            noise=mean(bamp);
            noiseSTD=std(bamp);
            tst=find(pks>noise+2*noiseSTD);
            if isempty(find(tst==1)) && isempty(find(tst==2))
                abscnt=abscnt+1;
                relIndex=[relIndex;i]; %relative to all matrices and togo
                %congrats, you might be a response!
                mback=mean(bmod);
                stdback=std(bmod);
                
                goodTime=find(mod>mback+stdback);
                response=mod(goodTime);
                goodContrast=int(goodTime);
                contrastResponse{abscnt,1}=goodTime;
                contrastResponse{abscnt,2}=response;
                contrastResponse{abscnt,3}=goodContrast;
                contrastResponse{abscnt,4}=mod-mean(bmod); %subtract background modulation
            else
                contrastTogo=setdiff(contrastTogo,togo(i));
            end
        else
            contrastTogo=setdiff(contrastTogo,togo(i));
        end
    else
        contrastTogo=setdiff(contrastTogo,togo(i));
    end
    
end
save(fullfile(superpath,'newAttempt','contrast',[Dtype,'_judgedContrast']),'contrastTogo','contrastResponse','togo')
%% plot contrast

clear
close all

superpath='f:\all_species_temp';

Dtype='wph';
cols='bgr';
version=1;

figure

for d=1:3
    load(fullfile(superpath,'newAttempt','contrast',[Dtype(d),'_judgedContrast']))
    coll=[]; av=[];
    for i=1:length(contrastResponse)
        %         contrastResponse{abscnt,1}=goodTime;
        %         contrastResponse{abscnt,2}=response;
        %         contrastResponse{abscnt,3}=goodContrast;
        
        time=contrastResponse{i,1};
        time=time(find(time<=8100));
        resp=contrastResponse{i,2};
        resp2=contrastResponse{i,4};
        fullResp=zeros(8100,1);
        fullResp(time)=resp./max(resp);
        %         fullResp2=resp2./max(resp2);
        fullResp2=resp2;
        %         fullResp(time)=resp;
        
        
        subplot(2,4,d)
        plot(fullResp,'color',[0.6 0.6 0.6])
        hold on
        coll=[coll;fullResp'];
        
        subplot(2,4,d+4)
        plot(fullResp2,'color',[0.6 0.6 0.6])
        hold on
        av=[av;fullResp2];
    end
    
    
    if version==2
        coll2=[]; std2=[];
        for c=1:size(coll,2)
            curr=coll(:,c);
            tmp=find(curr>0);
            coll2=[coll2 mean(curr(tmp))];
            std2=[std2 std(curr(tmp))];
        end
        subplot(2,4,d)
        plot(coll2,'-','color',cols(d))
        patch([1:length(coll2) length(coll2):-1:1],[coll2-std2 fliplr(coll2+std2)],cols(d),'facealpha',0.4)
        subplot(2,4,4)
        plot(coll2,'-','color',cols(d))
        hold on
    else
        coll1=mean(coll);
        std1=std(coll);
        subplot(2,4,d)
        plot(coll1,'-','color',cols(d))
        patch([1:length(coll1) length(coll1):-1:1],[coll1-std1 fliplr(coll1+std1)],cols(d),'facealpha',0.4)
        subplot(2,4,4)
        plot(coll1,'-','color',cols(d))
        hold on
    end
    
    avAll=median(av);
    stdAll=std(av);
    subplot(2,4,d+4)
    plot(avAll,'-','color',cols(d))
    hold on
    patch([1:length(avAll) length(avAll):-1:1],[avAll-stdAll fliplr(avAll+stdAll)],cols(d),'facealpha',0.4)
    subplot(2,4,8)
    plot(avAll,'-','color',cols(d))
    hold on
end
for s=1:8
    subplot(2,4,s)
    axis([0 8100 -inf inf])
end
%% ffflash firing rates lower sigma

load(fullfile(superpath,[Dtype,'_info']))
cd(superpath)

flash_rates=cell(length(info),1);
flash_raster=cell(length(info),10);
flash_ratesSmoo=cell(length(info),1);
flash_ratesSmoo2=cell(length(info),1);
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
    
    if ~isempty(regexp(Dtype,'w'))
        goodfiles=goodfiles(11:end);
    end
    
    allrate=[]; ffcnt=0;
    for f=goodfiles
        ffcnt=ffcnt+1;
        spikes=unit{1,2}{f,2};
        flash_raster{currID,ffcnt}=spikes;
        rate=MEA_spikerates(spikes,20,tt);
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
    
     allrate=[]; ffcnt=0;
    for f=goodfiles
        ffcnt=ffcnt+1;
        spikes=unit{1,2}{f,2};
        flash_raster{currID,ffcnt}=spikes;
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
    flash_ratesSmoo{currID,1}=mrate;
   
     allrate=[]; ffcnt=0;
    for f=goodfiles
        ffcnt=ffcnt+1;
        spikes=unit{1,2}{f,2};
        flash_raster{currID,ffcnt}=spikes;
        rate=MEA_spikerates(spikes,60,tt);
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
    flash_ratesSmoo2{currID,1}=mrate;
end

save(fullfile(superpath,[Dtype,'_flash20sigma.mat']),'flash_rates')
save(fullfile(superpath,[Dtype,'_flash40sigma.mat']),'flash_ratesSmoo')
save(fullfile(superpath,[Dtype,'_flash60sigma.mat']),'flash_ratesSmoo2')
save(fullfile(superpath,[Dtype,'_flashSpikes.mat']),'flash_raster')
%% calculate latency :(
clear
close all
path1='F:\all_species_temp';
path2='S:\data\Katja';
type='wph';
mouselist={'20150324b','20150424b','20150507a','20150507b'};
pathadd={'WT','Pig','Human'};

for t=1
    load(fullfile(path1,[type(t),'_flash20sigma']))
    load(fullfile(path1,[type(t),'_flashSpikes']))
    load(fullfile(path1,[type(t),'_info']))
    
    peaks=[]; lats=[];
    
    flash=[];  cnt=0; start=[];
    clear dates
    for i=1:length(info)
        if ~isempty(info{i,8})
            if info{i,8}~=0 && goodones(i)==1
                flash=[flash;i];
                cnt=cnt+1;
                dates{cnt}=info{i,1};
                start=[start;info{i,4}];
            end
        end
    end
    
    flashR=flash_rates(flash);
    flashS=flash_raster(flash,:);
    flashRtemp=cell(size(flashR,1),1);
    if t==1
        
        for f=1:length(flash)
            d=0; found=0;
            while found==0 && d<length(mouselist)
                d=d+1; found=0;
                if ~isempty(regexp(dates{f},mouselist{d}))
                    found=1;
                end
            end
            if found==1
                flashRtemp{f}=flashR{f}(17:end);
            else
                
                flashRtemp{f}=flashR{f}(33:end);
                
                
            end
            
        end
    else
        for f=1:length(flashR)
            
            flashRtemp{f}=flashR{f}(33:end);
            
        end
    end
    flashR=flashRtemp;
    
    for f=1:size(flashR,1)
        display([int2str(f) ,' out of ',int2str(size(flashR,1))])
        hekapath=fullfile(path2,pathadd{t},dates{f},'HEKA');
        hekalist=dir(fullfile(hekapath,'*phys'));
        m=start(f)-1; found=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(hekalist(m).name,'FFFlash'))
                found=1;
            end
        end
        prot=read_header_field_heka(hekapath,hekalist(m).name,'Stimulus Protocol');
        if prot(1)==0
            m=m+10; found=0;
            while found==0
                m=m+1;
                if ~isempty(regexp(hekalist(m).name,'FFFlash'))
                    found=1;
                end
            end
            prot=read_header_field_heka(hekapath,hekalist(m).name,'Stimulus Protocol');
        end
        prot=round(prot(:,1)); prot=prot(2:end-1);
        figure
        subplot(1,2,1)
        curspikes=flashS(f,:);
        
        for ce=1:size(flashS,2)
            if ~isempty(curspikes{ce})
                spikes=curspikes{ce};
                spikes=spikes-500;
                spikes=spikes(find(spikes>=0));
                spikes=spikes(find(spikes<length(flashR{f})-2000));
                for s=1:length(spikes)
                    plot([spikes(s) spikes(s)],[ce ce+0.8],'-')
                    hold on
                end
            end
        end
        for p=1:length(prot)
            plot([prot(p)-500 prot(p)-500],[0 size(flashS,2)+1],'-k')
        end
        axis([-inf inf 0 size(flashS,2)+1])
        
        subplot(1,2,2)
        plot(flashR{f}(500:end-2000));
        hold on
        for p=1:length(prot)
            plot([prot(p)-500 prot(p)-500],[0 max(flashR{f})],'-k')
        end
        
        %calculate latencies
        la=[]; pe=[];
        bckr=mean(flashR{f}(500:prot(1)));
        stback=std(flashR{f}(500:prot(1)));
        plot([0 length(flashR{f})-2500],[bckr+2*stback bckr+2*stback],'-g')
        for p=1:length(prot)-1
            curr=flashR{f}(prot(p):prot(p+1));
            curr=curr(80:end);
            [pp ll]=findpeaks(curr,'minpeakdistance',10);
            
            tmp=find(pp>bckr+2*stback); pp=pp(tmp); ll=ll(tmp);
            %             [pk loc]=max(curr);
            if ~isempty(pp)
                if length(pp)==1
                    pk=pp; loc=ll;
                else
                    in=[];
                    for j=1:length(pp)
                        mill=ll(j)-3; mall=ll(j)+3;
                        if mill<=0
                            mill=1;
                        end
                        if mall>length(curr)
                            mall=length(curr);
                        end
                        if curr(mill)<=pp(j) && curr(mall)<=pp(j)
                            in=[in j];
                        end
                        
                    end
                    pp=pp(in); ll=ll(in);
                    pp=pp(1:2);ll=ll(1:2);
                    
                    mi=min(curr(ll(1):ll(2)));
                    ma=min(pp);
                    if mi-bckr<(ma-bckr)*0.75
                        pk=pp(1); loc=ll(1);
                    elseif pp(2)>pp(1)
                        pk=pp(2); loc=ll(2);
                    else
                        pk=pp(1); loc=ll(1);
                    end
                end
                
                plot(prot(p)+79+loc-500,pk,'or')
                la=[la loc+79];
                pe=[pe pk];
            else
                la=[la 0];
                pe=[pe 0];
                
            end
        end
        lats=[lats;la];
        peaks=[peaks;pe];
        
        set(gcf,'position',[1601 1 1600 824])
        
        figure
        h = uibuttongroup('Position',[0 0 1 1]);
        % Create three radio buttons in the button group.
        u{1} = uicontrol('Style','checkbox','String','first',...
            'pos',[50 350 100 30],'parent',h,'HandleVisibility','off');
        u{2} = uicontrol('Style','checkbox','String','second',...
            'pos',[50 250 100 30],'parent',h,'HandleVisibility','off');
        u{3} = uicontrol('Style','checkbox','String','third',...
            'pos',[50 150 100 30],'parent',h,'HandleVisibility','off');
        u{4} = uicontrol('Style','checkbox','String','fourth',...
            'pos',[50 50 100 30],'parent',h,'HandleVisibility','off');
        ok=uicontrol('style','pushbutton','string','OK','position',  [150 200 50 50],'callback','uiresume');
        uiwait
        
        
        todo=[];
        for r=1:4
            tmp=get(u{r},'value');
            todo=[todo tmp];
        end
        todo=find(todo==1);
        
        if ~isempty(todo)
            figure(1)
            subplot(1,2,2)
            [x,y] = ginput(length(todo));
            peaks(f,todo)=y;
            lats(f,todo)=x+500-prot(todo);
            
            for tt=1:length(todo)
                if y(tt)>0
                    plot(x(tt),y(tt),'dm','markerfacecolor','m','markersize',10)
                end
                
            end
        end
        
        saveas(gcf,fullfile(path1,'latency',[type(t),'_',int2str(flash(f)),'.bmp']))
        close all
    end
    save(fullfile(path1,[type(t),'_latency']),'flash','flashR','peaks','lats')
end
%% calculate latency2 :(
clear
close all
path1='F:\all_species_temp';
path2='S:\data\Katja';
type='wph';
mouselist={'20150324b','20150424b','20150507a','20150507b'};
pathadd={'WT','Pig','Human'};

for t=1
    load(fullfile(path1,[type(t),'_flash60sigma']))
    load(fullfile(path1,[type(t),'_flashSpikes']))
    load(fullfile(path1,[type(t),'_info']))
    
  
    
    flash=[];  cnt=0; start=[];
    clear dates
    for i=1:length(info)
        if ~isempty(info{i,8})
            if info{i,8}~=0 && goodones(i)==1
                flash=[flash;i];
                cnt=cnt+1;
                dates{cnt}=info{i,1};
                start=[start;info{i,4}];
            end
        end
    end
    
    flashR=flash_ratesSmoo2(flash);
    flashS=flash_raster(flash,:);
    flashRtemp=cell(size(flashR,1),1);
    if t==1
        
        for f=1:length(flash)
            d=0; found=0;
            while found==0 && d<length(mouselist)
                d=d+1; found=0;
                if ~isempty(regexp(dates{f},mouselist{d}))
                    found=1;
                end
            end
            if found==1
                flashRtemp{f}=flashR{f}(17:end);
            else
                
                flashRtemp{f}=flashR{f}(33:end);
                
                
            end
            
        end
    else
        for f=1:length(flashR)
            
            flashRtemp{f}=flashR{f}(33:end);
            
        end
    end
    flashR=flashRtemp;
    peaks=zeros(size(flashR,1),4); lats=zeros(size(flashR,1),4);
    for f=1:size(flashR,1)
        display([int2str(f) ,' out of ',int2str(size(flashR,1))])
        hekapath=fullfile(path2,pathadd{t},dates{f},'HEKA');
        hekalist=dir(fullfile(hekapath,'*phys'));
        m=start(f)-1; found=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(hekalist(m).name,'FFFlash'))
                found=1;
            end
        end
        prot=read_header_field_heka(hekapath,hekalist(m).name,'Stimulus Protocol');
        if prot(1)==0
            m=m+10; found=0;
            while found==0
                m=m+1;
                if ~isempty(regexp(hekalist(m).name,'FFFlash'))
                    found=1;
                end
            end
            prot=read_header_field_heka(hekapath,hekalist(m).name,'Stimulus Protocol');
        end
        prot=round(prot(:,1)); prot=prot(2:end-1);
        figure
        subplot(1,2,1)
        curspikes=flashS(f,:);
        
        for ce=1:size(flashS,2)
            if ~isempty(curspikes{ce})
                spikes=curspikes{ce};
                spikes=spikes-500;
                spikes=spikes(find(spikes>=0));
                spikes=spikes(find(spikes<length(flashR{f})-2000));
                for s=1:length(spikes)
                    plot([spikes(s) spikes(s)],[ce ce+0.8],'-')
                    hold on
                end
            end
        end
        for p=1:length(prot)
            plot([prot(p)-500 prot(p)-500],[0 size(flashS,2)+1],'-k')
        end
        axis([-inf inf 0 size(flashS,2)+1])
        
        subplot(1,2,2)
        plot(flashR{f}(500:end-2000));
        hold on
        for p=1:length(prot)
            plot([prot(p)-500 prot(p)-500],[0 max(flashR{f})],'-k')
        end
        
        %calculate latencies
        la=[]; pe=[];
        bckr=mean(flashR{f}(500:prot(1)));
        stback=std(flashR{f}(500:prot(1)));
        plot([0 length(flashR{f})-2500],[bckr+2*stback bckr+2*stback],'-g')
        for p=1:length(prot)-1
            curr=flashR{f}(prot(p):prot(p+1));
            curr=curr(80:end);
            [pp ll]=findpeaks(curr,'minpeakdistance',10);
            
            tmp=find(pp>bckr+2*stback); pp=pp(tmp); ll=ll(tmp);
            %             [pk loc]=max(curr);
            if ~isempty(pp)
                if length(pp)==1
                    pk=pp; loc=ll;
                else
                    in=[];
                    for j=1:length(pp)
                        mill=ll(j)-3; mall=ll(j)+3;
                        if mill<=0
                            mill=1;
                        end
                        if mall>length(curr)
                            mall=length(curr);
                        end
                        if curr(mill)<=pp(j) && curr(mall)<=pp(j)
                            in=[in j];
                        end
                        
                    end
                    if length(in)>1
                    pp=pp(in); ll=ll(in);
                    pp=pp(1:2);ll=ll(1:2);
                    
                    mi=min(curr(ll(1):ll(2)));
                    ma=min(pp);
                    if mi-bckr<(ma-bckr)*0.75
                        pk=pp(1); loc=ll(1);
                    elseif pp(2)>pp(1)
                        pk=pp(2); loc=ll(2);
                    else
                        pk=pp(1); loc=ll(1);
                    end
                    else
                        pk=pp(in); loc=ll(in); 
                    end
                end
                pp=pk; ll=loc;
            else
               ll=0; pp=0;
                
            end
          
            if pp>=bckr+2*stback && pp>0
                peak75=find(curr(1:ll)-bckr>=(pp-bckr)*0.75);
                if peak75(1)==1
                    peak75=peak75(2:end);
                end
                tmp=find(curr(peak75-1)<curr(peak75));
                peak75=peak75(tmp);
                ll75=peak75(1);
                pp75=curr(ll75);
                plot(ll+79+prot(p)-500,pp,'ok','markerfacecolor','k','markersize',8)
                plot(ll75+79+prot(p)-500,pp75,'dr','markerfacecolor','r')
                la=[la ll75+79]; pe=[pe pp75-bckr];
                
            else
                la=[la 0];
                pe=[pe 0];
                
            end
        end
        lats(f,:)=la;
        peaks(f,:)=pe;
        
        set(gcf,'position',[1601 1 1600 824])
        
        figure
        h = uibuttongroup('Position',[0 0 1 1]);
        % Create three radio buttons in the button group.
        u{1} = uicontrol('Style','checkbox','String','first',...
            'pos',[50 350 100 30],'parent',h,'HandleVisibility','off');
        u{2} = uicontrol('Style','checkbox','String','second',...
            'pos',[50 250 100 30],'parent',h,'HandleVisibility','off');
        u{3} = uicontrol('Style','checkbox','String','third',...
            'pos',[50 150 100 30],'parent',h,'HandleVisibility','off');
        u{4} = uicontrol('Style','checkbox','String','fourth',...
            'pos',[50 50 100 30],'parent',h,'HandleVisibility','off');
         h2 = uibuttongroup('Position',[0 0 1 1]);
        % Create three radio buttons in the button group.
        v{1} = uicontrol('Style','checkbox','String','first',...
            'pos',[150 350 100 30],'parent',h,'HandleVisibility','off');
        v{2} = uicontrol('Style','checkbox','String','second',...
            'pos',[150 250 100 30],'parent',h,'HandleVisibility','off');
        v{3} = uicontrol('Style','checkbox','String','third',...
            'pos',[150 150 100 30],'parent',h,'HandleVisibility','off');
        v{4} = uicontrol('Style','checkbox','String','fourth',...
            'pos',[150 50 100 30],'parent',h,'HandleVisibility','off');
        ok=uicontrol('style','pushbutton','string','OK','position',  [250 200 50 50],'callback','uiresume');
        uiwait
        
        
        todo=[];
        for r=1:4
            tmp=get(u{r},'value');
            todo=[todo tmp];
        end
        todo=find(todo==1);
        
         todo2=[];
        for r=1:4
            tmp=get(v{r},'value');
            todo2=[todo2 tmp];
        end
        todo2=find(todo2==1);
        
        if ~isempty(todo)
            figure(1)
            subplot(1,2,2)
            [x,y] = ginput(length(todo));
            if ~isempty(todo2)
                hh = msgbox('now borders');
                figure(1)
            subplot(1,2,2)
            [xx,yy] = ginput(length(todo));
            end
            
            for tt=1:length(todo)
                if y(tt)>0
                    curr=flashR{f}(prot(todo(tt)):prot(todo(tt)+1));
                    curr=curr(80:end);
                    pp=y(tt); ll=x(tt)-79-prot(todo(tt))+500;
                    tst=find(todo2==todo(tt));
                    if ~isempty(tst)
                        starter=round(xx(tst)-prot(todo(tt))+500-80);
                    else
                        starter=1;
                    end
                    peak75=find(curr(starter:ll)-bckr>=(pp-bckr)*0.75);
                    if peak75(1)==1
                        peak75=peak75(2:end);
                    end
                     tmp=find(curr(peak75-1)<curr(peak75));
            peak75=peak75(tmp);
                ll75=peak75(1)+starter;
                    pp75=curr(ll75);
                    
                    plot(ll75+79+prot(todo(tt))-500,pp75,'^r','markerfacecolor','r','markersize',8)
                    lats(f,todo(tt))=ll75+79;
                    peaks(f,todo(tt))=pp75-bckr;
                else
                    text(prot(todo(tt))+200,max(flashR{f})-4,'none')
                    lats(f,todo(tt))=0;
                    peaks(f,todo(tt))=0;
                    
                end
                
            end
        end
        
        pause(1)
        saveas(gcf,fullfile(path1,'latency2',[type(t),'_',int2str(flash(f)),'.bmp']))
        close all
        if ~isempty(todo2)
        close(hh)
        end
    end
    save(fullfile(path1,[type(t),'_latency2']),'flash','flashR','peaks','lats')
end
%% calculate latency: check all delayed responses
clear
close all
path1='F:\all_species_temp';
path2='S:\data\Katja';
type='wph';
mouselist={'20150324b','20150424b','20150507a','20150507b'};
pathadd={'WT','Pig','Human'};

t=1;

load(fullfile(path1,[type(t),'_flash60sigma']))
load(fullfile(path1,[type(t),'_flashSpikes']))
load(fullfile(path1,[type(t),'_info']))
load(fullfile(path1,[type(t),'_latency2']))



for f=1:size(flashR,1)
    display([int2str(f) ,' out of ',int2str(size(flashR,1))])
    date=info{flash(f),1};
    hekapath=fullfile(path2,pathadd{t},date,'HEKA');
    hekalist=dir(fullfile(hekapath,'*phys'));
    
    start=info{flash(f),4};
    
    m=start-1; found=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(hekalist(m).name,'FFFlash'))
            found=1;
        end
    end
    prot=read_header_field_heka(hekapath,hekalist(m).name,'Stimulus Protocol');
    if prot(1)==0
        m=m+10; found=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(hekalist(m).name,'FFFlash'))
                found=1;
            end
        end
        prot=read_header_field_heka(hekapath,hekalist(m).name,'Stimulus Protocol');
    end
    prot=round(prot(:,1)); prot=prot(2:end-1);
    figure
    subplot(1,2,1)
    flashS=flash_raster(flash,:);
    curspikes=flashS(f,:);
    
    for ce=1:size(flashS,2)
        if ~isempty(curspikes{ce})
            spikes=curspikes{ce};
            spikes=spikes-500;
            spikes=spikes(find(spikes>=0));
            spikes=spikes(find(spikes<length(flashR{f})-2000));
            for s=1:length(spikes)
                plot([spikes(s) spikes(s)],[ce ce+0.8],'-')
                hold on
            end
        end
    end
    for p=1:length(prot)
        plot([prot(p)-500 prot(p)-500],[0 size(flashS,2)+1],'-k')
    end
    axis([-inf inf 0 size(flashS,2)+1])
    
    subplot(1,2,2)
    plot(flashR{f}(500:end-2000));
    hold on
    for p=1:length(prot)
        plot([prot(p)-500 prot(p)-500],[0 max(flashR{f})],'-k')
    end
    
    %calculate latencies
    currLat=lats(f,:);
    tmp=find(currLat>0);
    tmp2=find(currLat(tmp)<1000);
    condi=tmp(tmp2);
    
    if ~isempty(condi)
        la=zeros(1,4); pe=zeros(1,4);
        bckr=mean(flashR{f}(500:prot(1)));
        stback=std(flashR{f}(500:prot(1)));
        plot([0 length(flashR{f})-2500],[bckr+2*stback bckr+2*stback],'-g')
        for p=condi
            curr=flashR{f}(prot(p):prot(p+1));
            curr=curr(currLat(p)+120:end);
            [pp ll]=findpeaks(curr,'minpeakdistance',10);
            
            tmp=find(pp>bckr+2*stback); pp=pp(tmp); ll=ll(tmp);
            %             [pk loc]=max(curr);
            if ~isempty(pp)
                if length(pp)==1
                    pk=pp; loc=ll;
                else
                    in=[];
                    for j=1:length(pp)
                        mill=ll(j)-3; mall=ll(j)+3;
                        if mill<=0
                            mill=1;
                        end
                        if mall>length(curr)
                            mall=length(curr);
                        end
                        if curr(mill)<=pp(j) && curr(mall)<=pp(j)
                            in=[in j];
                        end
                        
                    end
                    if length(in)>1
                        pp=pp(in); ll=ll(in);
                        pp=pp(1:2);ll=ll(1:2);
                        
                        mi=min(curr(ll(1):ll(2)));
                        ma=min(pp);
                        if mi-bckr<(ma-bckr)*0.75
                            pk=pp(1); loc=ll(1);
                        elseif pp(2)>pp(1)
                            pk=pp(2); loc=ll(2);
                        else
                            pk=pp(1); loc=ll(1);
                        end
                    else
                        pk=pp(in); loc=ll(in);
                    end
                end
                pp=pk; ll=loc;
            else
                ll=0; pp=0;
                
            end
            
            if pp>=bckr+2*stback && pp>0
                peak75=find(curr(1:ll)-bckr>=(pp-bckr)*0.75);
                if peak75(1)==1
                    peak75=peak75(2:end);
                end
                tmp=find(curr(peak75-1)<curr(peak75));
                peak75=peak75(tmp);
                ll75=peak75(1);
                pp75=curr(ll75);
                plot(ll+currLat(p)+119+prot(p)-500,pp,'ok','markerfacecolor','k','markersize',8)
                plot(ll75+currLat(p)+119+prot(p)-500,pp75,'dr','markerfacecolor','r')
                la(p)=ll75+currLat(p)+119; pe(p)=pp75-bckr;
                
            else
                la(p)=0;
                pe(p)=0;
                
            end
        end
        lats(f,:)=la;
        peaks(f,:)=pe;
        
        set(gcf,'position',[1601 1 1600 824])
        
        
        
        figure
        h = uibuttongroup('Position',[0 0 1 1]);
        % Create three radio buttons in the button group.
        u{1} = uicontrol('Style','checkbox','String','first',...
            'pos',[50 350 100 30],'parent',h,'HandleVisibility','off');
        u{2} = uicontrol('Style','checkbox','String','second',...
            'pos',[50 250 100 30],'parent',h,'HandleVisibility','off');
        u{3} = uicontrol('Style','checkbox','String','third',...
            'pos',[50 150 100 30],'parent',h,'HandleVisibility','off');
        u{4} = uicontrol('Style','checkbox','String','fourth',...
            'pos',[50 50 100 30],'parent',h,'HandleVisibility','off');
        h2 = uibuttongroup('Position',[0 0 1 1]);
        % Create three radio buttons in the button group.
        v{1} = uicontrol('Style','checkbox','String','first',...
            'pos',[150 350 100 30],'parent',h,'HandleVisibility','off');
        v{2} = uicontrol('Style','checkbox','String','second',...
            'pos',[150 250 100 30],'parent',h,'HandleVisibility','off');
        v{3} = uicontrol('Style','checkbox','String','third',...
            'pos',[150 150 100 30],'parent',h,'HandleVisibility','off');
        v{4} = uicontrol('Style','checkbox','String','fourth',...
            'pos',[150 50 100 30],'parent',h,'HandleVisibility','off');
        ok=uicontrol('style','pushbutton','string','OK','position',  [250 200 50 50],'callback','uiresume');
        uiwait
        
        
        todo=[];
        for r=1:4
            tmp=get(u{r},'value');
            todo=[todo tmp];
        end
        todo=find(todo==1);
        
        todo2=[];
        for r=1:4
            tmp=get(v{r},'value');
            todo2=[todo2 tmp];
        end
        todo2=find(todo2==1);
        
        if ~isempty(todo)
            figure(1)
            subplot(1,2,2)
            [x,y] = ginput(length(todo));
            if ~isempty(todo2)
                hh = msgbox('now borders');
                figure(1)
                subplot(1,2,2)
                [xx,yy] = ginput(length(todo));
            end
            
            for tt=1:length(todo)
                if y(tt)>0
                    curr=flashR{f}(prot(todo(tt)):prot(todo(tt)+1));
                    curr=curr(currLat(tt)+120:end);
                    pp=y(tt); ll=x(tt)-currLat(tt)-119-prot(todo(tt))+500;
                    tst=find(todo2==todo(tt));
                    if ~isempty(tst)
                        starter=round(xx(tst)-prot(todo(tt))+500-119-currLat(tt));
                    else
                        starter=1;
                    end
                    peak75=find(curr(starter:ll)-bckr>=(pp-bckr)*0.75);
                    if peak75(1)==1
                        peak75=peak75(2:end);
                    end
                    tmp=find(curr(peak75-1)<curr(peak75));
                    peak75=peak75(tmp);
                    ll75=peak75(1)+starter;
                    pp75=curr(ll75);
                    
                    plot(ll75+currLat(tt)+119+prot(todo(tt))-500,pp75,'^r','markerfacecolor','r','markersize',8)
                    lats(f,todo(tt))=ll75+119+currLat(tt);
                    peaks(f,todo(tt))=pp75-bckr;
                else
                    text(prot(todo(tt))+200,max(flashR{f})-4,'none')
                    lats(f,todo(tt))=0;
                    peaks(f,todo(tt))=0;
                    
                end
                
            end
        end
        
        pause(1)
        saveas(gcf,fullfile(path1,'latency_del',[type(t),'_',int2str(flash(f)),'.bmp']))
        close all
        if ~isempty(todo2)
            close(hh)
        end
    else
        lats(f,:)=[0 0 0 0];
        peaks(f,:)=[0 0 0 0];
    end
end
%     save(fullfile(path1,[type(t),'_latency2']),'flash','flashR','peaks','lats')

%% remove empty lines
su=sum(lats,2);
tmp=find(su==0);
all=1:length(flash);
keep=setdiff(all,tmp);
flash=flash(keep);
lats=lats(keep,:);
peaks=peaks(keep,:);
flashR=flashR(keep);
%% plotting latencies
close all
clear

path1='F:\all_species_temp';
list='wph';
cols='bgr';
allLat=[]; type=[]; flashID=[];
for l=1:length(list);
    load(fullfile(path1,[list(l),'_latency2']))
    allLat=[allLat;lats];
    type=[type;zeros(length(lats),1)+l];
    flashID=[flashID;flash];
end

path2='F:\all_species_temp\final_pics';
figure
for i=1:4
    tmp=find(allLat(:,i)>0);
    data=allLat(tmp,i);
    me=median(data);
    st=std(data);
    data=[0;data;2000];
    
   subplot(2,2,i)
    [h int]=hist(data,40);
    bar(int(2:end-1),100/length(data)*h(2:end-1))
    hold on
    plot(me,max(100/length(data)*h(2:end-1))+2,'or')
    plot([me-st me+st],[max(100/length(data)*h(2:end-1))+2 max(100/length(data)*h(2:end-1))+2],'-r')
    axis([-100 2000 0 max(100/length(data)*h(2:end-1))+4])
title(['latency all species for step ',int2str(i),'; median=',num2str(round(me))])
ylabel('% of cells with flash-responses')
xlabel('latency [ms]')
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'latency_allSpecies_perCond.emf'))
close

close all
figure(1)
compl=reshape(allLat,length(allLat)*4,1);
tmp=find(compl>0);
data=compl(tmp);
me=median(data);
st=std(data);
data=[0;data;2000];
subplot(2,2,1)
figure
h=histfit(data,40,'kernel');
figure(1)
% [h int]=hist(data,40);
int=get(h(1),'xdata');
int=mean(int(2:3,:));
hi=get(h(1),'ydata');
hi=hi(2,:);
xfit=get(h(2),'xdata');
yfit=get(h(2),'ydata');
bar(int(2:end-1),100/length(data-2)*hi(2:end-1))
hold on
plot(xfit,100/length(data-2)*yfit,'-')
plot(me,max(100/length(data-2)*hi(2:end-1))+2,'or')
plot([me-st me+st],[max(100/length(data-2)*hi(2:end-1))+2 max(100/length(data-2)*hi(2:end-1))+2],'-r')
axis([-100 2000 0 max(100/length(data-2)*hi(2:end-1))+4])
title(['all','; median=',num2str(round(me))])
ylabel('% of cells with flash-responses')
xlabel('latency [ms]')

for t=1:length(list)
    tst=find(type==t);
    compl=reshape(allLat(tst,:),length(tst)*4,1);
    tmp=find(compl>0);
    data=compl(tmp);
    me=median(data);
    st=std(data);
    
    data=[0;data;2000];
    figure
h=histfit(data,40,'kernel');
figure(1)

    subplot(2,2,t+1)

   int=get(h(1),'xdata');
   keepint=reshape(int(1:2,:),length(int)*2,1);
keepint=sort(keepint);
keepint=[keepint;keepint(end)+diff(keepint(end-1:end));keepint(end)+diff(keepint(end-1:end))];
int=mean(int(2:3,:));
hi=get(h(1),'ydata');
keephi=[];
for x=1:length(hi)
    keephi=[keephi;hi(2:3,x)];
end
keephi(1:2)=0;
keephi(end-1:end)=0;
keephi=[0;keephi;0];
keephi=100/length(data-2)*keephi;
hi=hi(2,:);
xfit=get(h(2),'xdata');
yfit=get(h(2),'ydata');
bar(int(2:end-1),100/length(data-2)*hi(2:end-1))
hold on
plot(xfit,100/length(data-2)*yfit,'-')
plot(me,max(100/length(data-2)*hi(2:end-1))+2,'or')
plot([me-st me+st],[max(100/length(data-2)*hi(2:end-1))+2 max(100/length(data-2)*hi(2:end-1))+2],'-r')
axis([-100 2000 0 max(max(100/length(data-2)*yfit),max(100/length(data-2)*hi(2:end-1)))+4])
    title([list(t),' (all conditions)','; median=',num2str(round(me))])
    ylabel('% of cells with flash-responses')
    xlabel('latency [ms]')
    figure(100)
    plot(keepint,keephi,cols(t))
    hold on
end
figure(1)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'latency_acrossCond.emf'))
close
figure(100)
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'latency_acrossCond2.emf'))
close

for t=1:3
    figure
    for i=1:4
        tmp=find(type==t);
        la=allLat(tmp,:);
        tmp=find(la(:,i)>0);
        data=la(tmp,i);
        me=median(data);
        st=std(data);
        data=[0;data;2000];
        subplot(2,2,i)
        [h int]=hist(data,50);
        bar(int(2:end-1),100/length(data)*h(2:end-1))
        hold on
        plot(me,max(100/length(data)*h(2:end-1))+2,'or')
        plot([me-st me+st],[max(100/length(data)*h(2:end-1))+2 max(100/length(data)*h(2:end-1))+2],'-r')
        axis([-100 2000 0 max(100/length(data)*h(2:end-1))+4])
        title([list(t),', condition ',int2str(i),'; median=',num2str(round(me))])
        ylabel('% of cells with flash-responses')
        xlabel('latency [ms]')
    end
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile(path2,['latency_',list(t),'_perCond.emf']))
close
end

early=zeros(size(allLat,1),size(allLat,2));
all=zeros(size(allLat,1),size(allLat,2));
for i=1:4
    tmp=find(allLat(:,i)<=550);
    tmp2=find(allLat(tmp,i)>0);
early(tmp(tmp2),i)=1;
tmp=find(allLat(:,i)>0);
all(tmp,i)=1;
end

onoff=zeros(size(allLat,1),1);
anyresp=zeros(size(allLat,1),1);
for i=1:length(early)
    if isempty(find(early(i,:)>0));
        onoff(i)=0;
    elseif ~isempty(find(early(i,[1 4])==1)) && isempty(find(early(i,[2 3])==1))
        onoff(i)=-1;
    elseif isempty(find(early(i,[1 4])==1)) && ~isempty(find(early(i,[2 3])==1))
        onoff(i)=1;
    else
        onoff(i)=2;
    end
    
     if isempty(find(all(i,:)>0));
        anyresp(i)=0;
    elseif ~isempty(find(all(i,[1 4])==1)) && isempty(find(all(i,[2 3])==1))
        anyresp(i)=-1;
    elseif isempty(find(all(i,[1 4])==1)) && ~isempty(find(all(i,[2 3])==1))
        anyresp(i)=1;
    else
        anyresp(i)=2;
    end
end

figure
re=[-1 1 2 0];
for t=1:length(list)
    tst=find(type==t);
    num=[];
    for r=1:length(re)
    tmp=find(onoff(tst)==re(r));
    num=[num length(tmp)];
    end
    subplot(2,2,t)
    h=pie(num)
    hText = findobj(h,'Type','text'); % text handles
percentValues = get(hText,'String'); % percent values
str = {'off ';'on ';'onoff ';'none '}; % text
tmp=find(num>0);
combinedstrings = strcat(str(tmp),percentValues);
set(hText,{'String'},combinedstrings);
title(['polartiy based on early responses (<=550ms) of ',list(t),' (total: ',int2str(length(tst)),' cells)'])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'polarity_early.emf'))
close

figure
re=[-1 1 2];
for t=1:length(list)
    tst=find(type==t);
    num=[];
    for r=1:length(re)
    tmp=find(anyresp(tst)==re(r));
    num=[num length(tmp)];
    end
    subplot(2,2,t)
    h=pie(num)
    hText = findobj(h,'Type','text'); % text handles
percentValues = get(hText,'String'); % percent values
str = {'off ';'on ';'onoff ';'none '}; % text
tmp=find(num>0);
combinedstrings = strcat(str(tmp),percentValues);
set(hText,{'String'},combinedstrings);
title(['polartiy based on all responses of ',list(t),' (total: ',int2str(length(tst)),' cells)'])
end
set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(path2,'polarity_all.emf'))
close
    
%% ...just a short tester
figure
for i=1:4
    tmp=find(lats(:,i)>0);
    tmp2=find(peaks(tmp,i)>0);
    tmp=tmp(tmp2);
    data=lats(tmp,i);
    data=[0;data;2000];
   subplot(2,2,i)
    [h int]=hist(data,50);
    bar(int(2:end-1),100/length(data)*h(2:end-1))

end

%% --- here some good overview plots etc ---
%% POLARITY: check from DS stimulus

velTiming=[1001 20525 26811 31463 35264 38649 41833];
DStiming=[1001 6053 11121 15406 19708 24927 30162 35381 42000];


load(fullfile(superpath,[Dtype,'_DSblack']))
load(fullfile(superpath,[Dtype,'_DSwhite']))

timing=cell(length(DS_ratesB),2);
for i=1:length(DS_ratesB)
    if ~isempty(DS_ratesB{i}) && ~isempty(DS_ratesW{i})
        rateB=DS_ratesB{i};
        rateW=DS_ratesW{i};
        ps=zeros(8,2);
        for t=1:length(DStiming)-1
            currB=rateB(DStiming(t):DStiming(t+1));
            currW=rateW(DStiming(t):DStiming(t+1));
            [pB tB]=max(currB);
            [pW tW]=max(currW);
            ps(t,1)=tB;
            ps(t,2)=tW;
        end
        timing{i}=ps;
    end
end

pol=zeros(length(timing),1);
for t=1:length(timing)
    curr=timing{t};
    if ~isempty(curr)
        tmp=find(curr(:,1)>curr(:,2));
        if length(tmp)>5
            pol(t,1)=-1;
        elseif length(tmp)<3
            pol(t,1)=1;
        end
    end
end
load(fullfile(superpath,[Dtype,'_info']))
good=find(goodones==1);
flash=[];
for i=1:length(good)
    if ~isempty(info{good(i),8}) && isempty(find(info{good(i),8}==0))
        flash=[flash;good(i)];
    end
end

ds=[];
for i=1:length(good)
    if ~isempty(info{good(i),10})
        ds=[ds;good(i)];
    end
end
lf=[];
for i=1:length(good)
    if ~isempty(info{good(i),9})
        lf=[lf;good(i)];
    end
end
lfPol=cell2mat(info(lf,9));
% doublegood=intersect(ds,flash);

flashPol=cell2mat(info(flash,8));

cols='grb';
figure
overDS=[]; overLF=[];
for f=1:length(flashPol)
    id=flash(f);
    id2=find(ds==id);
    if ~isempty(id2)
        p1=flashPol(f);
        p2=pol(id2);
        if p1==1 && p2==1
            ok=1;
        elseif p1==-1 && p2==-1
            ok=1;
        elseif p1==2 && p2==-1
            ok=1;
        elseif p2==0
            ok=3;
        else
            ok=2;
            
        end
        plot(f,2,'o','markerfacecolor',cols(ok),'markeredgecolor',cols(ok))
        hold on
        overDS=[overDS; ok];
    else
        overDS=[overDS;-1];
    end
    id2=find(lf==id);
    if ~isempty(id2)
        p2=lfPol(id2);
        if p1==1 && p2==1
            ok=1;
        elseif p1==-1 && p2==-1
            ok=1;
        elseif p1==2 && p2==-1
            ok=1;
        elseif p2==0
            ok=3;
        else
            ok=2;
            
        end
        plot(f,1,'o','markerfacecolor',cols(ok),'markeredgecolor',cols(ok))
        hold on
        overLF=[overLF; ok];
    else
        overLF=[overLF;-1];
    end
end

tmp=find(overLF==2); %number of cells where flash and LF are contradictory
test=overDS(tmp);
tocheck=flash(tmp);
tocheck2=find(test==2);
tocheck2=tocheck(tocheck2);
dscells=[];
for ii=1:length(info)
    if info{ii,12}==1;
        dscells=[dscells;ii];
    end
end
tocheck3=setdiff(tocheck2,dscells); %number of contradictory cells where DS is also contradictory (and cell isn't DS)
%% POLARITY: check from 6vel stimulus


velTiming=[1001 20525 26811 31463 35264 38649 41833];
DStiming=[1001 6053 11121 15406 19708 24927 30162 35381 42000];


load(fullfile(superpath,[Dtype,'_vel6_new_Black']))
load(fullfile(superpath,[Dtype,'_vel6_new_White']))

timing=cell(length(vel6_newB),2);
for i=1:length(vel6_newB)
    if ~isempty(vel6_newB{i}) && ~isempty(vel6_newW{i})
        rateB=vel6_newB{i};
        rateW=vel6_newW{i};
        ps=zeros(8,2);
        for t=1:length(velTiming)-1
            currB=rateB(DStiming(t):DStiming(t+1));
            currW=rateW(DStiming(t):DStiming(t+1));
            [pB tB]=max(currB);
            [pW tW]=max(currW);
            ps(t,1)=tB;
            ps(t,2)=tW;
        end
        timing{i}=ps;
    end
end

pol=zeros(length(timing),1);
for t=1:length(timing)
    curr=timing{t};
    if ~isempty(curr)
        tmp=find(curr(:,1)>curr(:,2));
        if length(tmp)>4
            pol(t,1)=-1;
        elseif length(tmp)<3
            pol(t,1)=1;
        end
    end
end
load(fullfile(superpath,[Dtype,'_info']))
good=find(goodones==1);
flash=[];
for i=1:length(good)
    if ~isempty(info{good(i),8}) && isempty(find(info{good(i),8}==0))
        flash=[flash;good(i)];
    end
end

ds=[];
for i=1:length(good)
    if ~isempty(info{good(i),16})
        ds=[ds;good(i)];
    end
end
lf=[];
for i=1:length(good)
    if ~isempty(info{good(i),9})
        lf=[lf;good(i)];
    end
end
lfPol=cell2mat(info(lf,9));
% doublegood=intersect(ds,flash);

flashPol=cell2mat(info(flash,8));

cols='grb';
figure
overDS=[]; overLF=[];
for f=1:length(flashPol)
    
    id=flash(f);
    id2=find(ds==id);
    if ~isempty(id2)
        p1=flashPol(f);
        p2=pol(id2);
        if p1==1 && p2==1
            ok=1;
        elseif p1==-1 && p2==-1
            ok=1;
        elseif p1==2 && p2==-1
            ok=1;
        elseif p2==0
            ok=3;
        else
            ok=2;
            
        end
        plot(f,2,'o','markerfacecolor',cols(ok),'markeredgecolor',cols(ok))
        hold on
        overDS=[overDS; ok];
    else
        overDS=[overDS;-1];
    end
    id2=find(lf==id);
    if ~isempty(id2)
        p1=flashPol(f);
        p2=lfPol(id2);
        if p1==1 && p2==1
            ok=1;
        elseif p1==-1 && p2==-1
            ok=1;
        elseif p1==2 && p2==-1
            ok=1;
        elseif p2==0
            ok=3;
        else
            ok=2;
            
        end
        plot(f,1,'o','markerfacecolor',cols(ok),'markeredgecolor',cols(ok))
        hold on
        overLF=[overLF; ok];
    else
        overLF=[overLF;-1];
    end
end

tmp=find(overLF==2); %number of cells where flash and LF are contradictory
test=overDS(tmp);
tocheck=flash(tmp);
tocheck2=find(test==2);
tocheck2=tocheck(tocheck2);
% dscells=[];
% for ii=1:length(info)
%     if info{ii,12}==1;
%         dscells=[dscells;ii];
%     end
% end
% tocheck3=setdiff(tocheck2,dscells);
%% NICE PLOTS: DG normPeaks
% plot for each temporal/spatial frequency the spatial/temporal tuning
% curve based on normalized F1 peaks (normalized for each cell to its max.
% F1)
% blue=mouse, green=pig, red=human
% including STD (thin line), and 95%CI (thick line)
% version1: take always all cells which have been classified as
% "DG-positive" (at least 1 response to one DG combination)
% version2: take for each single compared value only the cells which
% respond
clear
close all
superpath='f:\all_species_temp';
version=1;

col1={'bbcc','ggkk','rrmm'};
col2={'bbbccc','gggkkk','rrrmmm'};
style1={'-','--','-','--'};
style2={'-','--','-.','-','--','-.'};
cols='bgr';
Dtype='wph';
spatials=[100 200 500 1000 2000 4000];
freqs=[1 2 4 8];
harmlist=[0.5 1 2 3];

collNorm=cell(3,1);
NofCT=zeros(3,6,4);
NofCS=zeros(3,4,6);
goodT=cell(3,6,4);
goodS=cell(3,4,6);
% NofC=zeros(3,24);
meanPeakT=zeros(3,6,4);
meanPeakS=zeros(3,4,6);
stT=zeros(3,6,4);
stS=zeros(3,4,6);
ciT=zeros(3,6,4);
ciS=zeros(3,4,6);
for d=1:3
    fcnt=0;
    load(fullfile(superpath,'newAttempt2','DG',[Dtype(d),'_normPeaks']))
    collNorm{d,1}=reshape(normPeaks(:,:,2),size(normPeaks,1),24);
    scnt=0;
    for s=1:4:24
        scnt=scnt+1;
        
        fcnt=fcnt+1;
        figure(fcnt)
        
        if version==2
            meanPeak=[]; st=[]; ci=[]; gcnt=0;
            for g=s:s+3
                gcnt=gcnt+1;
                good=find(normPeaks(:,g,2)>0);
                goodT{d,scnt,gcnt}=good;
                
                meanPeak=[meanPeak mean(normPeaks(good,g,2))];
                st=[st std(normPeaks(good,g,2))];
                ci=[ci 1.96*st(end)/sqrt(length(good))];
                NofCT(d,scnt,gcnt)=length(good);
            end
        else
            meanPeak=mean(normPeaks(:,s:s+3,2));
            st=std(normPeaks(:,s:s+3,2));
            ci=1.96*st/sqrt(size(normPeaks,1));
            for g=1:4
                NofCT(d,scnt,g)=size(normPeaks,1);
                goodT{d,scnt,g}=1:size(normPeaks,1);
            end
        end
        
        %         st=std(normPeaks(:,s:s+3,2));
        %         ci=1.96*st/sqrt(size(normPeaks,1));
        %             normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
        plot((freqs)+(d-1)*0.1,meanPeak,'color',cols(d),'linewidth',2)
        hold on
        for q=1:4
            plot([freqs(q)+(d-1)*0.1 freqs(q)+(d-1)*0.1],[meanPeak(q)-st(q) meanPeak(q)+st(q)],'color',cols(d))
            plot([freqs(q)+(d-1)*0.1 freqs(q)+(d-1)*0.1],[meanPeak(q)-ci(q) meanPeak(q)+ci(q)],'color',cols(d),'linewidth',5)
        end
        %         set(gca,'xtick',1:4)
        %         set(gca,'xticklabel',freqs)
        xlabel('frequency [Hz]')
        axis([0 9 -inf inf])
        title([int2str(spatials(scnt)),'um / f',num2str(harmlist(2))])
        
        
        meanPeakT(d,scnt,:)=meanPeak;
        stT(d,scnt,:)=st;
        ciT(d,scnt,:)=ci;
        
    end
    scnt=0;
    for s=1:4
        fcnt=fcnt+1;
        figure(fcnt)
        scnt=scnt+1;
        
        if version==2
            meanPeak=[]; st=[]; ci=[]; gcnt=0;
            for g=s:4:24
                gcnt=gcnt+1;
                good=find(normPeaks(:,g,2)>0);
                goodS{d,scnt,gcnt}=good;
                
                meanPeak=[meanPeak mean(normPeaks(good,g,2))];
                st=[st std(normPeaks(good,g,2))];
                ci=[ci 1.96*st(end)/sqrt(length(good))];
                NofCS(d,scnt,gcnt)=length(good);
            end
        else
            meanPeak=mean(normPeaks(:,s:4:24,2));
            st=std(normPeaks(:,s:4:24,2));
            ci=1.96*st/sqrt(size(normPeaks,1));
            for g=1:6
                NofCS(d,scnt,g)=size(normPeaks,1);
                goodS{d,scnt,g}=1:size(normPeaks,1);
            end
        end
        
        %         meanPeak=mean(normPeaks(:,s:4:24,2));%contains peaks normalized to each cell's f1
        %         st=std(normPeaks(:,s:4:24,2));
        %                 ci=1.96*st/sqrt(size(normPeaks,1));
        
        %             normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
        semilogx(([100 200 500 1000 2000 4000])+spatials*0.1,meanPeak,'color',cols(d),'linewidth',2)
        hold on
        %         patch([1:6 6:-1:1],[meanPeak-st fliplr(meanPeak+st)],cols(d),'facealpha',0.4)
        for q=1:6
            semilogx([spatials(q)+spatials(q)*0.1 spatials(q)+spatials(q)*0.1],[meanPeak(q)-st(q) meanPeak(q)+st(q)],'color',cols(d))
            semilogx([spatials(q)+spatials(q)*0.1 spatials(q)+spatials(q)*0.1],[meanPeak(q)-ci(q) meanPeak(q)+ci(q)],'color',cols(d),'linewidth',5)
            
        end
        
        %         set(gca,'xtick',1:6)
        %         set(gca,'xticklabel',spatials)
        xlabel('spatial period [um]')
        axis([50 5000 -inf inf])
        title([int2str(freqs(s)),'Hz / f',num2str(harmlist(2))])
        
        meanPeakS(d,scnt,:)=meanPeak;
        stS(d,scnt,:)=st;
        ciS(d,scnt,:)=ci;
        
    end
    
end



for d=1:3
    figure(100*d)
    for fcnt=1:6
        plot((freqs)+(fcnt-7)*0.1,reshape(meanPeakT(d,fcnt,:),4,1),'color',col2{d}(fcnt),'linewidth',2,'linestyle',style2{fcnt})
        hold on
    end
    legend(num2str(spatials'))
    for fcnt=1:6
        for q=1:4
            plot([freqs(q)+(fcnt-7)*0.1 freqs(q)+(fcnt-7)*0.1],[reshape(meanPeakT(d,fcnt,q),1,1)-reshape(stT(d,fcnt,q),1,1) reshape(meanPeakT(d,fcnt,q),1,1)+reshape(stT(d,fcnt,q),1,1)],'color',col2{d}(fcnt),'linestyle',style2{fcnt})
            plot([freqs(q)+(fcnt-7)*0.1 freqs(q)+(fcnt-7)*0.1],[reshape(meanPeakT(d,fcnt,q),1,1)-reshape(ciT(d,fcnt,q),1,1) reshape(meanPeakT(d,fcnt,q),1,1)+reshape(ciT(d,fcnt,q),1,1)],'color',col2{d}(fcnt),'linewidth',3,'linestyle',style2{fcnt})
            
        end
    end
    
    %     set(gca,'xtick',1:4)
    %     set(gca,'xticklabel',freqs)
    xlabel('frequency [Hz]')
    axis([0 9 -inf inf])
    title(['temporal tuning (',Dtype(d),')'])
    
    
    
    figure(1000*d)
    for fcnt=1:4
        semilogx((spatials)+spatials*0.05*(fcnt-3),reshape(meanPeakS(d,fcnt,:),6,1),'color',col1{d}(fcnt),'linewidth',2,'linestyle',style1{fcnt})
        hold on
    end
    legend(num2str(freqs'))
    for fcnt=1:4
        for q=1:6
            semilogx([spatials(q)+spatials(q)*0.05*(fcnt-3) spatials(q)+spatials(q)*0.05*(fcnt-3)],[reshape(meanPeakS(d,fcnt,q),1,1)-reshape(stS(d,fcnt,q),1,1) reshape(meanPeakS(d,fcnt,q),1,1)+reshape(stS(d,fcnt,q),1,1)],'color',col1{d}(fcnt),'linestyle',style1{fcnt})
            semilogx([spatials(q)+spatials(q)*0.05*(fcnt-3) spatials(q)+spatials(q)*0.05*(fcnt-3)],[reshape(meanPeakS(d,fcnt,q),1,1)-reshape(ciS(d,fcnt,q),1,1) reshape(meanPeakS(d,fcnt,q),1,1)+reshape(ciS(d,fcnt,q),1,1)],'color',col1{d}(fcnt),'linewidth',3,'linestyle',style1{fcnt})
        end
        
    end
    %     set(gca,'xtick',1:6)
    %     set(gca,'xticklabel',spatials)
    xlabel('spatial period [um]')
    axis([50 5000 -inf inf])
    title(['spatial tuning (',Dtype(d),')'])
end


for d=1:3
    tdata=collNorm{d};
    
    if version==2
        goods=reshape(goodT(d,:,:),24,1);
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            curr=curr(goods{m});
            meandata=[meandata mean(curr)];
        end
    else
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            meandata=[meandata mean(curr)];
        end
    end
    normthis=meandata/max(meandata);
    zi=reshape(normthis,length(unique(freqs)),length(unique(spatials)));
    zi=zi';
    zi=flipud(zi);
    figure(10000*d)
    set(gcf,'position',[1 31 1600 794])
    imagesc(1:length(unique(freqs)),1:length(unique(spatials)),zi)
    colormap(gray)
    set(gca,'YTick',1:length(unique(spatials)))
    set(gca,'XTick',1:length(unique(freqs)))
    set(gca,'XTickLabel',unique(freqs))
    xlabel('Hz')
    set(gca,'YTickLabel',fliplr(unique(spatials)))
    ylabel('period (um)')
end


statsS=zeros(3,4,6);%m-p,m-h, p-h
statsT=zeros(3,6,4); %6 spatial tunings, each 4 temporal freqs
cond=[1 2;1 3;2 3];
scnt=0;
for s=1:4:24 %go through temporal tuning
    scnt=scnt+1;
    sscnt=0;
    for ss=s:s+3
        sscnt=sscnt+1;
        for c=1:3
            if version==2
                data1=collNorm{cond(c,1)}(goodT{cond(c,1),scnt,sscnt},ss);
                data2=collNorm{cond(c,2)}(goodT{cond(c,2),scnt,sscnt},ss);
            else
                data1=collNorm{cond(c,1)}(:,ss);
                data2=collNorm{cond(c,2)}(:,ss);
            end
            p=ranksum(data1,data2);
            
            statsT(c,scnt,sscnt)=p; %comparison, figure,data point
        end
    end
end
scnt=0;
for s=1:4 %go through temporal tuning
    scnt=scnt+1;
    sscnt=0;
    for ss=s:4:24
        sscnt=sscnt+1;
        for c=1:3
            if version==2
                data1=collNorm{cond(c,1)}(goodS{cond(c,1),scnt,sscnt},ss);
                data2=collNorm{cond(c,2)}(goodS{cond(c,2),scnt,sscnt},ss);
            else
                data1=collNorm{cond(c,1)}(:,ss);
                data2=collNorm{cond(c,2)}(:,ss);
            end
            p=ranksum(data1,data2);
            
            statsS(c,scnt,sscnt)=p; %comparison, figure,data point
        end
    end
end



condlist={'w-p','w-h','p-h'};

fcnt=0;
for s=1:4:24
    
    fcnt=fcnt+1;
    figure(fcnt)
    yy=get(gca,'ytick');
    y=yy(end);
    icnt=0;
    for i=s:s+3
        icnt=icnt+1;
        for c=1:3
            curr=statsT(c,fcnt,icnt);
            curr=round(curr*1000)/1000;
            text(freqs(icnt),y+(y/8)*c,num2str(curr))
            text(freqs(icnt),y+(y/8)*(c+7),num2str(NofCT(c,fcnt,icnt)))
            if i==s
                text(0.5,y+(y/8)*c,condlist{c})
                text(0.5,y+(y/8)*(c+7),Dtype(c))
            end
        end
    end
    
    axis([0 9 yy(1) y+(y/8)*11])
end


% fcnt=0;
ssspatials=[100 300 600 1000 2000 4000];
for s=1:4
    
    fcnt=fcnt+1;
    figure(fcnt)
    yy=get(gca,'ytick');
    y=yy(end);
    icnt=0;
    for i=s:4:24
        icnt=icnt+1;
        for c=1:3
            curr=statsS(c,fcnt-6,icnt);
            curr=round(curr*1000)/1000;
            text(ssspatials(icnt),y+(y/8)*c,num2str(curr))
            text(ssspatials(icnt),y+(y/8)*(c+7),num2str(NofCS(c,fcnt-6,icnt)))
            if i==s
                text(10,y+(y/8)*c,condlist{c})
                text(10,y+(y/8)*(c+7),Dtype(c))
            end
        end
    end
    axis([50 5000 yy(1) y+(y/8)*11])
end





for f=1:6
    figure(f)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.emf']))
    
end
for f=1:4
    figure(f+6)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.emf']))
    
end

for d=1:3
    figure(100*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.emf']))
end
for d=1:3
    figure(10000*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_heatmap_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt2\all_stat',['stat_F1_heatmap_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
end
%% NICE PLOTS: DG absolute values
% plot for each temporal/spatial frequency the spatial/temporal tuning
% curve based on normalized F1 peaks (normalized for each cell to its max.
% F1)
% blue=mouse, green=pig, red=human
% including STD (thin line), and 95%CI (thick line)
% version1: take always all cells which have been classified as
% "DG-positive" (at least 1 response to one DG combination)
% version2: take for each single compared value only the cells which
% respond
clear
close all
superpath='f:\all_species_temp';
version=2;

col1={'bbcc','ggkk','rrmm'};
col2={'bbbccc','gggkkk','rrrmmm'};
style1={'-','--','-','--'};
style2={'-','--','-.','-','--','-.'};
cols='bgr';
Dtype='wph';
spatials=[100 200 500 1000 2000 4000];
freqs=[1 2 4 8];
harmlist=[0.5 1 2 3];

collNorm=cell(3,1);
NofCT=zeros(3,6,4);
NofCS=zeros(3,4,6);
goodT=cell(3,6,4);
goodS=cell(3,4,6);
% NofC=zeros(3,24);
meanPeakT=zeros(3,6,4);
meanPeakS=zeros(3,4,6);
stT=zeros(3,6,4);
stS=zeros(3,4,6);
ciT=zeros(3,6,4);
ciS=zeros(3,4,6);
for d=1:3
    fcnt=0;
    load(fullfile(superpath,'newAttempt2','DG',[Dtype(d),'_FFTpeaks']))
    collNorm{d,1}=reshape(FFTpeaks(:,:,2),size(FFTpeaks,1),24);
    scnt=0;
    for s=1:4:24
        scnt=scnt+1;
        
        fcnt=fcnt+1;
        figure(fcnt)
        
        if version==2
            meanPeak=[]; st=[]; ci=[]; gcnt=0;
            for g=s:s+3
                gcnt=gcnt+1;
                good=find(FFTpeaks(:,g,2)>0);
                goodT{d,scnt,gcnt}=good;
                
                meanPeak=[meanPeak mean(FFTpeaks(good,g,2))];
                st=[st std(FFTpeaks(good,g,2))];
                ci=[ci 1.96*st(end)/sqrt(length(good))];
                NofCT(d,scnt,gcnt)=length(good);
            end
        else
            meanPeak=mean(FFTpeaks(:,s:s+3,2));
            st=std(FFTpeaks(:,s:s+3,2));
            ci=1.96*st/sqrt(size(FFTpeaks,1));
            for g=1:4
                NofCT(d,scnt,g)=size(FFTpeaks,1);
                goodT{d,scnt,g}=1:size(FFTpeaks,1);
            end
        end
        
        %         st=std(normPeaks(:,s:s+3,2));
        %         ci=1.96*st/sqrt(size(normPeaks,1));
        %             normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
        plot((freqs)+(d-1)*0.1,meanPeak,'color',cols(d),'linewidth',2)
        hold on
        for q=1:4
            plot([freqs(q)+(d-1)*0.1 freqs(q)+(d-1)*0.1],[meanPeak(q)-st(q) meanPeak(q)+st(q)],'color',cols(d))
            plot([freqs(q)+(d-1)*0.1 freqs(q)+(d-1)*0.1],[meanPeak(q)-ci(q) meanPeak(q)+ci(q)],'color',cols(d),'linewidth',5)
        end
        %         set(gca,'xtick',1:4)
        %         set(gca,'xticklabel',freqs)
        xlabel('frequency [Hz]')
        axis([0 9 -inf inf])
        title([int2str(spatials(scnt)),'um / f',num2str(harmlist(2))])
        
        
        meanPeakT(d,scnt,:)=meanPeak;
        stT(d,scnt,:)=st;
        ciT(d,scnt,:)=ci;
        
    end
    scnt=0;
    for s=1:4
        fcnt=fcnt+1;
        figure(fcnt)
        scnt=scnt+1;
        
        if version==2
            meanPeak=[]; st=[]; ci=[]; gcnt=0;
            for g=s:4:24
                gcnt=gcnt+1;
                good=find(FFTpeaks(:,g,2)>0);
                goodS{d,scnt,gcnt}=good;
                
                meanPeak=[meanPeak mean(FFTpeaks(good,g,2))];
                st=[st std(FFTpeaks(good,g,2))];
                ci=[ci 1.96*st(end)/sqrt(length(good))];
                NofCS(d,scnt,gcnt)=length(good);
            end
        else
            meanPeak=mean(FFTpeaks(:,s:4:24,2));
            st=std(FFTpeaks(:,s:4:24,2));
            ci=1.96*st/sqrt(size(FFTpeaks,1));
            for g=1:6
                NofCS(d,scnt,g)=size(FFTpeaks,1);
                goodS{d,scnt,g}=1:size(FFTpeaks,1);
            end
        end
        
        %         meanPeak=mean(normPeaks(:,s:4:24,2));%contains peaks normalized to each cell's f1
        %         st=std(normPeaks(:,s:4:24,2));
        %                 ci=1.96*st/sqrt(size(normPeaks,1));
        
        %             normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
        plot((spatials)+(d-1)*20,meanPeak,'color',cols(d),'linewidth',2)
        hold on
        %         patch([1:6 6:-1:1],[meanPeak-st fliplr(meanPeak+st)],cols(d),'facealpha',0.4)
        for q=1:6
            plot([spatials(q)+(d-1)*20 spatials(q)+(d-1)*20],[meanPeak(q)-st(q) meanPeak(q)+st(q)],'color',cols(d))
            plot([spatials(q)+(d-1)*20 spatials(q)+(d-1)*20],[meanPeak(q)-ci(q) meanPeak(q)+ci(q)],'color',cols(d),'linewidth',5)
            
        end
        
        %         set(gca,'xtick',1:6)
        %         set(gca,'xticklabel',spatials)
        xlabel('spatial period [um]')
        axis([0 5000 -inf inf])
        title([int2str(freqs(s)),'Hz / f',num2str(harmlist(2))])
        
        meanPeakS(d,scnt,:)=meanPeak;
        stS(d,scnt,:)=st;
        ciS(d,scnt,:)=ci;
        
    end
    
end



for d=1:3
    figure(100*d)
    for fcnt=1:6
        plot((freqs)+(fcnt-7)*0.05,reshape(meanPeakT(d,fcnt,:),4,1),'color',col2{d}(fcnt),'linewidth',2,'linestyle',style2{fcnt})
        hold on
    end
    legend(num2str(spatials'))
    for fcnt=1:6
        for q=1:4
            plot([freqs(q)+(fcnt-7)*0.05 freqs(q)+(fcnt-7)*0.05],[reshape(meanPeakT(d,fcnt,q),1,1)-reshape(stT(d,fcnt,q),1,1) reshape(meanPeakT(d,fcnt,q),1,1)+reshape(stT(d,fcnt,q),1,1)],'color',col2{d}(fcnt),'linestyle',style2{fcnt})
            plot([freqs(q)+(fcnt-7)*0.05 freqs(q)+(fcnt-7)*0.05],[reshape(meanPeakT(d,fcnt,q),1,1)-reshape(ciT(d,fcnt,q),1,1) reshape(meanPeakT(d,fcnt,q),1,1)+reshape(ciT(d,fcnt,q),1,1)],'color',col2{d}(fcnt),'linewidth',3,'linestyle',style2{fcnt})
            
        end
    end
    
    %     set(gca,'xtick',1:4)
    %     set(gca,'xticklabel',freqs)
    xlabel('frequency [Hz]')
    axis([0 9 -inf inf])
    title(['temporal tuning (',Dtype(d),')'])
    
    
    
    figure(1000*d)
    for fcnt=1:4
        plot((spatials)+(fcnt-1)*20,reshape(meanPeakS(d,fcnt,:),6,1),'color',col1{d}(fcnt),'linewidth',2,'linestyle',style1{fcnt})
        hold on
    end
    legend(num2str(freqs'))
    for fcnt=1:4
        for q=1:6
            plot([spatials(q)+(fcnt-1)*20 spatials(q)+(fcnt-1)*20],[reshape(meanPeakS(d,fcnt,q),1,1)-reshape(stS(d,fcnt,q),1,1) reshape(meanPeakS(d,fcnt,q),1,1)+reshape(stS(d,fcnt,q),1,1)],'color',col1{d}(fcnt),'linestyle',style1{fcnt})
            plot([spatials(q)+(fcnt-1)*20 spatials(q)+(fcnt-1)*20],[reshape(meanPeakS(d,fcnt,q),1,1)-reshape(ciS(d,fcnt,q),1,1) reshape(meanPeakS(d,fcnt,q),1,1)+reshape(ciS(d,fcnt,q),1,1)],'color',col1{d}(fcnt),'linewidth',3,'linestyle',style1{fcnt})
        end
        
    end
    %      set(gca,'xtick',1:6)
    %         set(gca,'xticklabel',spatials)
    xlabel('spatial period [um]')
    axis([0 5000 -inf inf])
    title(['spatial tuning (',Dtype(d),')'])
end


for d=1:3
    tdata=collNorm{d};
    
    if version==2
        goods=reshape(goodT(d,:,:),24,1);
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            curr=curr(goods{m});
            meandata=[meandata mean(curr)];
        end
    else
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            meandata=[meandata mean(curr)];
        end
    end
    normthis=meandata/max(meandata);
    zi=reshape(normthis,length(unique(freqs)),length(unique(spatials)));
    zi=zi';
    zi=flipud(zi);
    figure(10000*d)
    set(gcf,'position',[1 31 1600 794])
    imagesc(1:length(unique(freqs)),1:length(unique(spatials)),zi)
    colormap(gray)
    set(gca,'YTick',1:length(unique(spatials)))
    set(gca,'XTick',1:length(unique(freqs)))
    set(gca,'XTickLabel',unique(freqs))
    xlabel('Hz')
    set(gca,'YTickLabel',fliplr(unique(spatials)))
    ylabel('period (um)')
end


statsS=zeros(3,4,6);%m-p,m-h, p-h
statsT=zeros(3,6,4); %6 spatial tunings, each 4 temporal freqs
cond=[1 2;1 3;2 3];
scnt=0;
for s=1:4:24 %go through temporal tuning
    scnt=scnt+1;
    sscnt=0;
    for ss=s:s+3
        sscnt=sscnt+1;
        for c=1:3
            if version==2
                data1=collNorm{cond(c,1)}(goodT{cond(c,1),scnt,sscnt},ss);
                data2=collNorm{cond(c,2)}(goodT{cond(c,2),scnt,sscnt},ss);
            else
                data1=collNorm{cond(c,1)}(:,ss);
                data2=collNorm{cond(c,2)}(:,ss);
            end
            p=ranksum(data1,data2);
            
            statsT(c,scnt,sscnt)=p; %comparison, figure,data point
        end
    end
end
scnt=0;
for s=1:4 %go through temporal tuning
    scnt=scnt+1;
    sscnt=0;
    for ss=s:4:24
        sscnt=sscnt+1;
        for c=1:3
            if version==2
                data1=collNorm{cond(c,1)}(goodS{cond(c,1),scnt,sscnt},ss);
                data2=collNorm{cond(c,2)}(goodS{cond(c,2),scnt,sscnt},ss);
            else
                data1=collNorm{cond(c,1)}(:,ss);
                data2=collNorm{cond(c,2)}(:,ss);
            end
            p=ranksum(data1,data2);
            
            statsS(c,scnt,sscnt)=p; %comparison, figure,data point
        end
    end
end



condlist={'w-p','w-h','p-h'};

fcnt=0;
for s=1:4:24
    
    fcnt=fcnt+1;
    figure(fcnt)
    yy=get(gca,'ytick');
    y=yy(end);
    icnt=0;
    for i=s:s+3
        icnt=icnt+1;
        for c=1:3
            curr=statsT(c,fcnt,icnt);
            curr=round(curr*1000)/1000;
            text(freqs(icnt),y+(y/8)*c,num2str(curr))
            text(freqs(icnt),y+(y/8)*(c+7),num2str(NofCT(c,fcnt,icnt)))
            if i==s
                text(0.5,y+(y/8)*c,condlist{c})
                text(0.5,y+(y/8)*(c+7),Dtype(c))
            end
        end
    end
    
    axis([0 9 yy(1) y+(y/8)*11])
end


ssspatials=[100 300 600 1000 2000 4000];
for s=1:4
    
    fcnt=fcnt+1;
    figure(fcnt)
    yy=get(gca,'ytick');
    y=yy(end);
    icnt=0;
    for i=s:4:24
        icnt=icnt+1;
        for c=1:3
            curr=statsS(c,fcnt-6,icnt);
            curr=round(curr*1000)/1000;
            text(ssspatials(icnt),y+(y/8)*c,num2str(curr))
            text(ssspatials(icnt),y+(y/8)*(c+7),num2str(NofCS(c,fcnt-6,icnt)))
            if i==s
                text(10,y+(y/8)*c,condlist{c})
                text(10,y+(y/8)*(c+7),Dtype(c))
            end
        end
    end
    axis([0 5000 yy(1) y+(y/8)*11])
end





for f=1:6
    figure(f)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_',int2str(spatials(f)),'um_absPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_',int2str(spatials(f)),'um_absPeaks_version',int2str(version),'.bmp']))
end
for f=1:4
    figure(f+6)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_',int2str(freqs(f)),'Hz_absPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_',int2str(freqs(f)),'Hz_absPeaks_version',int2str(version),'.bmp']))
end

for d=1:3
    figure(100*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_temporal_',Dtype(d),'_absPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_temporal_',Dtype(d),'_absPeaks_version',int2str(version),'.bmp']))
    figure(1000*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_spatial_',Dtype(d),'_absPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_spatial_',Dtype(d),'_absPeaks_version',int2str(version),'.bmp']))
end
for d=1:3
    figure(10000*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_heatmap_',Dtype(d),'_absPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['astat_F1_heatmap_',Dtype(d),'_absPeaks_version',int2str(version),'.bmp']))
end
%% NICE PLOTS: DG normPeaks: sum F1toF3
% plot for each temporal/spatial frequency the spatial/temporal tuning
% curve based on normalized F1 peaks (normalized for each cell to its max.
% F1)
% blue=mouse, green=pig, red=human
% including STD (thin line), and 95%CI (thick line)
% version1: take always all cells which have been classified as
% "DG-positive" (at least 1 response to one DG combination)
% version2: take for each single compared value only the cells which
% respond


% => version1: is better
clear
close all
superpath='f:\all_species_temp';
version=1;

col1={'bbcc','ggkk','rrmm'};
col2={'bbbccc','gggkkk','rrrmmm'};
style1={'-','--','-','--'};
style2={'-','--','-.','-','--','-.'};
cols='bgr';
Dtype='wph';
spatials=[100 200 500 1000 2000 4000];
freqs=[1 2 4 8];
harmlist=[0.5 1 2 3];

collNorm=cell(3,1);
NofCT=zeros(3,6,4);
NofCS=zeros(3,4,6);
goodT=cell(3,6,4);
goodS=cell(3,4,6);
% NofC=zeros(3,24);
meanPeakT=zeros(3,6,4);
meanPeakS=zeros(3,4,6);
stT=zeros(3,6,4);
stS=zeros(3,4,6);
ciT=zeros(3,6,4);
ciS=zeros(3,4,6);
for d=1:3
    fcnt=0;
    load(fullfile(superpath,'newAttempt','DG',[Dtype(d),'_normPeaks']))
    collNorm{d,1}=reshape(sum(normPeaks(:,:,2:4),3),size(normPeaks,1),24);
    scnt=0;
    for s=1:4:24
        scnt=scnt+1;
        
        fcnt=fcnt+1;
        figure(fcnt)
        
        if version==2
            meanPeak=[]; st=[]; ci=[]; gcnt=0;
            for g=s:s+3
                gcnt=gcnt+1;
                good=find(normPeaks(:,g,2)>0);
                goodT{d,scnt,gcnt}=good;
                
                meanPeak=[meanPeak mean(sum(normPeaks(good,g,2:4),3))];
                st=[st std(sum(normPeaks(good,g,2:4),3))];
                ci=[ci 1.96*st(end)/sqrt(length(good))];
                NofCT(d,scnt,gcnt)=length(good);
            end
        else
            meanPeak=mean(sum(normPeaks(:,s:s+3,2:4),3));
            st=std(sum(normPeaks(:,s:s+3,2:4),3));
            ci=1.96*st/sqrt(size(normPeaks,1));
            for g=1:4
                NofCT(d,scnt,g)=size(normPeaks,1);
                goodT{d,scnt,g}=1:size(normPeaks,1);
            end
        end
        
        %         st=std(normPeaks(:,s:s+3,2));
        %         ci=1.96*st/sqrt(size(normPeaks,1));
        %             normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
        plot((freqs)+(d-1)*0.1,meanPeak,'color',cols(d),'linewidth',2)
        hold on
        for q=1:4
            plot([freqs(q)+(d-1)*0.1 freqs(q)+(d-1)*0.1],[meanPeak(q)-st(q) meanPeak(q)+st(q)],'color',cols(d))
            plot([freqs(q)+(d-1)*0.1 freqs(q)+(d-1)*0.1],[meanPeak(q)-ci(q) meanPeak(q)+ci(q)],'color',cols(d),'linewidth',5)
        end
        %         set(gca,'xtick',1:4)
        %         set(gca,'xticklabel',freqs)
        xlabel('frequency [Hz]')
        axis([0 9 -inf inf])
        title([int2str(spatials(scnt)),'um / f',num2str(harmlist(2))])
        
        
        meanPeakT(d,scnt,:)=meanPeak;
        stT(d,scnt,:)=st;
        ciT(d,scnt,:)=ci;
        
    end
    scnt=0;
    for s=1:4
        fcnt=fcnt+1;
        figure(fcnt)
        scnt=scnt+1;
        
        if version==2
            meanPeak=[]; st=[]; ci=[]; gcnt=0;
            for g=s:4:24
                gcnt=gcnt+1;
                good=find(normPeaks(:,g,2)>0);
                goodS{d,scnt,gcnt}=good;
                
                meanPeak=[meanPeak mean(sum(normPeaks(good,g,2:4),3))];
                st=[st std(sum(normPeaks(good,g,2:4),3))];
                ci=[ci 1.96*st(end)/sqrt(length(good))];
                NofCS(d,scnt,gcnt)=length(good);
            end
        else
            meanPeak=mean(sum(normPeaks(:,s:4:24,2:4),3));
            st=std(sum(normPeaks(:,s:4:24,2:4),3));
            ci=1.96*st/sqrt(size(normPeaks,1));
            for g=1:6
                NofCS(d,scnt,g)=size(normPeaks,1);
                goodS{d,scnt,g}=1:size(normPeaks,1);
            end
        end
        
        %         meanPeak=mean(normPeaks(:,s:4:24,2));%contains peaks normalized to each cell's f1
        %         st=std(normPeaks(:,s:4:24,2));
        %                 ci=1.96*st/sqrt(size(normPeaks,1));
        
        %             normMean=(meanPeak-min(meanPeak))./(max(meanPeak)-min(meanPeak));
        semilogx(([100 200 500 1000 2000 4000])+spatials*0.1,meanPeak,'color',cols(d),'linewidth',2)
        hold on
        %         patch([1:6 6:-1:1],[meanPeak-st fliplr(meanPeak+st)],cols(d),'facealpha',0.4)
        for q=1:6
            semilogx([spatials(q)+spatials(q)*0.1 spatials(q)+spatials(q)*0.1],[meanPeak(q)-st(q) meanPeak(q)+st(q)],'color',cols(d))
            semilogx([spatials(q)+spatials(q)*0.1 spatials(q)+spatials(q)*0.1],[meanPeak(q)-ci(q) meanPeak(q)+ci(q)],'color',cols(d),'linewidth',5)
            
        end
        
        %         set(gca,'xtick',1:6)
        %         set(gca,'xticklabel',spatials)
        xlabel('spatial period [um]')
        axis([50 5000 -inf inf])
        title([int2str(freqs(s)),'Hz / f',num2str(harmlist(2))])
        
        meanPeakS(d,scnt,:)=meanPeak;
        stS(d,scnt,:)=st;
        ciS(d,scnt,:)=ci;
        
    end
    
end



for d=1:3
    figure(100*d)
    for fcnt=1:6
        plot((freqs)+(fcnt-7)*0.1,reshape(meanPeakT(d,fcnt,:),4,1),'color',col2{d}(fcnt),'linewidth',2,'linestyle',style2{fcnt})
        hold on
    end
    legend(num2str(spatials'))
    for fcnt=1:6
        for q=1:4
            plot([freqs(q)+(fcnt-7)*0.1 freqs(q)+(fcnt-7)*0.1],[reshape(meanPeakT(d,fcnt,q),1,1)-reshape(stT(d,fcnt,q),1,1) reshape(meanPeakT(d,fcnt,q),1,1)+reshape(stT(d,fcnt,q),1,1)],'color',col2{d}(fcnt),'linestyle',style2{fcnt})
            plot([freqs(q)+(fcnt-7)*0.1 freqs(q)+(fcnt-7)*0.1],[reshape(meanPeakT(d,fcnt,q),1,1)-reshape(ciT(d,fcnt,q),1,1) reshape(meanPeakT(d,fcnt,q),1,1)+reshape(ciT(d,fcnt,q),1,1)],'color',col2{d}(fcnt),'linewidth',3,'linestyle',style2{fcnt})
            
        end
    end
    
    %     set(gca,'xtick',1:4)
    %     set(gca,'xticklabel',freqs)
    xlabel('frequency [Hz]')
    axis([0 9 -inf inf])
    title(['temporal tuning (',Dtype(d),')'])
    
    
    
    figure(1000*d)
    for fcnt=1:4
        semilogx((spatials)+spatials*0.05*(fcnt-3),reshape(meanPeakS(d,fcnt,:),6,1),'color',col1{d}(fcnt),'linewidth',2,'linestyle',style1{fcnt})
        hold on
    end
    legend(num2str(freqs'))
    for fcnt=1:4
        for q=1:6
            semilogx([spatials(q)+spatials(q)*0.05*(fcnt-3) spatials(q)+spatials(q)*0.05*(fcnt-3)],[reshape(meanPeakS(d,fcnt,q),1,1)-reshape(stS(d,fcnt,q),1,1) reshape(meanPeakS(d,fcnt,q),1,1)+reshape(stS(d,fcnt,q),1,1)],'color',col1{d}(fcnt),'linestyle',style1{fcnt})
            semilogx([spatials(q)+spatials(q)*0.05*(fcnt-3) spatials(q)+spatials(q)*0.05*(fcnt-3)],[reshape(meanPeakS(d,fcnt,q),1,1)-reshape(ciS(d,fcnt,q),1,1) reshape(meanPeakS(d,fcnt,q),1,1)+reshape(ciS(d,fcnt,q),1,1)],'color',col1{d}(fcnt),'linewidth',3,'linestyle',style1{fcnt})
        end
        
    end
    %     set(gca,'xtick',1:6)
    %     set(gca,'xticklabel',spatials)
    xlabel('spatial period [um]')
    axis([50 5000 -inf inf])
    title(['spatial tuning (',Dtype(d),')'])
end


for d=1:3
    tdata=collNorm{d};
    
    if version==2
        goods=reshape(goodT(d,:,:),24,1);
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            curr=curr(goods{m});
            meandata=[meandata mean(curr)];
        end
    else
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            meandata=[meandata mean(curr)];
        end
    end
    normthis=meandata/max(meandata);
    zi=reshape(normthis,length(unique(freqs)),length(unique(spatials)));
    zi=zi';
    zi=flipud(zi);
    figure(10000*d)
    set(gcf,'position',[1 31 1600 794])
    imagesc(1:length(unique(freqs)),1:length(unique(spatials)),zi)
    colormap(gray)
    set(gca,'YTick',1:length(unique(spatials)))
    set(gca,'XTick',1:length(unique(freqs)))
    set(gca,'XTickLabel',unique(freqs))
    xlabel('Hz')
    set(gca,'YTickLabel',fliplr(unique(spatials)))
    ylabel('period (um)')
end


statsS=zeros(3,4,6);%m-p,m-h, p-h
statsT=zeros(3,6,4); %6 spatial tunings, each 4 temporal freqs
cond=[1 2;1 3;2 3];
scnt=0;
for s=1:4:24 %go through temporal tuning
    scnt=scnt+1;
    sscnt=0;
    for ss=s:s+3
        sscnt=sscnt+1;
        for c=1:3
            if version==2
                data1=collNorm{cond(c,1)}(goodT{cond(c,1),scnt,sscnt},ss);
                data2=collNorm{cond(c,2)}(goodT{cond(c,2),scnt,sscnt},ss);
            else
                data1=collNorm{cond(c,1)}(:,ss);
                data2=collNorm{cond(c,2)}(:,ss);
            end
            p=ranksum(data1,data2);
            
            statsT(c,scnt,sscnt)=p; %comparison, figure,data point
        end
    end
end
scnt=0;
for s=1:4 %go through temporal tuning
    scnt=scnt+1;
    sscnt=0;
    for ss=s:4:24
        sscnt=sscnt+1;
        for c=1:3
            if version==2
                data1=collNorm{cond(c,1)}(goodS{cond(c,1),scnt,sscnt},ss);
                data2=collNorm{cond(c,2)}(goodS{cond(c,2),scnt,sscnt},ss);
            else
                data1=collNorm{cond(c,1)}(:,ss);
                data2=collNorm{cond(c,2)}(:,ss);
            end
            p=ranksum(data1,data2);
            
            statsS(c,scnt,sscnt)=p; %comparison, figure,data point
        end
    end
end



condlist={'w-p','w-h','p-h'};

fcnt=0;
for s=1:4:24
    
    fcnt=fcnt+1;
    figure(fcnt)
    yy=get(gca,'ytick');
    y=yy(end);
    icnt=0;
    for i=s:s+3
        icnt=icnt+1;
        for c=1:3
            curr=statsT(c,fcnt,icnt);
            curr=round(curr*1000)/1000;
            text(freqs(icnt),y+(y/8)*c,num2str(curr))
            text(freqs(icnt),y+(y/8)*(c+7),num2str(NofCT(c,fcnt,icnt)))
            if i==s
                text(0.5,y+(y/8)*c,condlist{c})
                text(0.5,y+(y/8)*(c+7),Dtype(c))
            end
        end
    end
    
    axis([0 9 yy(1) y+(y/8)*11])
end


% fcnt=0;
ssspatials=[100 300 600 1000 2000 4000];
for s=1:4
    
    fcnt=fcnt+1;
    figure(fcnt)
    yy=get(gca,'ytick');
    y=yy(end);
    icnt=0;
    for i=s:4:24
        icnt=icnt+1;
        for c=1:3
            curr=statsS(c,fcnt-6,icnt);
            curr=round(curr*1000)/1000;
            text(ssspatials(icnt),y+(y/8)*c,num2str(curr))
            text(ssspatials(icnt),y+(y/8)*(c+7),num2str(NofCS(c,fcnt-6,icnt)))
            if i==s
                text(10,y+(y/8)*c,condlist{c})
                text(10,y+(y/8)*(c+7),Dtype(c))
            end
        end
    end
    axis([50 5000 yy(1) y+(y/8)*11])
end





for f=1:6
    figure(f)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.emf']))
    
end
for f=1:4
    figure(f+6)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.emf']))
    
end

for d=1:3
    figure(100*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.emf']))
    
    figure(1000*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.emf']))
end
for d=1:3
    figure(10000*d)
    set(gcf,'position',[1 31 1600 794])
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_heatmap_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
    saveas(gcf,fullfile('f:\all_species_temp\newAttempt\all_stat',['stat_allF_heatmap_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
end
%% look at human RGCs that respond to 100/200 um

superpath='E:\all_species_temp';

load(fullfile(superpath,'newAttempt','DG','h_normPeaks'))

interesting=[];
for i=1:length(normPeaks)
    if ~isempty(find( normPeaks(i,1:8,2)>0))
        interesting=[interesting;i];
    end
end

data=reshape(normPeaks(:,1:8,2),size(normPeaks,1),8);

datai=data(interesting,:);

ids=togo(interesting);
%% TABLE: which cell responds to what stimulus
% 1) flash, 2) LF, 3) DG, 4) chirp, 5) 6vel black, 6) 6vel white, 7) DS

clear
close all
superpath='f:\all_species_temp';

Dtype='wph';

% idlist=[8 9 0 19 16 17 10/11];

for d=1:3
    
    load(fullfile(superpath,[Dtype(d),'_info']))
    load(fullfile(superpath,'newAttempt','DG',[Dtype(d),'_normPeaks']))
    response=zeros(size(info,1),7);
    pols=zeros(size(info,1),2);
    pols=pols-2;
    dates=unique(info(:,1));
    did=zeros(size(info,1),1);
    
    
    for i=1:size(info,1)
        found=0; m=0;
        while found==0
            m=m+1;
            if ~isempty(regexp(dates{m},info{i,1}))
                found=1;
            end
        end
        did(i)=m;
        
        
        if ~isempty(info{i,8}) %flash
            if info{i,8}~=0 && ~isnan(info{i,8})
                response(i,1)=1;
                if info{i,8}==1;
                    pols(i,1)=1;
                elseif info{i,8}==-1
                    pols(i,1)=-1;
                elseif info{i,8}==2
                    pols(i,1)=0;
                end
            end
        end
        if ~isempty(info{i,9}) %LF
            if info{i,9}~=0 && ~isnan(info{i,9})
                response(i,2)=1;
                if info{i,9}==1;
                    pols(i,2)=1;
                elseif info{i,9}==-1
                    pols(i,2)=-1;
                end
            end
        end
        if ~isempty(info{i,19})
            if info{i,19}~=0 && ~isnan(info{i,19})
                response(i,4)=1;
            end
        end
        if ~isempty(info{i,16})
            if info{i,16}~=0 && ~isnan(info{i,16})
                response(i,5)=1;
            end
        end
        if ~isempty(info{i,17})
            if info{i,17}~=0 && ~isnan(info{i,17})
                response(i,6)=1;
            end
        end
        
        %         if ~isempty(info{i,10}) ||  ~isempty(info{i,11})
        %             if info{i,10}~=0 || info{i,11}~=0
        %                 if ~isnan(info{i,10}) || ~isnan(info{i,11})
        %                     response(i,7)=1;
        %                 end
        %             end
        %         end
        if ~isempty(info{i,12})
            if info{i,12}~=0 && ~isnan(info{i,12})
                response(i,7)=1;
            end
        end
        
        if d==2
            if ~isempty(find(did(i)==[1 2 3 7]));
                response(i,4)=-1;
            end
            if ~isempty(find(did(i)==7));
                response(i,3)=-1;
            end
            if ~isempty(find(did(i)==9));
                response(i,[3 5 6])=-1;
            end
        end
        
        if d==3
            if ~isempty(find(did(i)==[1 2 3]));
                response(i,4)=-1;
            end
        end
    end
    response(togo,3)=1;
    tmp=find(goodones==0);
    tst=find(response==-1);
    response(tmp,:)=-1;
    response(tst)=-1;
    
    
    
    di=diff(did);
    tmp=find(di~=0);
    tmp=[0; tmp];
    tmp=tmp+1;
    changingPoints=tmp;
    
    dates2=dates(did(changingPoints));
    
    %     if d==1
    %     tst=find(response==0);
    %     response(tst)=0.5;
    %     end
    
    
    % response: -1=not considered/not measured, 0=no response, 1=response
    cols=[0 0 0;0.7 0.7 0.7;1 0 0];
    figure
    for i=1:length(response)
        for j=1:7
            rectangle('position',[j,length(response)-i,1,1],'facecolor',cols(response(i,j)+2,:),'edgecolor',cols(response(i,j)+2,:))
            hold on
        end
    end
    changing=length(response)-changingPoints+1;
    set(gca,'xtick',1.5:7.5)
    set(gca,'xticklabel',{'flash','LF','DG','chirp','6velB','6velW','DS'})
    set(gca,'ytick',flipud(changing));
    set(gca,'yticklabel',flipud(dates2))
    set(gca, 'TickDir', 'out')
    set(gcf,'position',[1766 18 402 799])
    axis([1 8 0 length(response)])
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_responses.bmp']))
    saveas(gcf, fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_responses.fig']))
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_responses.emf']))
    
    
    %pols: 1=ON, 0=ONOFF, -1=OFF -2=none/not considered
    cols2=[0 0 0;0 0.3 1;0 1 0;1 0 0];
    figure
    for i=1:length(pols)
        for j=1:2
            rectangle('position',[j,length(response)-i,1,1],'facecolor',cols2(pols(i,j)+3,:),'edgecolor',cols2(pols(i,j)+3,:))
            hold on
        end
    end
    set(gca,'xtick',1.5:2.5)
    set(gca,'xticklabel',{'flash','LF'})
    set(gca,'ytick',flipud(changing));
    set(gca,'yticklabel',flipud(dates2))
    set(gca, 'TickDir', 'out')
    set(gcf,'position',[1766 18 153 799])
    axis([1 3 0 length(response)])
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_polarity.bmp']))
    saveas(gcf, fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_polarity.fig']))
    saveas(gcf,fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_polarity.emf']))
    
    %
    %     figure
    %     imagesc(response,[-1 1])
    %     colormap('copper')
    %     set(gca,'xtick',1:7)
    %     set(gca,'xticklabel',{'flash','LF','DG','chirp','6velB','6velW','DS'})
    %     set(gca,'ytick',changingPoints);
    %     set(gca,'yticklabel',dates2)
    %     set(gca, 'TickDir', 'out')
    %     set(gcf,'position',[1766 18 402 799])
    %     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_responses.bmp']))
    %     saveas(gcf, fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_responses.fig']))
    %     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_responses.emf']))
    %
    %      figure
    %     imagesc(pols,[-2 1])
    %     colormap('jet')
    %     set(gca,'xtick',1:2)
    %     set(gca,'xticklabel',{'flash','LF'})
    %     set(gca,'ytick',changingPoints);
    %     set(gca,'yticklabel',dates2)
    %     set(gca, 'TickDir', 'out')
    %     set(gcf,'position',[1766 18 153 799])
    %     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_polarity.bmp']))
    %     saveas(gcf, fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_polarity.fig']))
    %     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\stat_figs',[Dtype(d),'_polarity.emf']))
    names=info(:,2);
    save(fullfile('E:\all_species_temp',[Dtype(d),'_responses']),'names','response','pols','changingPoints','dates','dates2','did')
end

%% --- HERE STARTS CLUSTERING SPATIOTEMP ---
%% PC on LFs (from Field 2007)
% -> no real clusters (smoothing changes clustering quite a lot)

clear
% close all
superpath='E:\all_species_temp';

Dtype='wph';


for d=1:3
    
    load(fullfile(superpath,[Dtype(d),'_LFs']))
    load(fullfile(superpath,[Dtype(d),'_info']))
    go=[]; pol=[];
    for i=1:length(info)
        if goodones(i)==1 && info{i,9}~=0
            go=[go;i];
            pol=[pol;info{i,9}];
        end
    end
    tmp=find(pol==-1);
    pol(tmp)=0;
    data=LFs(:,go)';
    
    data2=[]; wind=5;
    for j=1:size(data,1)
        smoo=[];
        test=data(j,:);
        for t=1:500-wind
            smoo=[smoo mean(test(t:t+wind-1))];
        end
        data2=[data2;smoo];
    end
    
    x=cov(data2);
    [V,~]=eig(x);
    pc_vectors=V(:,end-2:end);
    pc1=data2*pc_vectors(:,end);
    pc2=data2*pc_vectors(:,end-1);
    pc3=data2*pc_vectors(:,end-2);
    
    cols=[0 0 1;1 0 0];
    %     for i=1:size(data,1)
    %         figure
    %         plot(data2(i,:),'color',cols(pol(i)+1,:))
    %         title(num2str(pc1(i)),'interpreter','none')
    %         saveas(gcf,fullfile('E:\all_species_temp\LF_test',[Dtype(d),'_',int2str(go(i)),'_LF.bmp']))
    %         close
    %     end
    %
    
    figure
    subplot(2,2,1)
    scatter(pc1,pc2,20,cols(pol+1,:))
    xlabel('pc1')
    ylabel('pc2')
    subplot(2,2,2)
    scatter(pc1,pc3,20,cols(pol+1,:))
    xlabel('pc1')
    ylabel('pc3')
    subplot(2,2,3)
    scatter(pc2,pc3,20,cols(pol+1,:))
    xlabel('pc2')
    ylabel('pc3')
    
end
%% !!!PC on DG Fourier transform
% => eventually I used newAttempt2: 5 PCs and center of mass
clear
close all

fr=2;
ids=[];
load(fullfile('F:\all_species_temp\newAttempt2\DG','w_normPeaks'))
data1=normPeaks(:,:,fr);
% data1=fliplr(reshape(data1,size(data1,1),4,6));
ids=[ids;repmat(1,length(data1),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','p_normPeaks'))
data2=normPeaks(:,:,fr);
% data2=fliplr(reshape(data2,size(data2,1),4,6));
ids=[ids;repmat(2,length(data2),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','h_normPeaks'))
data3=normPeaks(:,:,fr);
% data3=fliplr(reshape(data3,size(data3,1),4,6));
ids=[ids;repmat(3,length(data3),1)];

data=cat(1,data1,data2,data3);
x=cov(data);
[V,~]=eig(x);
pc_vectors=V(:,end-5:end);
pc1=data*pc_vectors(:,end);
pc2=data*pc_vectors(:,end-1);
pc3=data*pc_vectors(:,end-2);
pc4=data*pc_vectors(:,end-3);
pc5=data*pc_vectors(:,end-4);

cols=[0 0 1;0 1 0;1 0 0];
figure
subplot(2,3,1)
scatter(pc1,pc2,15,cols(ids,:))
subplot(2,3,2)
scatter(pc1,pc3,15,cols(ids,:))
subplot(2,3,3)
scatter(pc1,pc4,15,cols(ids,:))
subplot(2,3,4)
scatter(pc2,pc3,15,cols(ids,:))
subplot(2,3,5)
scatter(pc2,pc4,15,cols(ids,:))
subplot(2,3,6)
scatter(pc3,pc4,15,cols(ids,:))
pcs=[pc1 pc2 pc3 pc4 pc5];
tID=ids;

% add center of mass
cmass=[];
for d=1:length(data)
curr=fliplr(reshape(data(d,:),4,6));
curr=curr';
tmp=find(curr<0.2);
binaryImage = true(size(curr));
binaryImage(tmp)=0;
labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage, curr, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
cmass=[cmass;centerOfMass];
end
pcs=[pcs cmass];


save(fullfile('F:\all_species_temp\newAttempt2\DG_PCanalysis',['PC_F',int2str(fr-1)]),'pcs','tID')
% save(fullfile('F:\all_species_temp\newAttempt2\DG_PCanalysis',['PC_F',int2str(fr-1),'new']),'pcs','tID')

%% PC on DG Fourier transform: only human
clear
close all

fr=2;
ids=[];
load(fullfile('e:\all_species_temp\newAttempt2\DG','h_normPeaks'))
data=normPeaks(:,:,fr);
% data1=fliplr(reshape(data1,size(data1,1),4,6));
ids=[ids;repmat(3,length(data),1)];


x=cov(data);
[V,~]=eig(x);
pc_vectors=V(:,end-5:end);
pc1=data*pc_vectors(:,end);
pc2=data*pc_vectors(:,end-1);
pc3=data*pc_vectors(:,end-2);
pc4=data*pc_vectors(:,end-3);
pc5=data*pc_vectors(:,end-4);

cols=[0 0 1;0 1 0;1 0 0];
figure
subplot(2,3,1)
scatter(pc1,pc2,15,cols(ids,:))
subplot(2,3,2)
scatter(pc1,pc3,15,cols(ids,:))
subplot(2,3,3)
scatter(pc1,pc4,15,cols(ids,:))
subplot(2,3,4)
scatter(pc2,pc3,15,cols(ids,:))
subplot(2,3,5)
scatter(pc2,pc4,15,cols(ids,:))
subplot(2,3,6)
scatter(pc3,pc4,15,cols(ids,:))
pcs=[pc1 pc2 pc3 pc4 pc5];
tID=ids;

% add center of mass
cmass=[];
for d=1:length(data)
curr=fliplr(reshape(data(d,:),4,6));
curr=curr';
tmp=find(curr<0.2);
binaryImage = true(size(curr));
binaryImage(tmp)=0;
labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage, curr, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
cmass=[cmass;centerOfMass];
end
pcs=[pcs cmass];


save(fullfile('e:\all_species_temp\newAttemptOnlyHum\DG_PCanalysis',['PC_F',int2str(fr-1)]),'pcs','tID')
% save(fullfile('F:\all_species_temp\newAttempt2\DG_PCanalysis',['PC_F',int2str(fr-1),'new']),'pcs','tID')
%% PC on DG Fourier transform: only mouse
clear
close all

fr=2;
ids=[];
load(fullfile('e:\all_species_temp\newAttempt2\DG','w_normPeaks'))
data=normPeaks(:,:,fr);
% data1=fliplr(reshape(data1,size(data1,1),4,6));
ids=[ids;repmat(3,length(data),1)];


x=cov(data);
[V,~]=eig(x);
pc_vectors=V(:,end-5:end);
pc1=data*pc_vectors(:,end);
pc2=data*pc_vectors(:,end-1);
pc3=data*pc_vectors(:,end-2);
pc4=data*pc_vectors(:,end-3);
pc5=data*pc_vectors(:,end-4);

cols=[0 0 1;0 1 0;1 0 0];
figure
subplot(2,3,1)
scatter(pc1,pc2,15,cols(ids,:))
subplot(2,3,2)
scatter(pc1,pc3,15,cols(ids,:))
subplot(2,3,3)
scatter(pc1,pc4,15,cols(ids,:))
subplot(2,3,4)
scatter(pc2,pc3,15,cols(ids,:))
subplot(2,3,5)
scatter(pc2,pc4,15,cols(ids,:))
subplot(2,3,6)
scatter(pc3,pc4,15,cols(ids,:))
pcs=[pc1 pc2 pc3 pc4 pc5];
tID=ids;

% add center of mass
cmass=[];
for d=1:length(data)
curr=fliplr(reshape(data(d,:),4,6));
curr=curr';
tmp=find(curr<0.2);
binaryImage = true(size(curr));
binaryImage(tmp)=0;
labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage, curr, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
cmass=[cmass;centerOfMass];
end
pcs=[pcs cmass];


save(fullfile('e:\all_species_temp\newAttemptOnlyMouse\DG_PCanalysis',['PC_F',int2str(fr-1)]),'pcs','tID')
% save(fullfile('F:\all_species_temp\newAttempt2\DG_PCanalysis',['PC_F',int2str(fr-1),'new']),'pcs','tID')
%% PC on DG Fourier transform for newAttempt4: set <0.2 to 0
clear
close all

fr=2;
ids=[];
load(fullfile('F:\all_species_temp\newAttempt2\DG','w_normPeaks'))
data1=normPeaks(:,:,fr);
% data1=fliplr(reshape(data1,size(data1,1),4,6));
ids=[ids;repmat(1,length(data1),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','p_normPeaks'))
data2=normPeaks(:,:,fr);
% data2=fliplr(reshape(data2,size(data2,1),4,6));
ids=[ids;repmat(2,length(data2),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','h_normPeaks'))
data3=normPeaks(:,:,fr);
% data3=fliplr(reshape(data3,size(data3,1),4,6));
ids=[ids;repmat(3,length(data3),1)];

data=cat(1,data1,data2,data3);
tmp=find(data<0.2);
data(tmp)=0;
x=cov(data);
[V,~]=eig(x);
pc_vectors=V(:,end-5:end);
pc1=data*pc_vectors(:,end);
pc2=data*pc_vectors(:,end-1);
pc3=data*pc_vectors(:,end-2);
pc4=data*pc_vectors(:,end-3);
pc5=data*pc_vectors(:,end-4);

cols=[0 0 1;0 1 0;1 0 0];
figure
subplot(2,3,1)
scatter(pc1,pc2,15,cols(ids,:))
subplot(2,3,2)
scatter(pc1,pc3,15,cols(ids,:))
subplot(2,3,3)
scatter(pc1,pc4,15,cols(ids,:))
subplot(2,3,4)
scatter(pc2,pc3,15,cols(ids,:))
subplot(2,3,5)
scatter(pc2,pc4,15,cols(ids,:))
subplot(2,3,6)
scatter(pc3,pc4,15,cols(ids,:))
pcs=[pc1 pc2 pc3 pc4 pc5];
tID=ids;

% add center of mass
cmass=[];
for d=1:length(data)
curr=fliplr(reshape(data(d,:),4,6));
curr=curr';
tmp=find(curr<0.2);
binaryImage = true(size(curr));
binaryImage(tmp)=0;
labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage, curr, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
cmass=[cmass;centerOfMass];
end
pcs=[pcs cmass];


save(fullfile('F:\all_species_temp\newAttempt4\DG_PCanalysis',['PC_F',int2str(fr-1)]),'pcs','tID')
% save(fullfile('F:\all_species_temp\newAttempt2\DG_PCanalysis',['PC_F',int2str(fr-1),'new']),'pcs','tID')
%% version including angle of Gaussian
clear 
close all
path1='F:\all_species_temp\newAttempt\';
load(fullfile(path1,'DG_Gaussian','GaussianNew2'))
keep=[];
for g=1:3
now=Gauss(g,:,5);
now=reshape(now,300,1);
now2=now(1:datalength(g));
keep=[keep;now2];
end
load(fullfile(path1,'DG_PCanalysis','PC_F1new'))
pcs=[pcs keep];
save(fullfile(path1,'DG_PCanalysis','PC_F1newAng'),'pcs','tID')
%% Gaussian of DG heatmaps
clear
close all
plotting=1;
superpath='F:\all_species_temp';
Dtype='wph';
cols='bgr';
st=0;
fact=1.5;

Gauss=zeros(3,300,5);
datalength=zeros(3,1);

for ty=1:3
    load(fullfile(superpath,'newAttempt','DG',[Dtype(ty),'_normPeaks']))
    datalength(ty)=size(normPeaks,1);
    for ce=1:size(normPeaks,1)
        % ce=180;
        
        data=fliplr(reshape(normPeaks(ce,:,2),4,6));
        % data=[0 0 0 0 0 0;0 0 0 0 0 1;0 0 0 0 0 0.8;0 0 0 0 0 0];
        % data=[0 1 0 0 0 0;0 1 2 1 0 0;0 0 1 3 1 0;0 0 0 1 2 1];
        % data=[0.2 0.3 0.25 0.2;0.2 0.5 0.6 0.3;0.3 0.7 1 0.4;0.2 0.35 0.4 0.2];
        
        minData = mean(mean(data))+st*std(reshape(data,24,1));
        p = (data - minData);
        tmp=find(p<0);
        p(tmp)=0;
        sumData = sum(sum(p));
        p=p/sumData;
        
        [y,x] = meshgrid(1:6,1:4);
        xx=reshape(x,24,1);
        yy=reshape(y,24,1);
        pp=reshape(p,24,1);
        mx=dot(xx,pp); %centerpoint
        my=dot(yy,pp); %centerpoint
        
        a=dot((xx-mx).^2,pp);
        b=dot((xx-mx).*(yy-my),pp);
        c=dot((yy-my).^2,pp);
        co=[a,b;b,c];
        
        if plotting==1
            figure
            subplot(2,3,1);imagesc(x(:)',y(:)',data');
            colormap('gray')
            title('original data')
            set(gca,'yticklabel',[4000 2000 1000 500 200 100])
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',[1 2 4 8])
            
            subplot(2,3,2);imagesc(x(:)',y(:)',p');
            colormap('gray')
            title('normalized')
            set(gca,'yticklabel',[4000 2000 1000 500 200 100])
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',[1 2 4 8])
            subplot(2,3,5);imagesc(x(:)',y(:)',p');
            hold on
            colormap('gray')
            title('Gaussian')
            set(gca,'yticklabel',[4000 2000 1000 500 200 100])
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',[1 2 4 8])
        end
        
        if fact>=1
            D=fact*sqrt(eigs(co));
        else
            D=sqrt(eigs(co));%Eigenvalues = two main axes
        end
        [V,dd]=eigs(co);
        ma=max(max(dd));
        [ma1 ma2]=find(dd==ma);
        %         ma1
        V=V(ma1,:); v=V(1)/V(2);
        
        Ang=atan2(V(2),V(1));
        %     Ang=tan(1/v-1); %angle in radians
        %    if V(2)<V(1)
        %         Ang=-Ang;
        %         change=1;
        %    else
        %        change=0;
        %     end
        AngD=rad2deg(Ang); %angle in degrees
        %         if abs(AngD)>180
        %             AngD=abs(AngD)-180;
        %         end
        
        
        centerX=mx;%x-axis center;1=left
        centerY=7-my;%y-axis center; 1=bottom
        
        Gauss(ty,ce,1)=centerX;
        Gauss(ty,ce,2)=centerY;
        Gauss(ty,ce,3)=D(1);
        Gauss(ty,ce,4)=D(2);
        Gauss(ty,ce,5)=Ang;
        
        theta = 0 : 0.01 : 2*pi;
        xEll = D(1)/2 * cos(theta) * cos(Ang) - D(2)/2 * sin(theta) * sin(Ang) + centerX;
        yEll = D(2)/2 * sin(theta) * cos(Ang) + D(1)/2 * cos(theta) * sin(Ang) + centerY;
        
        if plotting==1
            subplot(2,3,4)
            plot(xEll, yEll, 'LineWidth', 2);
            hold on
            plot(centerX,centerY,'ob','markerfacecolor','b','markersize',2)
            text(centerX,centerY,[num2str(round(centerX*100)/100),'/',num2str(round(centerY*100)/100)])
            title([num2str(D(1)), ' / ',num2str(D(2)),' / ',int2str(AngD),'°'])
            axis([0.5 4.5 0.5 6.5])
            
            centerX2=mx;%x-axis center;1=left
            centerY2=my;%y-axis center; 1=bottom
            Ang2=-Ang;
            
            theta = 0 : 0.01 : 2*pi;
            xEll2 = D(1)/2 * cos(theta) * cos(Ang2) - D(2)/2 * sin(theta) * sin(Ang2) + centerX2;
            yEll2 = D(2)/2 * sin(theta) * cos(Ang2) + D(1)/2 * cos(theta) * sin(Ang2) + centerY2;
            subplot(2,3,5)
            plot(xEll2, yEll2, 'LineWidth', 2);
            hold on
            plot(centerX2,centerY2,'ob','markerfacecolor','b','markersize',2)
            set(gcf,'position',[1601 1 1600 824])
            saveas(gcf,fullfile(superpath,'newAttempt','DG_Gaussian',[Dtype(ty),'_',int2str(togo(ce)),'_Gauss.bmp']))
            close
            % [m1 m2]=min(xEll);
            % plot([centerX m1],[centerY,yEll(m2)],'-r')
            % [m3 m4]=max(xEll);
            % m5=abs(m3-m1)/2+min(m3,m1);
            % keep=find(xEll<=m5);
            % tmp=keep;
            % tmp2=find(tmp<=max(m2,m4));
            % tmp=tmp(tmp2);
            % tmp2=find(tmp>=min(m2,m4));
            % tmp=tmp(tmp2);
            % m6=tmp(1);
            % tmp=find(xEll<xEll(m6));
            % tmp2=find(tmp<min(m2,m4));
            % if ~isempty(tmp2)
            %     tmp=tmp(tmp2);
            %     m7=tmp(end);
            % else
            %     tmp2=find(tmp>max(m2,m4));
            %     tmp=tmp(1);
            %     m7=tmp(end);
            % end
            %
            % plot([centerX xEll(m6)],[centerY,yEll(m6)],'-r')
            % %   plot([centerX xEll(m7)],[centerY,yEll(m7)],'-r')
            
            
            
        end
        figure(200)
        subplot(2,2,ty)
        plot(xEll, yEll, 'LineWidth', 1,'color',cols(ty));
        hold on
        plot(centerX,centerY,'o','markerfacecolor',cols(ty),'markeredgecolor',cols(ty),'markersize',2)
        
        
    end
    figure(200)
    subplot(2,2,ty)
    axis([0.5 4.5 0.5 6.5])
    title(Dtype(ty))
    set(gca,'ytick',1:6)
    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',[1 2 4 8])
end
save(fullfile(superpath,'newAttempt','DG_Gaussian','Gaussian'),'Gauss', 'datalength')
set(gcf,'position',[1601 1 1600 824])
saveas(gcf,fullfile(superpath,'newAttempt','DG_Gaussian','Gaussian_overview.bmp'))
saveas(gcf,fullfile(superpath,'newAttempt','DG_Gaussian','Gaussian_overview.fig'))
% close
%% Gaussian of DG heatmaps: better
clear
close all
plotting=0;
superpath='F:\all_species_temp';
Dtype='wph';
cols='bgr';
st=0;
fact=1.5;

Gauss=zeros(3,300,5);
datalength=zeros(3,1);

for ty=1:3
    load(fullfile(superpath,'newAttempt2','DG',[Dtype(ty),'_normPeaks']))
    datalength(ty)=size(normPeaks,1);
    for ce=1:size(normPeaks,1)
        % for ce=1:5;
        
        data=fliplr(reshape(normPeaks(ce,:,2),4,6));
        
        minData = mean(mean(data))+st*std(reshape(data,24,1));
        p = (data - minData);
        tmp=find(p<0);
        p(tmp)=0;
        sumData = sum(sum(p));
        p=p/sumData;
        
        
        
        
        
        if plotting==1
            [y,x] = meshgrid(1:6,1:4);
            figure
            subplot(2,3,1);imagesc(x(:)',y(:)',data');
            colormap('gray')
            title('original data')
            set(gca,'yticklabel',[4000 2000 1000 500 200 100])
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',[1 2 4 8])
            
            subplot(2,3,2);imagesc(x(:)',y(:)',p');
            colormap('gray')
            title('normalized')
            set(gca,'yticklabel',[4000 2000 1000 500 200 100])
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',[1 2 4 8])
            subplot(2,3,5);imagesc(x(:)',y(:)',p');
            hold on
            colormap('gray')
            title('Gaussian')
            set(gca,'yticklabel',[4000 2000 1000 500 200 100])
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',[1 2 4 8])
        end
        
        
        [y,x] = meshgrid(6:-1:1,1:4);
        xx=reshape(x,24,1);
        yy=reshape(y,24,1);
        pp=reshape(p,24,1);
        mx=dot(xx,pp); %centerpoint
        my=dot(yy,pp); %centerpoint
        
        a=dot((xx-mx).^2,pp);
        b=dot((xx-mx).*(yy-my),pp);
        c=dot((yy-my).^2,pp);
        co=[a,b;b,c];
        
        if fact>=1
            D=fact*sqrt(eigs(co));
        else
            D=sqrt(eigs(co));%Eigenvalues = two main axes
        end
        [V,dd]=eigs(co);
        ma=max(max(dd));
        [ma1 ma2]=find(dd==ma);
        %         ma1
        V=V(ma1,:); v=V(1)/V(2);
        
        Ang=atan2(V(2),V(1)); %used
        %                Ang=atan2(V(1),V(2));
        AngD=rad2deg(Ang); %angle in degrees
        
        ymin=get(gca,'ylim');
        ymin=ymin(1);
        
        centerX=mx;%x-axis center;1=left
        centerY=my;%y-axis center; 1=bottom
        
         Gauss(ty,ce,1)=centerX;
        Gauss(ty,ce,2)=centerY;
        Gauss(ty,ce,3)=D(1);
        Gauss(ty,ce,4)=D(2);
        Gauss(ty,ce,5)=Ang;
        
        centerX2=mx;%x-axis center;1=left
        centerY2=7-my;%y-axis center; 1=bottom
        Ang2=-Ang;
        [X,Y]=meshgrid(0.5:0.01:4.5,0.5:0.01:6.5);
        
        theta = pi/2-Ang;
        %
        aa = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
        bb = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
        cc = sin(theta)^2/2/(D(1)/1.5)^2 + cos(theta)^2/2/(D(2)/1.5)^2;
        allClust = exp( - (aa*(Y-centerY).^2 + 2*bb*(Y-centerY).*(X-centerX2) + cc*(X-centerX2).^2)) ;
        if plotting==1
            
            
            subplot(2,3,4)
            
            contour(X,Y,allClust(:,:),1,'color','r')
            %                              contour(X,Y,flipud(allClust(:,:)),1,'color','r')
            %                     plot(xEll2, yEll2, 'LineWidth', 1,'color',col(tID(currID(c))));
            hold on
            plot(centerX,centerY,'or','markerfacecolor','r','markersize',2)
            
            text(centerX,centerY,[num2str(round(centerX*100)/100),'/',num2str(round(centerY*100)/100)])
            title([num2str(D(1)), ' / ',num2str(D(2)),' / ',int2str(AngD),'°'])
            axis([0.5 4.5 0.5 6.5])
            
            subplot(2,3,5)
            contour(X,Y,flipud(allClust(:,:)),1,'color','r')
            
            hold on
            plot(centerX,centerY2,'or','markerfacecolor','r','markersize',2)
            
            %     set(gcf,'position',[1601 1 1600 824])
            %     saveas(gcf,fullfile(superpath,'newAttempt','DG_Gaussian',[Dtype(ty),'_',int2str(togo(ce)),'_Gauss.bmp']))
            %     close
            
        end
        figure(200)
        subplot(2,2,ty)
        contour(X,Y,allClust(:,:),1,'color',cols(ty))
        hold on
        plot(centerX,centerY,'o','markerfacecolor',cols(ty),'markeredgecolor',cols(ty),'markersize',2)
        
        
        [y,x] = meshgrid([4 2 1 0.5 0.2 0.1],[1 2 3 8]);
        xx=reshape(x,24,1);
        yy=reshape(y,24,1);
        pp=reshape(p,24,1);
        mx=dot(xx,pp); %centerpoint
        my=dot(yy,pp); %centerpoint
        
        a=dot((xx-mx).^2,pp);
        b=dot((xx-mx).*(yy-my),pp);
        c=dot((yy-my).^2,pp);
        co=[a,b;b,c];
        
        if fact>=1
            D=fact*sqrt(eigs(co));
        else
            D=sqrt(eigs(co));%Eigenvalues = two main axes
        end
        [V,dd]=eigs(co);
        ma=max(max(dd));
        [ma1 ma2]=find(dd==ma);
        %         ma1
        V=V(ma1,:); v=V(1)/V(2);
        
        Ang=atan2(V(2),V(1));
        AngD=rad2deg(Ang); %angle in degrees
        
        centerX2=mx;%x-axis center;1=left
        centerY=my;%y-axis center; 1=bottom
        
        [X,Y]=meshgrid(0.0:0.02:8.5,0.0:0.010:4.100);
        
        
        centerX=mx;%x-axis center;1=left
        centerY=my*1000;%y-axis center; 1=bottom
        
        Gauss2(ty,ce,1)=centerX;
        Gauss2(ty,ce,2)=centerY;
        Gauss2(ty,ce,3)=D(1);
        Gauss2(ty,ce,4)=D(2)*1000;
        Gauss2(ty,ce,5)=Ang;
    end
    figure(200)
    subplot(2,2,ty)
    axis([0.5 4.5 0.5 6.5])
    title(Dtype(ty))
    set(gca,'ytick',1:6)
    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',[1 2 4 8])
end

save(fullfile(superpath,'newAttempt2','DG_Gaussian','GaussianNew2'),'Gauss', 'datalength','Gauss2')
set(gcf,'position',[1601 1 1600 824])
saveas(gcf,fullfile(superpath,'newAttempt2','DG_Gaussian','Gaussian_overviewNew.bmp'))
saveas(gcf,fullfile(superpath,'newAttempt2','DG_Gaussian','Gaussian_overviewNew.fig'))
% close
%% prepare F1 amplitudes for clustering
clear
close all

fr=2;
ids=[];
load(fullfile('F:\all_species_temp\newAttempt2\DG','w_normPeaks'))
data1=normPeaks(:,:,fr);
% data1=fliplr(reshape(data1,size(data1,1),4,6));
ids=[ids;repmat(1,length(data1),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','p_normPeaks'))
data2=normPeaks(:,:,fr);
% data2=fliplr(reshape(data2,size(data2,1),4,6));
ids=[ids;repmat(2,length(data2),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','h_normPeaks'))
data3=normPeaks(:,:,fr);
% data3=fliplr(reshape(data3,size(data3,1),4,6));
ids=[ids;repmat(3,length(data3),1)];

data=cat(1,data1,data2,data3);
data=data(:,5:end);
amps=data;
tID=ids;
save(fullfile('F:\all_species_temp\newAttempt2\DG_PCanalysis',['amp_F',int2str(fr-1)]),'amps','tID')
%% prepare F1 amplitudes for clustering: only Human
clear
close all

fr=2;
ids=[];
load(fullfile('e:\all_species_temp\newAttempt2\DG','h_normPeaks'))
data=normPeaks(:,:,fr);
% data1=fliplr(reshape(data1,size(data1,1),4,6));
ids=[ids;repmat(3,length(data),1)];

data=data(:,5:end);
amps=data;
tID=ids;
save(fullfile('e:\all_species_temp\newAttemptOnlyHum\DG_PCanalysis',['amp_F',int2str(fr-1)]),'amps','tID')
%% prepare other parameters for clustering
clear
close all

fr=2;
ids=[];
load(fullfile('F:\all_species_temp\newAttempt2\DG','w_normPeaks'))
data1=normPeaks(:,:,fr);
% data1=fliplr(reshape(data1,size(data1,1),4,6));
ids=[ids;repmat(1,length(data1),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','p_normPeaks'))
data2=normPeaks(:,:,fr);
% data2=fliplr(reshape(data2,size(data2,1),4,6));
ids=[ids;repmat(2,length(data2),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','h_normPeaks'))
data3=normPeaks(:,:,fr);
% data3=fliplr(reshape(data3,size(data3,1),4,6));
ids=[ids;repmat(3,length(data3),1)];

data=cat(1,data1,data2,data3);

%spat

spat=[];
for s=1:4:21
me=mean(data(:,s:s+3)');
spat=[spat me'];
end
spat=max(spat')';


%temp
temp=[];
for t=1:4
me=mean(data(:,t:4:24)');
temp=[temp me'];
end
temp=max(temp')';

%speed
starts=[3 4 8 12 16 20];
le=[2 3 3 3 2 1];
speed=[]; scnt=0;
for s=starts
    scnt=scnt+1;
me=mean(data(:,s:3:s+le(scnt)*3)');
speed=[speed me'];
end
speed=max(speed')';

%CoM
cmass=[];
for d=1:length(data)
curr=fliplr(reshape(data(d,:),4,6));
curr=curr';
tmp=find(curr<0.2);
binaryImage = true(size(curr));
binaryImage(tmp)=0;
labeledImage = bwlabel(binaryImage);
measurements = regionprops(labeledImage, curr, 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
cmass=[cmass;centerOfMass];
end

su=sum(data')';

par=[spat temp speed cmass ];
% par=[spat temp speed su cmass ];

tID=ids;
save(fullfile('F:\all_species_temp\newAttempt5\DG_PCanalysis',['par_F',int2str(fr-1),'new']),'par','tID')
%% prepare other parameters for clustering 2
clear
close all

fr=2;
ids=[];
load(fullfile('F:\all_species_temp\newAttempt2\DG','w_normPeaks'))
data1=normPeaks(:,:,fr);
% data1=fliplr(reshape(data1,size(data1,1),4,6));
ids=[ids;repmat(1,length(data1),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','p_normPeaks'))
data2=normPeaks(:,:,fr);
% data2=fliplr(reshape(data2,size(data2,1),4,6));
ids=[ids;repmat(2,length(data2),1)];
load(fullfile('F:\all_species_temp\newAttempt2\DG','h_normPeaks'))
data3=normPeaks(:,:,fr);
% data3=fliplr(reshape(data3,size(data3,1),4,6));
ids=[ids;repmat(3,length(data3),1)];

data=cat(1,data1,data2,data3);

%spat

spat=[];
for s=1:4:21
me=mean(data(:,s:s+3)');
spat=[spat me'];
end
spat=spat(:,2:end);
% spat=spat./repmat(max(spat')',1,5);


%temp
temp=[];
for t=1:4
me=mean(data(:,t:4:24)');
temp=[temp me'];
end
% temp=temp./repmat(max(temp')',1,4);

%speed
starts=[3 4 8 12 16 20];
le=[2 3 3 3 2 1];
speed=[]; scnt=0;
for s=starts
    scnt=scnt+1;
me=mean(data(:,s:3:s+le(scnt)*3)');
speed=[speed me'];
end
speed=max(speed')';

su=sum(data')';


means=[spat temp speed su];
means=[spat temp speed ]; %new

tID=ids;
save(fullfile('F:\all_species_temp\newAttempt6\DG_PCanalysis',['means_F',int2str(fr-1),'new']),'means','tID')
%% !!clustering kmeans
clear
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% type='GaussianNew2';
% path1='F:\all_species_temp\newAttempt\DG_PCanalysis';
type='PC_F1';
% path1='F:\all_species_temp\newAttempt6\DG_PCanalysis';
path1='E:\all_species_temp\newAttemptOnlyMouse\DG_PCanalysis';
% type='means_F1';

outlier=1;

load(fullfile(path1,type))
savename='clust_Kmeans_';
close all

orig=pcs;
gau=0;
% orig=Gauss;
% gau=1;
% orig=amps;
% gau=0;

if gau==1
    
    orig2clust=[]; tID=zeros(sum(datalength),1);
    last=1;
    for i=1:3
        curr=reshape(Gauss(i,:,:),300,5); %centerX, centerY, x-axis, y-axis, angle
        curr=curr(1:datalength(i),:);
        orig2clust=[orig2clust;curr];
        tID(last:last+datalength(i)-1)=repmat(i,datalength(i),1);
        last=last+datalength(i);
    end
    tmp=find(orig2clust(:,5)<0);
    tmp2=abs(-deg2rad(180)-orig2clust(tmp,5));
    orig2clust(tmp,5)=tmp2;
    keepOrig=orig2clust;
else
    orig2clust=orig;
end
if outlier==1
    for p=1:size(orig2clust,2)-1
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
    end
curr=orig2clust(:,end);
[tmp tmp2]=sort(curr);
norm=(curr-min(curr))/(tmp(end-1)+0.2-min(curr));
norm(tmp2(end))=1;
orig2clust(:,end)=norm;
else
for p=1:size(orig2clust,2)
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
end
end
opts=statset('maxiter',1000);
ID=[];
for k=18:30
    [idx c]=kmeans(orig2clust,k,'emptyaction','drop','start','cluster','replicates',10000,'options',opts);
    k
    ID=[ID idx];
end

if gau==1
    ang=rad2deg(keepOrig(:,5));
    keepOrig(:,5)=ang;
else
    keepOrig=orig2clust;
end
save(fullfile(path1,[savename,type]),'orig2clust','keepOrig','ID','tID')
%% !!clustering fuzzy c-means
clear
% path1='F:\all_species_temp\newAttempt3\DG_Gaussian';
% type='GaussianNew2';
% path1='F:\all_species_temp\newAttempt\DG_PCanalysis';
type='means_F1new';
path1='F:\all_species_temp\newAttempt6\DG_PCanalysis';
% type='amp_F1';

load(fullfile(path1,type))
savename='clust_Cmeans_';
close all
outlier=0;
orig=means;
gau=0;
% orig=Gauss;
% gau=1;
% orig=amps;
% gau=0;

if gau==1
    
    orig2clust=[]; tID=zeros(sum(datalength),1);
    last=1;
    for i=1:3
        curr=reshape(Gauss(i,:,:),300,5); %centerX, centerY, x-axis, y-axis, angle
        curr=curr(1:datalength(i),:);
        orig2clust=[orig2clust;curr];
        tID(last:last+datalength(i)-1)=repmat(i,datalength(i),1);
        last=last+datalength(i);
    end
    tmp=find(orig2clust(:,5)<0);
    tmp2=abs(-deg2rad(180)-orig2clust(tmp,5));
    orig2clust(tmp,5)=tmp2;
    keepOrig=orig2clust;
else
    orig2clust=orig;
end
if outlier==1
    for p=1:size(orig2clust,2)-1
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
    end
curr=orig2clust(:,end);
[tmp tmp2]=sort(curr);
norm=(curr-min(curr))/(tmp(end-1)+0.2-min(curr));
norm(tmp2(end))=1;
orig2clust(:,end)=norm;
else
for p=1:size(orig2clust,2)
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
end
end

opts=[2;10000;nan;0];
ID=[];
for k=2:30
    [c,idx,obj_fcn]=fcm(orig2clust,k,opts);
 [tmp pos]=max(idx);
    k
    ID=[ID pos'];
end

if gau==1
    ang=rad2deg(keepOrig(:,5));
    keepOrig(:,5)=ang;
else
    keepOrig=orig2clust;
end
save(fullfile(path1,[savename,type]),'orig2clust','keepOrig','ID','tID')
%% !!clustering fuzzy GK
clear
% path1='F:\all_species_temp\newAttempt3\DG_Gaussian';
% type='GaussianNew2';
% path1='F:\all_species_temp\newAttempt\DG_PCanalysis';
type='means_F1new';
path1='F:\all_species_temp\newAttempt6\DG_PCanalysis';
% type='amp_F1';

addpath(genpath('C:\Program Files\MATLAB\R2014a\toolbox\FUZZCLUST'))
load(fullfile(path1,type))
savename='clust_fuzzyGK_';
close all
outlier=0;
orig=means;
gau=0;
% orig=Gauss;
% gau=1;
% orig=amps;
% gau=0;

if gau==1
    
    orig2clust=[]; tID=zeros(sum(datalength),1);
    last=1;
    for i=1:3
        curr=reshape(Gauss(i,:,:),300,5); %centerX, centerY, x-axis, y-axis, angle
        curr=curr(1:datalength(i),:);
        orig2clust=[orig2clust;curr];
        tID(last:last+datalength(i)-1)=repmat(i,datalength(i),1);
        last=last+datalength(i);
    end
    tmp=find(orig2clust(:,5)<0);
    tmp2=abs(-deg2rad(180)-orig2clust(tmp,5));
    orig2clust(tmp,5)=tmp2;
    keepOrig=orig2clust;
else
    orig2clust=orig;
end
if outlier==1
    for p=1:size(orig2clust,2)-1
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
    end
curr=orig2clust(:,end);
[tmp tmp2]=sort(curr);
norm=(curr-min(curr))/(tmp(end-1)+0.2-min(curr));
norm(tmp2(end))=1;
orig2clust(:,end)=norm;
else
for p=1:size(orig2clust,2)
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
end
end

ID=[]; valids=[];
for c=2:30
param.c=c;
param.m=2;
param.e=1e-6;
param.ro=ones(1,param.c);
param.val=3;
data.X=orig2clust;
result = GKclust(data,param);
[ma id]=max(result.data.f');
ID=[ID id'];
%validation
va=[];
param.val=2;
val = validity(result,data,param);
va=[va real(val.validity.SC) real(val.validity.S) real(val.validity.XB)];
param.val=3;
val = validity(result,data,param);
va=[va val.validity.DI val.validity.ADI];
valids=[valids va'];
end

if gau==1
    ang=rad2deg(keepOrig(:,5));
    keepOrig(:,5)=ang;
else
    keepOrig=orig2clust;
end
save(fullfile(path1,[savename,type]),'orig2clust','keepOrig','ID','tID','valids')
%% plot clusters
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
% type='GaussianNew2_corr';
version=1;
type='PC_F1new';
% clustM='clust_MixGau_';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))


% sorted=Tmat;
sorted=ID;
% sorted=data;
gau=0;

clust2plot=14;

% if version==2
%      keepOrig(:,2)=keepOrig(:,2)/1000;
%                      keepOrig(:,4)=keepOrig(:,4)/1000;
% end

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0];
col='bgr';
col2=[0 0 1;0 1 0;1 0 0];
theta = 0 : 0.01 : 2*pi;


for i=clust2plot-1
    
    u=unique(sorted(:,i));
    %     if gau==1
    subs=length(u)+1;
    %     else
    %         subs=length(u);
    %     end
    if subs<=2
        x=1; y=2;
    elseif subs<=4
        x=2; y=2;
    elseif subs<=6
        x=2; y=3;
    elseif subs<=9
        x=3; y=3;
    elseif subs<=12
        x=3; y=4;
    elseif subs<=15
        x=3; y=5;
    elseif subs<=20
        x=4;y=5;
    elseif subs<=25
        x=5;y=5;
    else
        x=5; y=6;
    end
    for q=1:length(u)
        currID=find(sorted(:,i)==u(q));
        
        if gau==1
            for c=1:length(currID)
                
                
                if version==2
                    
                    [X,Y]=meshgrid(0.0:0.02:8.5,0.0:0.010:4.100);
                    theta=pi/2-deg2rad(keepOrig(currID(c),5));
                    aa = cos(theta)^2/2/(keepOrig(currID(c),3)/1.5)^2 + sin(theta)^2/2/(keepOrig(currID(c),4)/1000/1.5)^2;
                    bb = -sin(2*theta)/4/(keepOrig(currID(c),3)/1.5)^2 + sin(2*theta)/4/(keepOrig(currID(c),4)/1000/1.5)^2 ;
                    cc = sin(theta)^2/2/(keepOrig(currID(c),3)/1.5)^2 + cos(theta)^2/2/(keepOrig(currID(c),4)/1000/1.5)^2;
                    allClust = exp( - (aa*(Y-keepOrig(currID(c),2)/1000).^2 + 2*bb*(Y-keepOrig(currID(c),2)/1000).*(X-keepOrig(currID(c),1)) + cc*(X-keepOrig(currID(c),1)).^2)) ;
                    
                    [X,Y]=meshgrid(0.0:0.02:8.5,0.0:10:4100);
                else
                    [X,Y]=meshgrid(0:0.02:5,0:0.02:7);
                    theta = -deg2rad(keepOrig(currID(c),5));
                    
                    aa = cos(theta)^2/2/(keepOrig(currID(c),3)/1.5)^2 + sin(theta)^2/2/(keepOrig(currID(c),4)/1.5)^2;
                    bb = -sin(2*theta)/4/(keepOrig(currID(c),3)/1.5)^2 + sin(2*theta)/4/(keepOrig(currID(c),4)/1.5)^2 ;
                    cc = sin(theta)^2/2/(keepOrig(currID(c),3)/1.5)^2 + cos(theta)^2/2/(keepOrig(currID(c),4)/1.5)^2;
                    allClust = exp( - (aa*(X-keepOrig(currID(c),1)).^2 + 2*bb*(X-keepOrig(currID(c),1)).*(Y-keepOrig(currID(c),2)) + cc*(Y-keepOrig(currID(c),2)).^2)) ;
                    
                end
                
                
                %
                %                 xEll2 = keepOrig(currID(c),3)/2 * cos(theta) * cos(keepOrig(currID(c),5)) - keepOrig(currID(c),4)/2 * sin(theta) * sin(keepOrig(currID(c),5)) + keepOrig(currID(c),1);
                %                 yEll2 = keepOrig(currID(c),4)/2 * sin(theta) * cos(keepOrig(currID(c),5)) + keepOrig(currID(c),3)/2 * cos(theta) * sin(keepOrig(currID(c),5)) + keepOrig(currID(c),2);
                
                if plotting==1
                    figure(i)
                    subplot(x,y,q)
                    contour(X,Y,allClust(:,:),1,'color',col(tID(currID(c))))
                    %                     plot(xEll2, yEll2, 'LineWidth', 1,'color',col(tID(currID(c))));
                    hold on
                    
                    plot(keepOrig(currID(c),1),keepOrig(currID(c),2),'o','markerfacecolor',col(tID(currID(c))),'markeredgecolor',col(tID(currID(c))),'markersize',2)
                    
                    
                    subplot(x,y,subs)
                    %                     plot(xEll2, yEll2, 'LineWidth', 1,'color',colors(u(q),:));
                    contour(X,Y,allClust(:,:),1,'color',colors(u(q),:))
                    hold on
                    
                    plot(keepOrig(currID(c),1),keepOrig(currID(c),2),'o','markerfacecolor',colors(u(q),:),'markeredgecolor',colors(u(q),:),'markersize',2)
                    
                    
                end
            end
            if plotting==1
                figure(i)
                subplot(x,y,q)
                %         axis([0.5 4.5 0.5 6.5])
                if version==2
                    axis([0 8.5 0 4100])
                else
                    axis([-0.5 5.5 0.5 7])
                    set(gca,'ytick',1:6)
                    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
                    set(gca,'xtick',1:4)
                    set(gca,'xticklabel',[1 2 4 8])
                end
            end
            
            NofC=zeros(3,1);
            for t=1:3
                tmp=find(tID(currID)==t);
                NofC(t)=length(tmp);
                
                if plotting==1
                    figure(100+i)
                    subplot(x,y,q)
                    
                    scatter(keepOrig(currID(tmp),3),keepOrig(currID(tmp),4),col(t))
                    hold on
                    plot(mean(keepOrig(currID(tmp),3)),mean(keepOrig(currID(tmp),4)),'*k','markersize',5)
                    
                    xlabel('x-axis of Gauss')
                    ylabel('y-axis of Gauss')
                    if version==2
                        plot([0 8],[0 4000],'-','color',[0.6 0.6 0.6])
                        axis([0 9 0 4100])
                    else
                        plot([0 2],[0 2],'-','color',[0.6 0.6 0.6])
                        axis([0 3 0 2])
                    end
                    title(['cluster ',int2str(u(q))])
                    
                    figure(1000+i)
                    subplot(x,y,q)
                    scatter(keepOrig(currID(tmp),1),keepOrig(currID(tmp),5),col(t))
                    hold on
                    plot(mean(keepOrig(currID(tmp),1)),mean(keepOrig(currID(tmp),5)),'*k','markersize',5)
                    plot([0.5 8.5],[90 90],'-','color',[0.6 0.6 0.6])
                    plot([0.5 8.5],[45 45],'-','color',[0.6 0.6 0.6])
                    plot([0.5 8.5],[135 135],'-','color',[0.6 0.6 0.6])
                    xlabel('x-center')
                    ylabel('angle')
                    
                    title(['cluster ',int2str(u(q))])
                    if version==2
                        axis([0.5 9.5 0 180])
                    else
                        set(gca,'xtick',1:4)
                        set(gca,'xticklabel',[1 2 4 8])
                        axis([0.5 4.5 0 180])
                    end
                    
                    figure(10000+i)
                    subplot(x,y,q)
                    scatter(keepOrig(currID(tmp),5),keepOrig(currID(tmp),2),col(t))
                    hold on
                    plot(mean(keepOrig(currID(tmp),5)),mean(keepOrig(currID(tmp),2)),'*k','markersize',5)
                    
                    ylabel('y-center')
                    xlabel('angle')
                    title(['cluster ',int2str(u(q))])
                    
                    if version==2
                        plot([90 90],[0 4100],'-','color',[0.6 0.6 0.6])
                        plot([45 45],[0 4100],'-','color',[0.6 0.6 0.6])
                        plot([135 135],[0 4100],'-','color',[0.6 0.6 0.6])
                        axis([0 180 0 4100 ])
                    else
                        plot([90 90],[0.5 6.5],'-','color',[0.6 0.6 0.6])
                        plot([45 45],[0.5 6.5],'-','color',[0.6 0.6 0.6])
                        plot([135 135],[0.5 6.5],'-','color',[0.6 0.6 0.6])
                        set(gca,'ytick',1:6)
                        set(gca,'yticklabel',[100 200 500 1000 2000 4000])
                        axis([0 180 0.5 6.5 ])
                    end
                    
                    
                end
            end
            if plotting==1
                figure(i)
                subplot(x,y,q)
                title(['cluster ',int2str(u(q)),' (m=',int2str(NofC(1)),'; ',int2str(100/length(find(tID==1))*NofC(1)),'% / p=',int2str(NofC(2)),'; ',int2str(100/length(find(tID==2))*NofC(2)),'% / h=',int2str(NofC(3)),'; ',int2str(100/length(find(tID==3))*NofC(3)),'%)'])
                
                
                figure(100000+i)
                subplot(x,y,q)
                %                 xEll3 = mean(keepOrig(currID,3))/2 * cos(theta) * cos(mean(keepOrig(currID,5))) - mean(keepOrig(currID,4))/2 * sin(theta) * sin(mean(keepOrig(currID,5))) + mean(keepOrig(currID,1));
                %                 yEll3 = mean(keepOrig(currID,4))/2 * sin(theta) * cos(mean(keepOrig(currID,5))) + mean(keepOrig(currID,3))/2 * cos(theta) * sin(mean(keepOrig(currID,5))) + mean(keepOrig(currID,2));
                %                 plot(xEll3, yEll3, 'LineWidth', 1,'color',colors(q,:));
                if version==2
                    [X,Y]=meshgrid(0.0:0.02:8.5,0.0:0.010:4.100);
                    theta = pi/2-deg2rad(mean(keepOrig(currID,5)));
                    aa = cos(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + sin(theta)^2/2/(mean(keepOrig(currID,4))/1000/1.5)^2;
                    bb = -sin(2*theta)/4/(mean(keepOrig(currID,3))/1.5)^2 + sin(2*theta)/4/(mean(keepOrig(currID,4))/1000/1.5)^2 ;
                    cc = sin(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + cos(theta)^2/2/(mean(keepOrig(currID,4))/1000/1.5)^2;
                    allClust = exp( - (aa*(Y-mean(keepOrig(currID,2))/1000).^2 + 2*bb*(Y-mean(keepOrig(currID,2))/1000).*(X-mean(keepOrig(currID,1))) + cc*(X-mean(keepOrig(currID,1))).^2)) ;
                    [X,Y]=meshgrid(0.0:0.02:8.5,0.0:10:4100);
                else
                    [X,Y]=meshgrid(0:0.02:5,0:0.02:7);
                    theta = -deg2rad(mean(keepOrig(currID,5)));
                    aa = cos(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + sin(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
                    bb = -sin(2*theta)/4/(mean(keepOrig(currID,3))/1.5)^2 + sin(2*theta)/4/(mean(keepOrig(currID,4))/1.5)^2 ;
                    cc = sin(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + cos(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
                    allClust = exp( - (aa*(X-mean(keepOrig(currID,1))).^2 + 2*bb*(X-mean(keepOrig(currID,1))).*(Y-mean(keepOrig(currID,2))) + cc*(Y-mean(keepOrig(currID,2))).^2)) ;
                end
                
                contour(X,Y,allClust(:,:),1,'color',colors(q,:))
                hold on
                plot(mean(keepOrig(currID,1)),mean(keepOrig(currID,2)),'o','markerfacecolor',colors(q,:),'markeredgecolor',colors(q,:),'markersize',2)
                title(['cluster ',int2str(u(q)),': av. Gaussian'])
                if version==2
                    axis([0 9 0 4100])
                else
                    axis([-0.5 5.5 0.5 6.5])
                    set(gca,'ytick',1:6)
                    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
                    set(gca,'xtick',1:4)
                    set(gca,'xticklabel',[1 2 4 8])
                end
                
                subplot(x,y,subs)
                contour(X,Y,allClust(:,:),1,'color',colors(q,:))
                hold on
                plot(mean(keepOrig(currID,1)),mean(keepOrig(currID,2)),'o','markerfacecolor',colors(q,:),'markeredgecolor',colors(q,:),'markersize',2)
                title('av. Gaussian')
                if version==2
                    axis([0 9 0 4100])
                else
                    axis([-0.5 5.5 0.5 6.5])
                    set(gca,'ytick',1:6)
                    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
                    set(gca,'xtick',1:4)
                    set(gca,'xticklabel',[1 2 4 8])
                end
                
            end
        else %PC
             NofC=zeros(3,1);
            for t=1:3
                tmp=find(tID(currID)==t);
                NofC(t)=length(tmp);
            end
            if plotting==1
                figure(i)
                subplot(x,y,q)
                scatter(keepOrig(currID,1),keepOrig(currID,2), 15,col2(tID(currID),:));
                hold on
                xlabel('PC1')
                ylabel('PC2')
                title(['cluster ',int2str(u(q)),' (m=',int2str(NofC(1)),'; ',int2str(100/length(find(tID==1))*NofC(1)),'% / p=',int2str(NofC(2)),'; ',int2str(100/length(find(tID==2))*NofC(2)),'% / h=',int2str(NofC(3)),'; ',int2str(100/length(find(tID==3))*NofC(3)),'%)'])
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,2))-0.1 max(keepOrig(:,2))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,1),keepOrig(currID,2), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,2))-0.1 max(keepOrig(:,2))+0.1])
                figure(100+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,1),keepOrig(currID,3), 15,col2(tID(currID),:));
                hold on
                xlabel('PC1')
                ylabel('PC3')
                title(['cluster ',int2str(u(q))])
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,3))-0.1 max(keepOrig(:,3))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,1),keepOrig(currID,3), 15,col2(tID(currID),:));
                hold on
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,3))-0.1 max(keepOrig(:,3))+0.1])
                figure(1000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,1),keepOrig(currID,4), 15,col2(tID(currID),:));
                hold on
                xlabel('PC1')
                ylabel('PC4')
                title(['cluster ',int2str(u(q))])
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,4))-0.1 max(keepOrig(:,4))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,1),keepOrig(currID,4), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,4))-0.1 max(keepOrig(:,4))+0.1])
                figure(10000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,1),keepOrig(currID,5), 15,col2(tID(currID),:));
                hold on
                xlabel('PC1')
                ylabel('PC5')
                title(['cluster ',int2str(u(q))])
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,5))-0.1 max(keepOrig(:,5))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,1),keepOrig(currID,5), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,5))-0.1 max(keepOrig(:,5))+0.1])
                figure(100000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,2),keepOrig(currID,3), 15,col2(tID(currID),:));
                hold on
                xlabel('PC2')
                ylabel('PC3')
                title(['cluster ',int2str(u(q))])
                axis([min(keepOrig(:,2))-0.1 max(keepOrig(:,2))+0.1 min(keepOrig(:,3))-0.1 max(keepOrig(:,3))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,2),keepOrig(currID,3), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,2))-0.1 max(keepOrig(:,2))+0.1 min(keepOrig(:,3))-0.1 max(keepOrig(:,3))+0.1])
                figure(1000000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,4),keepOrig(currID,5), 15,col2(tID(currID),:));
                hold on
                xlabel('PC4')
                ylabel('PC5')
                title(['cluster ',int2str(u(q))])
                axis([min(keepOrig(:,4))-0.1 max(keepOrig(:,4))+0.1 min(keepOrig(:,5))-0.1 max(keepOrig(:,5))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,4),keepOrig(currID,5), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,4))-0.1 max(keepOrig(:,4))+0.1 min(keepOrig(:,5))-0.1 max(keepOrig(:,5))+0.1])
                if size(keepOrig,2)>5
                     figure(10000000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,6),keepOrig(currID,7), 15,col2(tID(currID),:));
                hold on
                xlabel('CoM (x-coord)')
                ylabel('CoM (y-coord)')
                title(['cluster ',int2str(u(q))])
                axis([min(keepOrig(:,6))-0.1 max(keepOrig(:,6))+0.1 min(keepOrig(:,7))-0.1 max(keepOrig(:,7))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,6),keepOrig(currID,7), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,6))-0.1 max(keepOrig(:,6))+0.1 min(keepOrig(:,7))-0.1 max(keepOrig(:,7))+0.1])
                end
            end
        end
        
    end
    if plotting==1 && gau==1
        figure(i)
        subplot(x,y,subs)
        title('overview')
        if version==2
            axis([0 9 0 4100])
        else
            axis([0.5 4.5 0.5 6.5])
            set(gca,'ytick',1:6)
            set(gca,'yticklabel',[100 200 500 1000 2000 4000])
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',[1 2 4 8])
        end
        
    end
end
%% plot clusters: after decision
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
% type='GaussianNew2_corr';
version=1;
type='PC_F1new_sel';
% clustM='clust_MixGau_';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))

ctype=1;

sorted=ID;

gau=0;


colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0];
col='bgr';
col2=[0 0 1;0 1 0;1 0 0];
theta = 0 : 0.01 : 2*pi;

clusteringResult=[1 11 12;2 9 5;3 14 13;4 7 2;5 10 1;6 13 9;7 8 7;8 6 10;9 8 14;10 5 3;11 12 4;12 3 11;13 5 8;14 10 3];

   clusts=clusteringResult(:,ctype);
   u=unique(clusts);
   missing=1:max(sorted);
   missing=setdiff(missing,u);
   clusts=[clusts;missing'];
   cluname=[1:max(sorted) repmat(0,1,length(missing))];
i=1;
%     subs=length(clusts)+1;
subs=16;

    if subs<=2
        x=1; y=2;
    elseif subs<=4
        x=2; y=2;
    elseif subs<=6
        x=2; y=3;
    elseif subs<=9
        x=3; y=3;
    elseif subs<=12
        x=3; y=4;
    elseif subs<=15
        x=3; y=5;
    elseif subs<=20
        x=4;y=5;
    elseif subs<=25
        x=5;y=5;
    else
        x=5; y=6;
    end
    for q=1:length(clusts)
        currID=find(sorted==clusts(q));
        
        if gau==1
            for c=1:length(currID)
                
                
                if version==2
                    
                    [X,Y]=meshgrid(0.0:0.02:8.5,0.0:0.010:4.100);
                    theta=pi/2-deg2rad(keepOrig(currID(c),5));
                    aa = cos(theta)^2/2/(keepOrig(currID(c),3)/1.5)^2 + sin(theta)^2/2/(keepOrig(currID(c),4)/1000/1.5)^2;
                    bb = -sin(2*theta)/4/(keepOrig(currID(c),3)/1.5)^2 + sin(2*theta)/4/(keepOrig(currID(c),4)/1000/1.5)^2 ;
                    cc = sin(theta)^2/2/(keepOrig(currID(c),3)/1.5)^2 + cos(theta)^2/2/(keepOrig(currID(c),4)/1000/1.5)^2;
                    allClust = exp( - (aa*(Y-keepOrig(currID(c),2)/1000).^2 + 2*bb*(Y-keepOrig(currID(c),2)/1000).*(X-keepOrig(currID(c),1)) + cc*(X-keepOrig(currID(c),1)).^2)) ;
                    
                    [X,Y]=meshgrid(0.0:0.02:8.5,0.0:10:4100);
                else
                    [X,Y]=meshgrid(0:0.02:5,0:0.02:7);
                    theta = -deg2rad(keepOrig(currID(c),5));
                    
                    aa = cos(theta)^2/2/(keepOrig(currID(c),3)/1.5)^2 + sin(theta)^2/2/(keepOrig(currID(c),4)/1.5)^2;
                    bb = -sin(2*theta)/4/(keepOrig(currID(c),3)/1.5)^2 + sin(2*theta)/4/(keepOrig(currID(c),4)/1.5)^2 ;
                    cc = sin(theta)^2/2/(keepOrig(currID(c),3)/1.5)^2 + cos(theta)^2/2/(keepOrig(currID(c),4)/1.5)^2;
                    allClust = exp( - (aa*(X-keepOrig(currID(c),1)).^2 + 2*bb*(X-keepOrig(currID(c),1)).*(Y-keepOrig(currID(c),2)) + cc*(Y-keepOrig(currID(c),2)).^2)) ;
                    
                end
                
                
                %
                %                 xEll2 = keepOrig(currID(c),3)/2 * cos(theta) * cos(keepOrig(currID(c),5)) - keepOrig(currID(c),4)/2 * sin(theta) * sin(keepOrig(currID(c),5)) + keepOrig(currID(c),1);
                %                 yEll2 = keepOrig(currID(c),4)/2 * sin(theta) * cos(keepOrig(currID(c),5)) + keepOrig(currID(c),3)/2 * cos(theta) * sin(keepOrig(currID(c),5)) + keepOrig(currID(c),2);
                
                if plotting==1
                    figure(i)
                    subplot(x,y,q)
                    contour(X,Y,allClust(:,:),1,'color',col(tID(currID(c))))
                    %                     plot(xEll2, yEll2, 'LineWidth', 1,'color',col(tID(currID(c))));
                    hold on
                    
                    plot(keepOrig(currID(c),1),keepOrig(currID(c),2),'o','markerfacecolor',col(tID(currID(c))),'markeredgecolor',col(tID(currID(c))),'markersize',2)
                    
                    
                    subplot(x,y,subs)
                    %                     plot(xEll2, yEll2, 'LineWidth', 1,'color',colors(u(q),:));
                    contour(X,Y,allClust(:,:),1,'color',colors(q,:))
                    hold on
                    
                    plot(keepOrig(currID(c),1),keepOrig(currID(c),2),'o','markerfacecolor',colors(q,:),'markeredgecolor',colors(q,:),'markersize',2)
                    
                    
                end
            end
            if plotting==1
                figure(i)
                subplot(x,y,q)
                %         axis([0.5 4.5 0.5 6.5])
                if version==2
                    axis([0 8.5 0 4100])
                else
                    axis([-0.5 5.5 0.5 7])
                    set(gca,'ytick',1:6)
                    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
                    set(gca,'xtick',1:4)
                    set(gca,'xticklabel',[1 2 4 8])
                end
            end
            
            NofC=zeros(3,1);
            for t=1:3
                tmp=find(tID(currID)==t);
                NofC(t)=length(tmp);
                
                if plotting==1
                    figure(100+i)
                    subplot(x,y,q)
                    
                    scatter(keepOrig(currID(tmp),3),keepOrig(currID(tmp),4),col(t))
                    hold on
                    plot(mean(keepOrig(currID(tmp),3)),mean(keepOrig(currID(tmp),4)),'*k','markersize',5)
                    
                    xlabel('x-axis of Gauss')
                    ylabel('y-axis of Gauss')
                    if version==2
                        plot([0 8],[0 4000],'-','color',[0.6 0.6 0.6])
                        axis([0 9 0 4100])
                    else
                        plot([0 2],[0 2],'-','color',[0.6 0.6 0.6])
                        axis([0 3 0 2])
                    end
                    title(['cluster ',int2str(cluname(q))])
                    
                    figure(1000+i)
                    subplot(x,y,q)
                    scatter(keepOrig(currID(tmp),1),keepOrig(currID(tmp),5),col(t))
                    hold on
                    plot(mean(keepOrig(currID(tmp),1)),mean(keepOrig(currID(tmp),5)),'*k','markersize',5)
                    plot([0.5 8.5],[90 90],'-','color',[0.6 0.6 0.6])
                    plot([0.5 8.5],[45 45],'-','color',[0.6 0.6 0.6])
                    plot([0.5 8.5],[135 135],'-','color',[0.6 0.6 0.6])
                    xlabel('x-center')
                    ylabel('angle')
                    
                    title(['cluster ',int2str(cluname(q))])
                    if version==2
                        axis([0.5 9.5 0 180])
                    else
                        set(gca,'xtick',1:4)
                        set(gca,'xticklabel',[1 2 4 8])
                        axis([0.5 4.5 0 180])
                    end
                    
                    figure(10000+i)
                    subplot(x,y,q)
                    scatter(keepOrig(currID(tmp),5),keepOrig(currID(tmp),2),col(t))
                    hold on
                    plot(mean(keepOrig(currID(tmp),5)),mean(keepOrig(currID(tmp),2)),'*k','markersize',5)
                    
                    ylabel('y-center')
                    xlabel('angle')
                    title(['cluster ',int2str(cluname(q))])
                    
                    if version==2
                        plot([90 90],[0 4100],'-','color',[0.6 0.6 0.6])
                        plot([45 45],[0 4100],'-','color',[0.6 0.6 0.6])
                        plot([135 135],[0 4100],'-','color',[0.6 0.6 0.6])
                        axis([0 180 0 4100 ])
                    else
                        plot([90 90],[0.5 6.5],'-','color',[0.6 0.6 0.6])
                        plot([45 45],[0.5 6.5],'-','color',[0.6 0.6 0.6])
                        plot([135 135],[0.5 6.5],'-','color',[0.6 0.6 0.6])
                        set(gca,'ytick',1:6)
                        set(gca,'yticklabel',[100 200 500 1000 2000 4000])
                        axis([0 180 0.5 6.5 ])
                    end
                    
                    
                end
            end
            if plotting==1
                figure(i)
                subplot(x,y,q)
                title(['cluster ',int2str(cluname(q)),' (m=',int2str(NofC(1)),'; ',int2str(100/length(find(tID==1))*NofC(1)),'% / p=',int2str(NofC(2)),'; ',int2str(100/length(find(tID==2))*NofC(2)),'% / h=',int2str(NofC(3)),'; ',int2str(100/length(find(tID==3))*NofC(3)),'%)'])
                
                
                figure(100000+i)
                subplot(x,y,q)
                %                 xEll3 = mean(keepOrig(currID,3))/2 * cos(theta) * cos(mean(keepOrig(currID,5))) - mean(keepOrig(currID,4))/2 * sin(theta) * sin(mean(keepOrig(currID,5))) + mean(keepOrig(currID,1));
                %                 yEll3 = mean(keepOrig(currID,4))/2 * sin(theta) * cos(mean(keepOrig(currID,5))) + mean(keepOrig(currID,3))/2 * cos(theta) * sin(mean(keepOrig(currID,5))) + mean(keepOrig(currID,2));
                %                 plot(xEll3, yEll3, 'LineWidth', 1,'color',colors(q,:));
                if version==2
                    [X,Y]=meshgrid(0.0:0.02:8.5,0.0:0.010:4.100);
                    theta = pi/2-deg2rad(mean(keepOrig(currID,5)));
                    aa = cos(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + sin(theta)^2/2/(mean(keepOrig(currID,4))/1000/1.5)^2;
                    bb = -sin(2*theta)/4/(mean(keepOrig(currID,3))/1.5)^2 + sin(2*theta)/4/(mean(keepOrig(currID,4))/1000/1.5)^2 ;
                    cc = sin(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + cos(theta)^2/2/(mean(keepOrig(currID,4))/1000/1.5)^2;
                    allClust = exp( - (aa*(Y-mean(keepOrig(currID,2))/1000).^2 + 2*bb*(Y-mean(keepOrig(currID,2))/1000).*(X-mean(keepOrig(currID,1))) + cc*(X-mean(keepOrig(currID,1))).^2)) ;
                    [X,Y]=meshgrid(0.0:0.02:8.5,0.0:10:4100);
                else
                    [X,Y]=meshgrid(0:0.02:5,0:0.02:7);
                    theta = -deg2rad(mean(keepOrig(currID,5)));
                    aa = cos(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + sin(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
                    bb = -sin(2*theta)/4/(mean(keepOrig(currID,3))/1.5)^2 + sin(2*theta)/4/(mean(keepOrig(currID,4))/1.5)^2 ;
                    cc = sin(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + cos(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
                    allClust = exp( - (aa*(X-mean(keepOrig(currID,1))).^2 + 2*bb*(X-mean(keepOrig(currID,1))).*(Y-mean(keepOrig(currID,2))) + cc*(Y-mean(keepOrig(currID,2))).^2)) ;
                end
                
                contour(X,Y,allClust(:,:),1,'color',colors(q,:))
                hold on
                plot(mean(keepOrig(currID,1)),mean(keepOrig(currID,2)),'o','markerfacecolor',colors(q,:),'markeredgecolor',colors(q,:),'markersize',2)
                title(['cluster ',int2str(cluname(q)),': av. Gaussian'])
                if version==2
                    axis([0 9 0 4100])
                else
                    axis([-0.5 5.5 0.5 6.5])
                    set(gca,'ytick',1:6)
                    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
                    set(gca,'xtick',1:4)
                    set(gca,'xticklabel',[1 2 4 8])
                end
                
                subplot(x,y,subs)
                contour(X,Y,allClust(:,:),1,'color',colors(q,:))
                hold on
                plot(mean(keepOrig(currID,1)),mean(keepOrig(currID,2)),'o','markerfacecolor',colors(q,:),'markeredgecolor',colors(q,:),'markersize',2)
                title('av. Gaussian')
                if version==2
                    axis([0 9 0 4100])
                else
                    axis([-0.5 5.5 0.5 6.5])
                    set(gca,'ytick',1:6)
                    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
                    set(gca,'xtick',1:4)
                    set(gca,'xticklabel',[1 2 4 8])
                end
                
            end
        else %PC
            
            if plotting==1
                figure(i)
                subplot(x,y,q)
                scatter(keepOrig(currID,1),keepOrig(currID,2), 15,col2(tID(currID),:));
                hold on
                xlabel('PC1')
                ylabel('PC2')
                title(['cluster ',int2str(cluname(q)),' (orig: ',int2str(clusts(q)),')'])
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,2))-0.1 max(keepOrig(:,2))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,1),keepOrig(currID,2), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,2))-0.1 max(keepOrig(:,2))+0.1])
                figure(100+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,1),keepOrig(currID,3), 15,col2(tID(currID),:));
                hold on
                xlabel('PC1')
                ylabel('PC3')
                title(['cluster ',int2str(cluname(q)),' (orig: ',int2str(clusts(q)),')'])
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,3))-0.1 max(keepOrig(:,3))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,1),keepOrig(currID,3), 15,col2(tID(currID),:));
                hold on
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,3))-0.1 max(keepOrig(:,3))+0.1])
                figure(1000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,1),keepOrig(currID,4), 15,col2(tID(currID),:));
                hold on
                xlabel('PC1')
                ylabel('PC4')
                title(['cluster ',int2str(cluname(q)),' (orig: ',int2str(clusts(q)),')'])
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,4))-0.1 max(keepOrig(:,4))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,1),keepOrig(currID,4), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,4))-0.1 max(keepOrig(:,4))+0.1])
                figure(10000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,1),keepOrig(currID,5), 15,col2(tID(currID),:));
                hold on
                xlabel('PC1')
                ylabel('PC5')
               title(['cluster ',int2str(cluname(q)),' (orig: ',int2str(clusts(q)),')'])
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,5))-0.1 max(keepOrig(:,5))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,1),keepOrig(currID,5), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,1))-0.1 max(keepOrig(:,1))+0.1 min(keepOrig(:,5))-0.1 max(keepOrig(:,5))+0.1])
                figure(100000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,2),keepOrig(currID,3), 15,col2(tID(currID),:));
                hold on
                xlabel('PC2')
                ylabel('PC3')
                title(['cluster ',int2str(cluname(q)),' (orig: ',int2str(clusts(q)),')'])
                axis([min(keepOrig(:,2))-0.1 max(keepOrig(:,2))+0.1 min(keepOrig(:,3))-0.1 max(keepOrig(:,3))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,2),keepOrig(currID,3), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,2))-0.1 max(keepOrig(:,2))+0.1 min(keepOrig(:,3))-0.1 max(keepOrig(:,3))+0.1])
                figure(1000000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,4),keepOrig(currID,5), 15,col2(tID(currID),:));
                hold on
                xlabel('PC4')
                ylabel('PC5')
                title(['cluster ',int2str(cluname(q)),' (orig: ',int2str(clusts(q)),')'])
                axis([min(keepOrig(:,4))-0.1 max(keepOrig(:,4))+0.1 min(keepOrig(:,5))-0.1 max(keepOrig(:,5))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,4),keepOrig(currID,5), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,4))-0.1 max(keepOrig(:,4))+0.1 min(keepOrig(:,5))-0.1 max(keepOrig(:,5))+0.1])
                if size(keepOrig,2)>5
                     figure(10000000+i)
                subplot(x,y,q)
                scatter(keepOrig(currID,6),keepOrig(currID,7), 15,col2(tID(currID),:));
                hold on
                xlabel('CoM (x-coord)')
                ylabel('CoM (y-coord)')
                title(['cluster ',int2str(cluname(q)),' (orig: ',int2str(clusts(q)),')'])
                axis([min(keepOrig(:,6))-0.1 max(keepOrig(:,6))+0.1 min(keepOrig(:,7))-0.1 max(keepOrig(:,7))+0.1])
                subplot(x,y,subs)
                scatter(keepOrig(currID,6),keepOrig(currID,7), 15,colors(q,:));
                hold on
                axis([min(keepOrig(:,6))-0.1 max(keepOrig(:,6))+0.1 min(keepOrig(:,7))-0.1 max(keepOrig(:,7))+0.1])
                end
            end
        end
        
    end
    if plotting==1 && gau==1
        figure(i)
        subplot(x,y,subs)
        title('overview')
        if version==2
            axis([0 9 0 4100])
        else
            axis([0.5 4.5 0.5 6.5])
            set(gca,'ytick',1:6)
            set(gca,'yticklabel',[100 200 500 1000 2000 4000])
            set(gca,'xtick',1:4)
            set(gca,'xticklabel',[1 2 4 8])
        end
        
    end

%% statistics of clusters
clear
close all
path1='F:\all_species_temp\newAttempt\DG_Gaussian';
type='GaussianNew2_corr';
% type='PC_F1';
% clustM='clust_Hierarch_';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))

uclu=unique(data);
stats=zeros(length(uclu),10);



for u=1:length(uclu)
    curr=find(data==uclu(u));
    curr=keepOrig(curr,:);
    
    if u==6 || u==10
        tmp=find(curr(:,5)<90);
        curr(tmp,5)=180-curr(tmp,5);
    end
    me=mean(curr);
    stats(u,1:5)=me;
    st=std(curr);
    stats(u,6:10)=st;
end
%% similarity matrices
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='F:\all_species_temp\newAttempt\DG_PCanalysis';
% type='Gaussian';
type='PC_F1new';
% clustM='clust_Hierarch_';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))

sorted=ID;
% sorted=Tmat;



I=orig2clust; %contains parameters that went into clustering (?)


D = squareform( pdist(I) );       %'# euclidean distance between rows of I
M = exp(-(1/10) * D);              %# similarity matrix between rows of I
figure
imagesc(M)
title(['Similarity matrix non-clustered data'])

% for each parameter
for q=1:size(I,2)
    II=I(:,q);
    DD = squareform( pdist(II) );       %'# euclidean distance between rows of I
    MM = exp(-(1/10) * DD);
    MM2=linkage(MM);
    figure
    [t1 t2 t3]=dendrogram(MM2,size(MM,1));
    MM3=MM(t3,t3);
    imagesc(MM3)
end


figure
M2=linkage(M);
[t1 t2 t3]=dendrogram(M2,size(M,1));
% close
M3=M(t3,t3);
imagesc(M3)
title(['Similarity matrix sorted non-clustered data'])


for clust=1:size(sorted,2)
    data=sorted(:,clust);
    [a idx1]=sort(data);
    idx=repmat(idx1,1,size(orig2clust,2));
    coll3=orig2clust(idx);
    I=coll3;
    D = squareform( pdist(I) );       %'# euclidean distance between columns of I
    M = exp(-(1/10) * D);              %# similarity matrix between columns of I
    figure
    imagesc(M)
    changes=1;
    for u=2:length(a)
        if a(u)~=a(u-1)
            changes=[changes u];
        end
    end
    set(gca,'xtick',changes)
    set(gca,'xticklabel',unique(a))
    set(gca,'ytick',changes)
    set(gca,'yticklabel',unique(a))
end
%% !!new silhouette attempt
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='e:\all_species_temp\newAttemptOnlyM\DG_PCanalysis'; %newAttempt3: corrected for outlier
% type='GaussianNew2';
type='PC_F1';
% clustM='clust_MixGau_';
clustM='clust_Kmeans_';
% clustM='clust_fuzzyGK_';
load(fullfile(path1,[clustM,type]))

% sorted=data;
sorted=ID;
% sorted=Tmat;


close all
counting=0; maxval=0; colors={'b','k'};
clear collBSS

sis=[]; ai=[]; bi=[];
for clust=1:size(sorted,2)
    counting=counting+1;
    data=sorted(:,clust);
    
    if size(sorted,2)==1
        uclu=unique(data);
        for u=1:length(uclu)
            curr=find(data==uclu(u));
            curr2=orig2clust(curr,:);
            
            if u==6 || u==10
                tmp=find(curr2(:,5)<1.5);
                curr2(tmp,5)=pi-curr2(tmp,5);
                orig2clust(curr,5)=curr2(:,5);
            end
            
        end
    end
    

   
   
  ais=[]; bis=[];
    for c=1:size(orig2clust,1)
        cpoint=orig2clust(c,:);
        others=1:size(orig2clust,1);
        others=setdiff(others,c)';
        osame=find(data(others)==data(c));
        osame=others(osame);
        odiff=find(data(others)~=data(c));
        odiff=others(odiff);
        
        
        disame=sqrt(sum((repmat(cpoint,length(osame),1)-orig2clust(osame,:)).^2));
        disame=mean(disame);
        didiff=sqrt(sum((repmat(cpoint,length(odiff),1)-orig2clust(odiff,:)).^2));
        didiff=min(didiff);
        
        ais=[ais;disame];
        bis=[bis;didiff];
    end
    
    ai=[ai ais]; bi=[bi bis];
    si=(bis-ais)./max([ais bis],[],2);
   sis=[sis si];

end
sis=round(sis*100)/100;

figure
for a=1:size(ai,2)
    plot(a+1-0.1,mean(ai(:,a)),'ob','markersize',10,'markerfacecolor','b')
    hold on
    plot(a+1+0.1,mean(bi(:,a)),'or','markersize',10,'markerfacecolor','r')
    if a==1
        legend('mean eucl. dist within','min eucl. dist to others')
    end
    plot(a+1-0.1,ai(:,a),'ob','markersize',4) % mean eucl. dist. within clust
    hold on
    plot(a+1+0.1,bi(:,a),'or','markersize',4) % min eucl. dist to others
    
end

figure
if size(sorted,2)==1
    hist(sis)
    title([int2str(length(unique(sorted))),' clusters'])
    hold on
    plot([0 0],[0 200])
    axis([-inf inf 0 200])
else
    for s=1:size(sis,2)
        subplot(5,6,s)
        hist(sis(:,s))
        title([int2str(s+1),' clusters'])
        hold on
        plot([0 0],[0 200])
        axis([-inf inf 0 200])
    end
end
figure
col=colormap('jet'); col=col(1:2:end,:); ra=randperm(32); col=col(ra,:);
if size(sorted,2)==1
    cl=1;
    [clu,idx]=sort(sorted(:,:));
    currsis=sis(idx,:);
    for p=1:size(sis,1)
        plot([p p],[0 currsis(p)],'-','color',col(clu(p),:))
        hold on
        cl2=clu(p);
        if cl~=cl2
            disp([int2str(cl2),'/',int2str(p)])
        end
        cl=cl2;
    end
    axis([0 size(sis,1) -1 1 ])
    title([int2str(length(unique(sorted))),' clusters'])
else
    for s=1:size(sis,2)
        subplot(5,6,s)
        [clu,idx]=sort(sorted(:,s));
        currsis=sis(idx,s);
        for p=1:size(sis,1)
            plot([p p],[0 currsis(p)],'-','color',col(clu(p),:))
            hold on
        end
        axis([0 size(sis,1) -1 1 ])
        title([int2str(s+1),' clusters'])
    end
end
figure
plot(2:size(sis,2)+1,mean(sis,1),'o')
hold on
for s=1:size(sis,2)
    plot([s+1 s+1],[min(sis(:,s)) max(sis(:,s))],'-')
end
plot([1 31],[0.5 0.5],'-k')
plot([1 31],[0.75 0.75],'-k')
plot([1 31],[0 0],'-k')
title('mean silhouette +/- max and min')
axis([1 31 -1 1])
%% new silhouette attempt 2
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='e:\all_species_temp\newAttempt3\DG_PCanalysis';
% type='GaussianNew2';
type='PC_F1new';
% clustM='clust_MixGau_';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))

% sorted=data;
sorted=ID;
% sorted=Tmat;


close all
counting=0; maxval=0; colors={'b','k'};
clear collBSS

sis=[]; ai=[]; bi=[];
for clust=1:size(sorted,2)
    counting=counting+1;
    data=sorted(:,clust);
    
  
    
    uclu=unique(data);
   
     allsis=[];
    for c=1:size(orig2clust,1)
        cpoint=orig2clust(c,:);
        bi=[]; ais=[]; bis=[];
        for uc=1:length(uclu)
            datapoints=find(data==uc);
            if uc==data(c)
                datapoints=setdiff(datapoints,c);
            end
            dist=mean(sum((repmat(cpoint,length(datapoints),1)-orig2clust(datapoints,:)).^2,2)); %mean distance to all points in cluster
            
            if uc==data(c)
               ais=dist;
            else
                bi=[bi;dist];
            end
        end
        bis=min(bi);
  
%     for c=1:size(orig2clust,1)
%         cpoint=orig2clust(c,:);
%         others=1:size(orig2clust,1);
%         others=setdiff(others,c)';
%         osame=find(data(others)==data(c));
%         osame=others(osame);
%         odiff=find(data(others)~=data(c));
%         odiff=others(odiff);
%         
%         
%         disame=sqrt(sum((repmat(cpoint,length(osame),1)-orig2clust(osame,:)).^2));
%         disame=mean(disame);
%         didiff=sqrt(sum((repmat(cpoint,length(odiff),1)-orig2clust(odiff,:)).^2));
%         didiff=min(didiff);
%         
%         ais=[ais;disame];
%         bis=[bis;didiff];
%     end
    
%     ai=[ai ais]; bi=[bi bis];
    si=(bis-ais)./max([ais bis],[],2);
    allsis=[allsis; si];
    end
    sis=[sis allsis];
end
sis=round(sis*100)/100;

figure
for a=1:size(ai,2)
    plot(a+1-0.1,mean(ai(:,a)),'ob','markersize',10,'markerfacecolor','b')
    hold on
    plot(a+1+0.1,mean(bi(:,a)),'or','markersize',10,'markerfacecolor','r')
    if a==1
        legend('mean eucl. dist within','min eucl. dist to others')
    end
    plot(a+1-0.1,ai(:,a),'ob','markersize',4) % mean eucl. dist. within clust
    hold on
    plot(a+1+0.1,bi(:,a),'or','markersize',4) % min eucl. dist to others
    
end

figure
if size(sorted,2)==1
    hist(sis)
    title([int2str(length(unique(sorted))),' clusters'])
    hold on
    plot([0 0],[0 200])
    axis([-inf inf 0 200])
else
    for s=1:size(sis,2)
        subplot(5,6,s)
        hist(sis(:,s))
        title([int2str(s+1),' clusters'])
        hold on
        plot([0 0],[0 200])
        axis([-inf inf 0 200])
    end
end
figure
col=colormap('jet'); col=col(1:2:end,:); ra=randperm(32); col=col(ra,:);
if size(sorted,2)==1
    cl=1;
    [clu,idx]=sort(sorted(:,:));
    currsis=sis(idx,:);
    for p=1:size(sis,1)
        plot([p p],[0 currsis(p)],'-','color',col(clu(p),:))
        hold on
        cl2=clu(p);
        if cl~=cl2
            disp([int2str(cl2),'/',int2str(p)])
        end
        cl=cl2;
    end
    axis([0 size(sis,1) -1 1 ])
    title([int2str(length(unique(sorted))),' clusters'])
else
    for s=1:size(sis,2)
        subplot(5,6,s)
        [clu,idx]=sort(sorted(:,s));
        currsis=sis(idx,s);
        for p=1:size(sis,1)
            plot([p p],[0 currsis(p)],'-','color',col(clu(p),:))
            hold on
        end
        axis([0 size(sis,1) -1 1 ])
        title([int2str(s+1),' clusters'])
    end
end
figure
plot(2:size(sis,2)+1,mean(sis,1),'o')
hold on
for s=1:size(sis,2)
    plot([s+1 s+1],[min(sis(:,s)) max(sis(:,s))],'-')
end
plot([1 31],[0.5 0.5],'-k')
plot([1 31],[0.75 0.75],'-k')
plot([1 31],[0 0],'-k')
title('mean silhouette +/- max and min')
axis([1 31 -1 1])
%% BSS and WSS
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='F:\all_species_temp\newAttempt2\DG_PCanalysis';
% type='Gaussian';
type='PC_F1';
% clustM='clust_Hierarch_';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))

sorted=ID;
% sorted=Tmat;


close all
counting=0; maxval=0; colors={'b','k'};
clear collBSS
figure(1)
for clust=1:size(sorted,2)
    counting=counting+1;
    figure(1)
    subplot(6,6,counting)
    data=sorted(:,clust);
    [a idx1]=sort(data);
    changes=1;
    for u=2:length(a)
        if a(u)~=a(u-1)
            
            changes=[changes u];
        end
    end
    changes=[changes size(orig2clust,1)+1];
    BSS=0; overallmean=mean(orig2clust); maxval=0; collWSS=0;
    for j=1:length(changes)-1
        currdata=orig2clust(idx1(changes(j):changes(j+1)-1),:);
        currmean=mean(currdata);
        WSS=0;
        for s=1:size(currdata,1)
            currcurr=currdata(s,:);
            WSS=WSS+(currcurr-currmean).^2;
        end
        collWSS=collWSS+sum(WSS);
        maxval=max(maxval,max(WSS));
        figure(1)
        subplot(6,6,counting)
        plot(repmat(j,length(WSS),1),WSS,'ko-','markersize',1.5)
        hold on
        BSS=BSS+size(currdata,1)*(overallmean-currmean).^2;
    end
    collBSS(counting,1:size(BSS,2))=BSS; %between clusters sum of squares (shoud be bit)
    axis([0 max(a)+1 0 maxval+0.1])
    
    title([int2str(max(data)),' clusters / ',num2str(collWSS)]) %witin cluster sum of squares (shoud be small)
    figure(10)
    plot(max(a), BSS,'o')
    hold on
    
    
    %silhouette coefficient
    figure(3)
    subplot(6,6,counting)
    silhouette(orig2clust(idx1,:),a)
    sil=silhouette(orig2clust(idx1,:),a);%high value = well matched to own cluster, poorly matched to others
    
    %        figure(clust+100)
    %        silhouette(coll2(idx1),a);
    figure(2)
    collsil=0; clust2take=0;
    subplot(6,6,counting)
    for b=1:length(changes)-1
        medsil=mean(sil(changes(b):changes(b+1)-1));
        if ~isnan(medsil)
            colors=circshift(colors,[0,-1]);
            clust2take=clust2take+1;
            plot(changes(b):changes(b+1)-1,sil(changes(b):changes(b+1)-1),'.','color',colors{1})
            hold on
            collsil=collsil+abs(medsil);
            plot(changes(b)+(changes(b+1)-changes(b))/2,medsil,'or','markerfacecolor','r','markersize',3)
        end
    end
    
    title([int2str(max(data)),' clust/ av. sil.: ',num2str(collsil/clust2take)])
    
    figure(4)
    plot(max(data),mean(sil),'o')
    hold on
    
    
end
%% hierarchical clustering
clear
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% type='Gaussian';
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
type='PC_F1new';

load(fullfile(path1,type))
savename='clust_Hierarch_';
close all

orig=pcs;
gau=0;
% orig=Gauss;
% gau=1;

if gau==1
    orig2clust=[]; tID=zeros(sum(datalength),1);
    last=1;
    for i=1:3
        curr=reshape(Gauss(i,:,:),300,5); %centerX, centerY, x-axis, y-axis, angle
        curr=curr(1:datalength(i),:);
        orig2clust=[orig2clust;curr];
        tID(last:last+datalength(i)-1)=repmat(i,datalength(i),1);
        last=last+datalength(i);
    end
    tmp=find(orig2clust(:,5)<0);
    tmp2=abs(-deg2rad(180)-orig2clust(tmp,5));
    orig2clust(tmp,5)=tmp2;
else
    orig2clust=orig;
end

Y=pdist(orig2clust);
tst=squareform(Y);
Z = linkage(Y,'centroid');

Tmat=[];
for cut=2:35
    T=cluster(Z,'maxclust',cut);
    Tmat=[Tmat T];
    
end

%cut=14 seems to be optimal (and everything above)
figure
dendrogram(Z)

if gau==1
    keepOrig=orig2clust;
    ang=rad2deg(orig2clust(:,5));
    orig2clust(:,5)=ang;
else
    keepOrig=orig2clust;
end
ID=Tmat;
save(fullfile(path1,[savename,type]),'orig2clust','keepOrig','ID','tID')
%% F2 of ONOFF cells in human
load(fullfile('F:\all_species_temp\newAttempt\DG','h_normPeaks'))
load(fullfile('F:\all_species_temp\h_info'))

good=find(goodones==1);
fullfield=cell2mat(info(good,8));
full=find(fullfield==2);
onoff=good(full);
full=find(fullfield==1);
on=good(full);

onoffID=[];
for o=1:length(onoff)
    tmp=find(togo==onoff(o));
    onoffID=[onoffID;tmp];
end

onID=[];
for o=1:length(on)
    tmp=find(togo==on(o));
    onID=[onID;tmp];
end

f2ONOFF=normPeaks(onoffID,:,3);
f2ON=normPeaks(onID,:,3);

f2meanONOFF=mean(f2ONOFF)';
f2meanON=mean(f2ON)';

compMean=[f2meanONOFF f2meanON];
diMean=diff(fliplr(compMean)')';
spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
compMean=[compMean diMean spat' freq'];

f2meanONOFF2=[]; Noff=[];
for o=1:size(f2ONOFF,2)
    tmp=find(f2ONOFF(:,o)>0);
    f2meanONOFF2=[f2meanONOFF2;mean(f2ONOFF(tmp,o))];
    Noff=[Noff;100/size(f2ONOFF,1)*length(tmp)];
end
f2meanON2=[];Non=[];
for o=1:size(f2ON,2)
    tmp=find(f2ON(:,o)>0);
    f2meanON2=[f2meanON2;mean(f2ON(tmp,o))];
    Non=[Non;100/size(f2ON,1)*length(tmp)];
end


compMean2=[f2meanONOFF2 f2meanON2];
diMean2=diff(fliplr(compMean2)')';
spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
compMean2=[compMean2 diMean2 spat' freq' Noff Non];
compMean=[compMean Noff Non];

load(fullfile('F:\all_species_temp','newAttempt','DG','h_DGFFT_binary'))
keep=FFTabs;
load(fullfile('F:\all_species_temp','newAttempt','DG','h_DGFFT_freq'))

totake=15;
figure
cnt=0;
for o=1:length(onoffID)
    
    cnt=cnt+1;
    subplot(4,4,cnt)
    plot(FFTabs{onoffID(o),6,totake}(1:181),keep{onoffID(o),6,totake}(1:181))
end


figure
cnt=0;o=0;
while cnt<16
    o=o+1;
    cnt=cnt+1;
    subplot(4,4,cnt)
    plot(FFTabs{onID(o),6,totake}(1:181),keep{onID(o),6,totake}(1:181))
end

%% ... Mixture of Gaussians
orig=Gauss;
gau=1;


if gau==1
    
    orig2clust=[]; tID=zeros(sum(datalength),1);
    last=1;
    for i=1:3
        curr=reshape(Gauss(i,:,:),300,5); %centerX, centerY, x-axis, y-axis, angle
        curr=curr(1:datalength(i),:);
        orig2clust=[orig2clust;curr];
        tID(last:last+datalength(i)-1)=repmat(i,datalength(i),1);
        last=last+datalength(i);
    end
    tmp=find(orig2clust(:,5)<0);
    tmp2=abs(-deg2rad(180)-orig2clust(tmp,5));
    orig2clust(tmp,5)=tmp2;
    keepOrig=orig2clust;
else
    orig2clust=orig;
end
for p=1:size(orig2clust,2)
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
end

for k=2:30
    GMModel = fitgmdist(orig2clust,k);
end


figure
y = [zeros(size(orig2clust,1),1)];
h = gscatter(orig2clust(:,1),orig2clust(:,2),y);
hold on
ezcontour(@(x1,x2,x3,x4,x5)pdf(GMModel,[x1 x2 x3 x4 x5]),get(gca,{'XLim','YLim'}))


%% ...just for testing
load fisheriris
classes = unique(species)
[~,score] = pca(meas,'NumComponents',2); %data for clustering
GMModels = cell(3,1); % Preallocation
options = statset('MaxIter',1000);
rng(1); % For reproducibility

for j = 1:3
    GMModels{j} = fitgmdist(score,j,'Options',options);
    fprintf('\n GM Mean for %i Component(s)\n',j)
    Mu = GMModels{j}.mu
end

figure
for j = 1:3
    subplot(2,2,j)
    gscatter(score(:,1),score(:,2),species)
    h = gca;
    hold on
    ezcontour(@(x1,x2)pdf(GMModels{j},[x1 x2]),...
        [h.XLim h.YLim],100)
    title(sprintf('GM Model - %i Component(s)',j));
    xlabel('1st principal component');
    ylabel('2nd principal component');
    if(j ~= 3)
        legend off;
    end
    hold off
end
g = legend;
g.Position = [0.7 0.25 0.1 0.1];
%% !!clustering mixture of gaussians
clear
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% type='GaussianNew2';
% path1='F:\all_species_temp\newAttempt\DG_PCanalysis';
% type='PC_F2';

path1='F:\all_species_temp\newAttempt4\DG_PCanalysis';
type='PC_F1new';

load(fullfile(path1,type))
savename='clust_MixGau_';
close all
outlier=1;
% orig=amps;
gau=0;
orig=pcs;
% gau=0;
% orig=Gauss;
% gau=1;


if gau==1
    
    orig2clust=[]; tID=zeros(sum(datalength),1);
    last=1;
    for i=1:3
        curr=reshape(Gauss(i,:,:),300,5); %centerX, centerY, x-axis, y-axis, angle
        curr=curr(1:datalength(i),:);
        orig2clust=[orig2clust;curr];
        tID(last:last+datalength(i)-1)=repmat(i,datalength(i),1);
        last=last+datalength(i);
    end
    tmp=find(orig2clust(:,5)<0);
    tmp2=abs(-deg2rad(180)-orig2clust(tmp,5));
    orig2clust(tmp,5)=tmp2;
    keepOrig=orig2clust;
else
    orig2clust=orig;
end
if outlier==1
    for p=1:size(orig2clust,2)-1
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
    end
curr=orig2clust(:,end);
[tmp tmp2]=sort(curr);
norm=(curr-min(curr))/(tmp(end-1)+0.2-min(curr));
norm(tmp2(end))=1;
orig2clust(:,end)=norm;
else
for p=1:size(orig2clust,2)
    curr=orig2clust(:,p);
    norm=(curr-min(curr))/(max(curr)-min(curr));
    orig2clust(:,p)=norm;
end
end

data=orig2clust;
mn = min(data); mx = max(data);

ID=[]; pvals=cell(29,1); gmms=cell(29,1);
for K=2:30
    % fit a GMM model
    gmm = fitgmdist(data, K, 'Options',statset('MaxIter',1000),'covtype','diagonal','replicates',1000);
    [clustInd,~,p] = cluster(gmm, data);
    ID=[ID clustInd];
    pvals{K-1}=p;
    gmmx{K-1}=gmm;
    K
end

if gau==1
    ang=rad2deg(keepOrig(:,5));
    keepOrig(:,5)=ang;
else
    keepOrig=orig2clust;
end
save(fullfile(path1,[savename,type]),'orig2clust','keepOrig','ID','tID','gmmx')

%% plotting of mix gauss
d1=1; d2=2;

for kk=1:length(gmmx)
    gmm=gmmx{kk};
    K=gmm.NComponents;
    allClust=zeros(101,101,K);
    
    for c=1:K
        xcnt=0;
        for x=0.00:0.01:1
            ycnt=0; xcnt=xcnt+1;
            for y=0.00:0.01:1
                ycnt=ycnt+1;
                allClust(xcnt,ycnt,c)=gmm.PComponents(c)*exp(-((x-gmm.mu(c,d1))^2/(2*gmm.Sigma(1,d1,c))+(y-gmm.mu(c,d2))^2/(2*gmm.Sigma(1,d2,c))));
                allClust(xcnt,ycnt,c)=exp(-((x-gmm.mu(c,d1))^2/(2*gmm.Sigma(1,d1,c))+(y-gmm.mu(c,d2))^2/(2*gmm.Sigma(1,d2,c))));
                
            end
        end
    end
    
    figure
    for c=1:K
        subplot(4,4,c)
        imagesc(allClust(:,:,c))
    end
    subplot(4,4,16)
    imagesc(mean(allClust,3))
    
    
    figure
    [X,Y]=meshgrid(0:0.01:1,0:0.01:1);
    for c=1:K
        contour(X,Y,allClust(:,:,c),3)
        hold on
    end
end
%% !!explained variance
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='F:\all_species_temp\newAttempt6\DG_PCanalysis';
% type='GaussianNew2';
type='means_F1';
clustM='clust_fuzzyGK_';
% clustM='clust_Cmeans_';
load(fullfile(path1,[clustM,type]))

% sorted=data;
sorted=ID;
% sorted=Tmat;

grandMean=mean(orig2clust);


counting=0; between=[];within=[];
totVars=[];
for clust=1:size(sorted,2)
    counting=counting+1;
    data=sorted(:,clust);
    
    if size(sorted,2)==1
        uclu=unique(data);
        for u=1:length(uclu)
            curr=find(data==uclu(u));
            curr2=orig2clust(curr,:);
            
            if u==6 || u==10
                tmp=find(curr2(:,5)<1.5);
                curr2(tmp,5)=pi-curr2(tmp,5);
                orig2clust(curr,5)=curr2(:,5);
            end
            
        end
    end
    uclust=unique(data);
    
    totVar=[];
    for c=1:length(uclust)
        now=find(data==c);
        curr=sum(orig2clust(now,:)-repmat(grandMean,length(now),1).^2,1);
        totVar=[totVar;curr];
    end
    totVar= sum(totVar)/(size(data,1)-1);
    totVars=[totVars;totVar];
    
    vars=[];
    for c=1:length(uclust)
        now=find(data==c);
        curr=length(now).*(mean(orig2clust(now,:))-grandMean).^2;
        vars=[vars;curr];
    end
    vv=[];
    for v=1:size(vars,2)
        tmp=find(~isnan(vars(:,v)));
        va=sum(vars(tmp,v))/(length(uclust)-1);
        vv=[vv va];
    end
    vars=vv;
    between=[between;vars];
    
    vars=[];
    for c=1:length(uclust)
        now=find(data==c);
        curr=sum((orig2clust(now,:)-repmat(mean(orig2clust(now,:)),length(now),1)).^2,1);
        vars=[vars;curr];
    end
      vv=[];
    for v=1:size(vars,2)
        tmp=find(~isnan(vars(:,v)));
        va=sum(vars(tmp,v))/(size(data,1)-length(uclust));
        vv=[vv va];
    end
    vars=vv;
   
    within=[within;vars];
end


figure
    plot(2:30,sum(between')./sum(within'),'o-')
    title('F-test var(between)/var(within)')

figure
subplot(1,2,1)
plot(2:30,sum(within'),'o-')
title('within variance')
subplot(1,2,2)
plot(2:30,sum(between'),'or-')
title('between variance')
%% !!fuzzy validity indices
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='F:\all_species_temp\newAttempt6\DG_PCanalysis';
% type='GaussianNew2';
type='means_F1';
clustM='clust_fuzzyGK_';
% clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))

% sorted=data;
sorted=ID;
% sorted=Tmat;
tilist={'partition index','separation index','Xie and Beni index','Dunn index','alternative Dunn index'};
figure
for s=1:5
    subplot(2,3,s)
    plot(2:30,valids(s,:),'o-')
    title(tilist{s})
end
%% !!only once!!!! make clustering numbers from 1 to x (fuzzy GK)
clear
close all
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';

load(fullfile(path1,'clust_fuzzyGK_PC_F1new'))

oclust=cell(29,1);
for i=1:size(ID,2)
    un=unique(ID(:,i));
    oclust{i}=un;
    new=1:length(un);
    data=ID(:,i);
    data2=ID(:,i);
    for u=1:length(new)
        tmp=find(data==un(u));
        data2(tmp)=u;
    end
    ID(:,i)=data2;
end
save(fullfile(path1,'clust_fuzzyGK_PC_F1new'),'keepOrig','ID','tID','orig2clust','valids','oclust')
%% !!check BIC for MixGau
figure
for g=1:length(gmmx)
    curr=gmmx{g};
    bic=curr.BIC;
    aic=curr.AIC;
    plot(g+1,bic,'o')
    hold on
    plot(g+1,aic,'or')
    if g==1
        legend('bic','aic')
    end
end
%% compare mixGau and Kmeans
clear
close all
path1='F:\all_species_temp\newAttempt\DG_PCanalysis';
type='PC_F1new';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))
sorted=ID;

clustM='clust_MixGau_';
load(fullfile(path1,[clustM,type]))
sorted2=ID;


clu=13;
clu2=15;
data=sorted(:,clu);
data2=sorted2(:,clu2);

figure
for c=1:clu+1
    tst=find(data==c);
    tst2=(data2(tst));
    u2=unique(tst2);
    nums=[];
    for u=1:length(u2)
        tmp=find(tst2==u2(u));
        nums=[nums;length(tmp)];
        plot(c,u2(u),'o')
        hold on
    end
    [~,id]=max(nums);
    plot(c,u2(id),'dr','markersize',10,'markerfacecolor','r')
end
%% adjust clustering result: Kmeans, 19 clust
clear
close all
path1='F:\all_species_temp\newAttempt\DG_Gaussian';
type='Gaussian';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))
sorted=ID;

clu=18;
data=sorted(:,clu);

tmp=find(data==13);
data(tmp)=5;
tmp=find(data==17);
data(tmp)=8;

unclust=unique(data);
unclust=setdiff(unclust,14);

% newclust=[];
% todist=find(data==14);
% for t=1:length(todist)
%     meandists=[];
%     for u=1:length(unclust)
%         tmp=find(data==unclust(u));
%         dist=sqrt(sum((repmat(orig2clust(todist(t),:),length(tmp),1)-orig2clust(tmp,:)).^2));
%         sudist=sum(dist);
%         meandists=[meandists;dist];
%     end
%     meandists=sum(meandists');
%     [~,id]=min(meandists);
%     newclust=[newclust;unclust(id)];
% %     mis=[];
% %     for m=1:5
% %         [~,val]=sort(meandists(:,m));
% %         mis=[mis val];
% %     end
% %     mis=mis(1:3,:);
% %     unew=unique(mis);
% %     nu=[];
% %     for un=1:length(unew);
% %         tst=find(mis==unew(un));
% %         nu=[nu;length(tst)];
% %     end
% %     [tmp,id]=max(nu);
% %     totake=unclust(unew(id));
% %     newclust=[newclust;totake];
% end
% data(todist)=newclust;

un=unique(data);
data2=data;

for n=1:length(un)
    tmp=find(data==un(n));
    data2(tmp)=n;
end
data=data2;

save(fullfile('F:\all_species_temp\newAttempt\DG_Gaussian',[clustM,type,'_corr']),'data','orig2clust','keepOrig','tID')

%% adjust clustering result: Kmeans, 16 clust
clear
close all
path1='F:\all_species_temp\newAttempt\DG_Gaussian';
type='GaussianNew2';
clustM='clust_Kmeans_';
load(fullfile(path1,[clustM,type]))
sorted=ID;

clu=15;
data=sorted(:,clu);

tmp=find(data==9);
data(tmp)=6;
tmp=find(data==13);
data(tmp)=11;
tmp=find(data==15);
data(tmp)=12;




un=unique(data);
data2=data;

for n=1:length(un)
    tmp=find(data==un(n));
    data2(tmp)=n;
end
data=data2;

save(fullfile('F:\all_species_temp\newAttempt\DG_Gaussian',[clustM,type,'_corr']),'data','orig2clust','keepOrig','tID')
%% plot mean heatmaps for cluster result
clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
path1='e:\all_species_temp\newAttempt3\DG_PCanalysis';
file='clust_Kmeans_PC_F1';
load(fullfile(path1,file))
sel=1; clust=12; %actual clust
gau=0;

if sel==1
    data=ID(:,clust-1);
end

superpath='e:\all_species_temp';
Dtype='wph';
normDG=[];
for ty=3
    load(fullfile(superpath,'newAttempt2','DG',[Dtype(ty),'_normPeaks']))
    normDG=[normDG;normPeaks];
end

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0.2 0.2 0.2; 0 1 1; 1 1 0;0.5 0 0];

clu=unique(data);
for c=1:length(clu)
    currID=find(data==clu(c));
    if gau==1
    figure(1)
    
    for s=1:length(clu)
        subplot(5,6,s)
        [X,Y]=meshgrid(0:0.02:5,0:0.02:7);
        theta = -deg2rad(mean(keepOrig(currID,5)));
        aa = cos(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + sin(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
        bb = -sin(2*theta)/4/(mean(keepOrig(currID,3))/1.5)^2 + sin(2*theta)/4/(mean(keepOrig(currID,4))/1.5)^2 ;
        cc = sin(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + cos(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
        allClust = exp( - (aa*(X-mean(keepOrig(currID,1))).^2 + 2*bb*(X-mean(keepOrig(currID,1))).*(Y-mean(keepOrig(currID,2))) + cc*(Y-mean(keepOrig(currID,2))).^2)) ;
        
        contour(X,Y,allClust(:,:),1,'color',[0.9 0.9 0.9])
        hold on
        
        
    end
    
    
   subplot(4,5,c)
    [X,Y]=meshgrid(0:0.02:5,0:0.02:7);
    theta = -deg2rad(mean(keepOrig(currID,5)));
    aa = cos(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + sin(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
    bb = -sin(2*theta)/4/(mean(keepOrig(currID,3))/1.5)^2 + sin(2*theta)/4/(mean(keepOrig(currID,4))/1.5)^2 ;
    cc = sin(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + cos(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
    allClust = exp( - (aa*(X-mean(keepOrig(currID,1))).^2 + 2*bb*(X-mean(keepOrig(currID,1))).*(Y-mean(keepOrig(currID,2))) + cc*(Y-mean(keepOrig(currID,2))).^2)) ;
    
    contour(X,Y,allClust(:,:),1,'color',colors(c,:))
    hold on
    plot(mean(keepOrig(currID,1)),mean(keepOrig(currID,2)),'o','markerfacecolor',colors(c,:),'markeredgecolor',colors(c,:),'markersize',2)
    title(['cluster ',int2str(clu(c)),': av. Gaussian'])
    
    axis([-0.5 5.5 0.5 6.5])
    set(gca,'ytick',1:6)
    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',[1 2 4 8])
    
    subplot(4,5,length(clu)+1)
    contour(X,Y,allClust(:,:),1,'color',colors(c,:))
    hold on
    axis([-0.5 5.5 0.5 6.5])
    set(gca,'ytick',1:6)
    set(gca,'yticklabel',[100 200 500 1000 2000 4000])
    set(gca,'xtick',1:4)
    set(gca,'xticklabel',[1 2 4 8])
    end
    
    figure(2)
    subplot(4,5,c)
    curr=normDG(currID,:,2);
    curr=median(curr,1);
    curr=fliplr(reshape(curr,4,6));
    [y,x] = meshgrid(1:6,1:4);
    imagesc(x(:)',y(:)',curr');
    colormap('gray')
    title([int2str(c),' (mean; ',int2str(length(currID)),'cells)'])
     figure(3)
    subplot(4,5,c)
    curr=normDG(currID,:,2);
    curr=std(curr,0,1);
    curr=fliplr(reshape(curr,4,6));
    [y,x] = meshgrid(1:6,1:4);
    imagesc(curr'*100,[0 100]);
    colormap('gray')
       title([int2str(c),' (std)'])
end
%% !!plot all heatmaps per cluster
clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
% file='clust_MixGau_GaussianNew2';
path1='e:\all_species_temp\newAttempt3\DG_PCanalysis';
file='clust_Kmeans_PC_F1new_corr';
load(fullfile(path1,file))
sel=0; clust=14;


if sel==1
    data=ID(:,clust-1);
end


superpath='e:\all_species_temp';
Dtype='wph';
normDG=[];
for ty=1:3
    load(fullfile(superpath,'newAttempt6','DG',[Dtype(ty),'_normPeaks']))
    normDG=[normDG;normPeaks];
end

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0.2 0.2 0.2; 0 1 1; 1 1 0;0.5 0 0];

clu=unique(data);
for c=1:length(clu)
    currID=find(data==clu(c));
     su=0;
    figure
    for i=1:length(currID)
   
    if su<25
        su=su+1;
        else
            su=1;
            figure
    end
    subplot(5,5,su)
    curr=normDG(currID(i),:,2);
   
    curr=fliplr(reshape(curr,4,6));
    [y,x] = meshgrid(1:6,1:4);
    imagesc(x(:)',y(:)',curr');
    colormap('gray')
    title(['clust ',int2str(c),' / ',int2str(currID(i))])
    end
end
%%!! make one normPeaks matrix
path1='F:\all_species_temp\newAttempt2\DG';
type='wph';
normP=[];
for t=1:3
    load(fullfile(path1,[type(t),'_normPeaks']))
    curr=reshape(normPeaks(:,:,2),size(normPeaks,1),24);
    normP=[normP;curr];
end
save(fullfile(path1,'normPeaks'),'normP')
%% contour plots
clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
load(fullfile('F:\all_species_temp\newAttempt3\DG','normPeaks'))
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
file='clust_Kmeans_PC_F1new';
load(fullfile(path1,file))
sel=1; clust=15;
gau=0;

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0;1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0 0 0; 0 1 1; 1 1 0];
col='bgr';

if sel==1
    data=ID(:,clust-1);
end
%
for d=1:length(data)
    

dat=fliplr(reshape(normP(d,:),4,6));
[X,Y]=meshgrid(1:4,1:6);
% [C h]=contour(X,Y,flipud(dat'));
% hold on
figure
[C h]=contour(X,Y,flipud(dat'),[0.6 0.6],'color',col(tID(d)),'linewidth',1);
close
tmp=find(C(1,:)>=1);
tmp2=find(C(1,tmp)<=4);
tmp=tmp(tmp2);
xax=max(C(1,tmp))-min(C(1,tmp));
tmp=find(C(2,:)>=1);
tmp2=find(C(2,tmp)<=6);
tmp=tmp(tmp2);
yax=max(C(2,tmp))-min(C(2,tmp));
xax=1/4*xax;
yax=1/6*yax;
% per=0.85;
per=1-min([xax yax])*0.6+diff([xax yax])*0.6;
% per=1-min([xax yax])*0.6;
if per>=1
    per=0.99;
elseif per<=0
    per=0.1;
end

figure(100)
subplot(4,5,data(d))
[C h]=contour(X,Y,flipud(dat'),[per per],'color',col(tID(d)),'linewidth',1);
hold on
title(int2str(data(d)))
% close
% figure(100)
% 
% [so id]=sort(C(1,end-20:end),'descend');
% %   id=22-(id); so=fliplr(C(1,end-20:end)); so=so(id);
% tmp=find(so==1);tmp=id(tmp);
% tmp2=find(C(2,end-20:end)==1);
% %   tmp2=20-tmp2;
% tmp=intersect(tmp,tmp2);
% tst=find(id==tmp);
% tst2=find(so(tst+1:end)<1);
% if length(tst2)<2
%     [so id]=sort(C(1,end-30:end),'descend');
%     %   id=22-(id); so=fliplr(C(1,end-20:end)); so=so(id);
%     tmp=find(so==1);tmp=id(tmp);
%     tmp2=find(C(2,end-30:end)==1);
%     %   tmp2=20-tmp2;
%     tmp=intersect(tmp,tmp2);
%     tst=find(id==tmp);
%     tst2=find(so(tst+1:end)<1);
%     tst2=tst2([1 2]);
%     un=unique(so(tst+tst2));
%     val=[];
%     for uu=1:length(un)
%         tt=find(so(tst+tst2)==un(uu));
%         now=min(id(tst+tst2(tt)));
%         val=[val now];
%     end
%     id1=val(1)+length(C)-31;
%     id2=val(2)+length(C)-32;
% else
%     if length(tst2==2)
%     tst2=tst2(end-1:end);
%     else
%       tst2=tst2([1 2]);   
%     end
%     id1=id(tst+tst2(2))+length(C)-21;
%     id2=id(tst+tst2(1))+length(C)-22;
%     if id2<id1
%         [so id]=sort(C(1,end-30:end),'descend');
%         %   id=22-(id); so=fliplr(C(1,end-20:end)); so=so(id);
%         tmp=find(so==1);tmp=id(tmp);
%         tmp2=find(C(2,end-30:end)==1);
%         %   tmp2=20-tmp2;
%         tmp=intersect(tmp,tmp2);
%         tst=find(id==tmp);
%         tst2=find(so(tst+1:end)<1);
%         un=unique(so(tst+tst2));
%         val=[];
%         for uu=1:length(un)
%             tt=find(so(tst+tst2)==un(uu));
%             now=min(id(tst+tst2(tt)));
%             val=[val now];
%         end
%         id1=val(1)+length(C)-31;
%         id2=val(2)+length(C)-32;
%     else
%         un=unique(so(tst+tst2));
%         val=[];
%         for uu=1:length(un)
%             tt=find(so(tst+tst2)==un(uu));
%             now=min(id(tst+tst2(tt)));
%             val=[val now];
%         end
%         id1=val(1)+length(C)-21;
%         id2=val(2)+length(C)-22;
%     end
% end
% 
% 
% 
% m=id1-1; found=0;
% while found==0
%     m=m+1;
%     if abs(diff(C(2,m:m+1)))<C(2,id1)/4
%         found=1;
%     end
% end
% id1=m;

% subplot(4,5,data(d))
% plot(C(1,id1:id2),C(2,id1:id2),'color',col(tID(d)),'linewidth',1)
% hold on
% plot(C(1,id1:id2),C(2,id1:id2),'k','linewidth',3)
% axis([1 4 1 6])
end
%% !!check how similar different clustering algorithms are
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';

load(fullfile(path1,'clust_Kmeans_PC_F1new'))
Kclust=ID;


load(fullfile(path1,'clust_Cmeans_PC_F1new'))
Cclust=ID;


load(fullfile(path1,'clust_MixGau_PC_F1new'))
MGclust=ID;
load(fullfile(path1,'clust_fuzzyGK_PC_F1new'))
GKclust=ID;

takeClust=14;

currK=Kclust(:,takeClust-1);
currC=Cclust(:,takeClust-1);
currMG=MGclust(:,takeClust-1);
currGK=GKclust(:,17-1);

Kmat=zeros(length(currK),length(currK));
for i=1:length(currK)
    clust=currK(i);
    tmp=find(currK==clust);
    Kmat(i,tmp)=1;
end

Cmat=zeros(length(currK),length(currK));
for i=1:length(currK)
    clust=currC(i);
    tmp=find(currC==clust);
    Cmat(i,tmp)=1;
end


MGmat=zeros(length(currK),length(currK));
for i=1:length(currK)
    clust=currMG(i);
    tmp=find(currMG==clust);
    MGmat(i,tmp)=1;
end

% compare Kmat with Cmat

KC=zeros(length(Kmat),1);
for i=1:length(Kmat)
    tmp=find(Kmat(i,:)==1);
    tmp2=find(Cmat(i,tmp)==1);
    per=100/length(tmp)*length(tmp2);
    KC(i)=per;
    
    
end



corrClust=[];
for c=1:takeClust
    tmp=find(currK==c);
    overlap=[];
    for cc=1:takeClust
        tmp2=find(currC==cc);
        ov=intersect(tmp,tmp2);
        overlap=[overlap;length(ov)];
    end
    [ma id]=max(overlap);
    corrClust=[corrClust;id];
   
end


corrclustKC=[]; corrvarsKC=[];overlapKC=[];
for c=1:takeClust
    now=find(currK==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currC==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustKC=[corrclustKC;id];
    corrvarsKC=[corrvarsKC;mi];
     now2=find(currC==id);
     ov=intersect(now,now2);
     overlapKC=[overlapKC;100/length(now)*length(ov)];
end

corrclustKG=[]; corrvarsKG=[];overlapKG=[];
for c=1:takeClust
    now=find(currK==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currGK==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustKG=[corrclustKG;id];
    corrvarsKG=[corrvarsKG;mi];
     now2=find(currGK==id);
     ov=intersect(now,now2);
     overlapKG=[overlapKG;100/length(now)*length(ov)];
end


corrclustKM=[]; corrvarsKM=[];overlapKM=[];
for c=1:takeClust
    now=find(currK==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currMG==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustKM=[corrclustKM;id];
    corrvarsKM=[corrvarsKM;mi];
     now2=find(currMG==id);
    ov=intersect(now,now2);
     overlapKM=[overlapKM;100/length(now)*length(ov)];
end

corrclustGC=[]; corrvarsGC=[]; overlapGC=[];
for c=1:takeClust
    now=find(currGK==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currC==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustGC=[corrclustGC;id];
    corrvarsGC=[corrvarsGC;mi];
    now2=find(currC==id);
    ov=intersect(now,now2);
     overlapGC=[overlapGC;100/length(now)*length(ov)];
end


corrclustGM=[]; corrvarsGM=[]; overlapGM=[];
for c=1:takeClust
    now=find(currGK==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currMG==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustGM=[corrclustGM;id];
    corrvarsGM=[corrvarsGM;mi];
    now2=find(currMG==id);
    ov=intersect(now,now2);
     overlapGM=[overlapGM;100/length(now)*length(ov)];
end


corrclustCG=[]; corrvarsCG=[];overlapCG=[];
for c=1:takeClust
    now=find(currC==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currGK==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustCG=[corrclustCG;id];
    corrvarsCG=[corrvarsCG;mi];
     now2=find(currGK==id);
     ov=intersect(now,now2);
     overlapCG=[overlapCG;100/length(now)*length(ov)];
end


corrclustMG=[]; corrvarsMG=[];overlapMG=[];
for c=1:takeClust
    now=find(currMG==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currGK==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustMG=[corrclustMG;id];
    corrvarsMG=[corrvarsMG;mi];
     now2=find(currGK==id);
     ov=intersect(now,now2);
     overlapMG=[overlapMG;100/length(now)*length(ov)];
end



corrclustCM=[]; corrvarsCM=[]; overlapCM=[];
for c=1:takeClust
    now=find(currC==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currMG==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustCM=[corrclustCM;id];
    corrvarsCM=[corrvarsCM;mi];
    now2=find(currMG==id);
    ov=intersect(now,now2);
     overlapCM=[overlapCM;100/length(now)*length(ov)];
end

   
corrclustMK=[]; corrvarsMK=[]; overlapMK=[];
for c=1:takeClust
    now=find(currMG==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currK==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustMK=[corrclustMK;id];
    corrvarsMK=[corrvarsMK;mi];
    now2=find(currK==id);
    ov=intersect(now,now2);
     overlapMK=[overlapMK;100/length(now)*length(ov)];
end

corrclustGK=[]; corrvarsGK=[]; overlapGK=[];
for c=1:takeClust
    now=find(currGK==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currK==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustGK=[corrclustGK;id];
    corrvarsGK=[corrvarsGK;mi];
    now2=find(currK==id);
    ov=intersect(now,now2);
     overlapGK=[overlapGK;100/length(now)*length(ov)];
end

corrclustMC=[]; corrvarsMC=[]; overlapMC=[];
for c=1:takeClust
    now=find(currMG==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currC==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustMC=[corrclustMC;id];
    corrvarsMC=[corrvarsMC;mi];
    now2=find(currC==id);
    ov=intersect(now,now2);
     overlapMC=[overlapMC;100/length(now)*length(ov)];
end

corrclustCK=[]; corrvarsCK=[]; overlapCK=[];
for c=1:takeClust
    now=find(currC==c);
    currMean=mean(orig2clust(now,:));
    vars=[];
    for cc=1:takeClust
        now2=find(currK==cc);
        curr=sum((mean(orig2clust(now2,:))-currMean).^2);
        vars=[vars;curr];
    end
    [mi id]=min(vars);
    corrclustCK=[corrclustCK;id];
    corrvarsCK=[corrvarsCK;mi];
    now2=find(currK==id);
    ov=intersect(now,now2);
     overlapCK=[overlapCK;100/length(now)*length(ov)];
end
%% find perfect overlaps
clear
close all
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';

load(fullfile(path1,'clust_Kmeans_PC_F1new'))
Kclust=ID;


load(fullfile(path1,'clust_Cmeans_PC_F1new'))
Cclust=ID;


load(fullfile(path1,'clust_MixGau_PC_F1new'))
MGclust=ID;
load(fullfile(path1,'clust_fuzzyGK_PC_F1new'))
GKclust=ID;

takeClust=14;

currK=Kclust(:,takeClust-1);
currC=Cclust(:,takeClust-1);
currMG=MGclust(:,takeClust-1);
currGK=GKclust(:,17-1);

optClust=zeros(14,14,14,14);

for k=1:14
    tmp1=find(currK==k);
    for c=1:14
        tmp2=find(currC==c);
        for m=1:14
            tmp3=find(currMG==m);
            for g=1:14
                
                tmp4=find(currGK==g);
                now=unique([tmp1;tmp2;tmp3; tmp4]);
                comp=sum(sum((orig2clust(now,:)-repmat(mean(orig2clust(now,:)),length(now),1)).^2,1))/length(now);
                optClust(k,c,m,g)=comp;
            end
        end
    end
end

finalClust=zeros(14,4);
finalVar=zeros(14,1);
for cl=1:14
    now=reshape(optClust(cl,:,:,:),14,14,14);
    mi=min(min(min(now)));
    
    [M3 idx3] = min(now, [], 3);
    [M2 idx2] = min(M3, [], 2);
    [M1 idx1] = min(M2, [], 1);
    idx1;
    idx2=idx2(idx1);
    idx3=idx3(idx1,idx2);
    finalClust(cl,1)=cl;
    finalClust(cl,2)=idx1;
    finalClust(cl,3)=idx2;
    finalClust(cl,4)=idx3;
    finalVar(cl)=min(now(:));
end
%% !!check clustering results
clear
close all

% clusteringResult=[1 1 6;2 13 8;3 2 9;4 15 3;5 10 13;6 7 3;7 11 1;8 5 12;9 9 5;10 4 7;11 3 16;12 6 12;13 16 14;14 14 11;15 2 15;16 15 2];
% path1='F:\all_species_temp\newAttempt2\DG_PCanalysis';

clusteringResult=[[1 11 12;2 9 5;3 14 13;4 7 2;5 10 1;6 13 9;7 8 7;8 6 10;9 8 14;10 5 3;11 12 4;12 3 11;13 5 8;14 10 3]];
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';

load(fullfile(path1,'clust_Kmeans_PC_F1new'))
Kclust=ID;


load(fullfile(path1,'clust_Cmeans_PC_F1new'))
Cclust=ID;


load(fullfile(path1,'clust_MixGau_PC_F1new'))
MGclust=ID;

load(fullfile(path1,'clust_fuzzyGK_PC_F1new'))
GKclust=ID;

takeClust=14;

currK=Kclust(:,takeClust-1);
currC=Cclust(:,takeClust-1);
currMG=MGclust(:,takeClust-1);
currGK=GKclust(:,takeClust-1);

grandMean=sum(sum((orig2clust-repmat(mean(orig2clust),length(orig2clust),1)).^2,1))/length(orig2clust);

for c=1:takeClust
    
    tmp=find(currK==c);
    tmp2=find(currC==clusteringResult(c,2));
    tmp3=find(currMG==clusteringResult(c,3));
    tmp4=find(currGK==clusteringResult(c,4));
    
%      vars=[];
%     for c=1:length(uclust)
%         now=find(data==c);
%         curr=length(now).*(mean(orig2clust(now,:))-grandMean).^2;
%         vars=[vars;curr];
%     end
%     vars=sum(vars)/(length(uclust)-1);
%     between=[between;vars];
%     
%     vars=[];
%     for c=1:length(uclust)
%         now=find(data==c);
%         curr=sum((orig2clust(now,:)-repmat(mean(orig2clust(now,:)),length(now),1)).^2,1);
%         vars=[vars;curr];
%     end
%     vars=sum(vars)/(size(data,1)-length(uclust));
%     within=[within;vars];
    now=unique([tmp;tmp2;tmp3]);
comp=sum(sum((orig2clust(now,:)-repmat(mean(orig2clust(now,:)),length(now),1)).^2,1))/length(now);    
comp0=sum(sum((orig2clust(tmp,:)-repmat(mean(orig2clust(tmp,:)),length(tmp),1)).^2,1))/length(tmp);   
    comp1=sum(mean(orig2clust(tmp,:))-mean(orig2clust(tmp2,:)).^2)/(length(tmp)+length(tmp2))*100;
    comp2=sum(mean(orig2clust(tmp,:))-mean(orig2clust(tmp3,:)).^2)/(length(tmp)+length(tmp3))*100;
    comp3=sum(mean(orig2clust(tmp2,:))-mean(orig2clust(tmp3,:)).^2)/(length(tmp3)+length(tmp2))*100;
    
    figure(1)
    subplot(4,4,c)
    plot(orig2clust(tmp,1),orig2clust(tmp,2),'or','markerfacecolor','r')
    hold on
%     title([int2str(c),' (',num2str(comp1),'/',num2str(comp2),'/',num2str(comp3),')'])
title([int2str(c),' (',num2str(100/grandMean*comp),'/',num2str(100/grandMean*comp0),')'])
    plot(orig2clust(tmp2,1),orig2clust(tmp2,2),'ob','markerfacecolor','b')
    plot(orig2clust(tmp3,1),orig2clust(tmp3,2),'og','markerfacecolor','g')
    
    
    figure(2)
    subplot(4,4,c)
    plot(orig2clust(tmp,3),orig2clust(tmp,4),'or','markerfacecolor','r')
    hold on
    title([int2str(c),' (',num2str(comp1),'/',num2str(comp2),'/',num2str(comp3),')'])
    plot(orig2clust(tmp2,3),orig2clust(tmp2,4),'ob','markerfacecolor','b')
    plot(orig2clust(tmp3,3),orig2clust(tmp3,4),'og','markerfacecolor','g')
    
    figure(3)
    subplot(4,4,c)
    plot(orig2clust(tmp,6),orig2clust(tmp,7),'or','markerfacecolor','r')
    hold on
    title([int2str(c),' (',num2str(comp1),'/',num2str(comp2),'/',num2str(comp3),')'])
    plot(orig2clust(tmp2,6),orig2clust(tmp2,7),'ob','markerfacecolor','b')
    plot(orig2clust(tmp3,6),orig2clust(tmp3,7),'og','markerfacecolor','g')
    
    int1=intersect(tmp,tmp2);
    int2=intersect(int1,tmp3);%these are in all 3
    
    
    
    co1=intersect(tmp,tmp2);
    co1=setdiff(co1,int2); % in Kmeans and Cmean
    
    co2=intersect(tmp,tmp3);
    co2=setdiff(co2,int2); % in Kmeans and MG
    
    co3=intersect(tmp3,tmp2);
    co3=setdiff(co3,int2); % in Cmeans and MG
    
    figure(1)
    plot(orig2clust(int2,1),orig2clust(int2,2),'ok','markerfacecolor','k')
    plot(orig2clust(co1,1),orig2clust(co1,2),'om','markerfacecolor','m')
    plot(orig2clust(co2,1),orig2clust(co2,2),'oc','markerfacecolor','c')
    plot(orig2clust(co3,1),orig2clust(co3,2),'oy','markerfacecolor','y')
    plot(orig2clust(int2,1),orig2clust(int2,2),'ok','markerfacecolor','k')
    axis([0 1 0 1])
    
    figure(2)
    plot(orig2clust(int2,3),orig2clust(int2,4),'ok','markerfacecolor','k')
    plot(orig2clust(co1,3),orig2clust(co1,4),'om','markerfacecolor','m')
    plot(orig2clust(co2,3),orig2clust(co2,4),'oc','markerfacecolor','c')
    plot(orig2clust(co3,3),orig2clust(co3,4),'oy','markerfacecolor','y')
    plot(orig2clust(int2,3),orig2clust(int2,4),'ok','markerfacecolor','k')
    axis([0 1 0 1])
    
    
    
    figure(3)
    plot(orig2clust(int2,6),orig2clust(int2,7),'ok','markerfacecolor','k')
    plot(orig2clust(co1,6),orig2clust(co1,7),'om','markerfacecolor','m')
    plot(orig2clust(co2,6),orig2clust(co2,7),'oc','markerfacecolor','c')
    plot(orig2clust(co3,6),orig2clust(co3,7),'oy','markerfacecolor','y')
    plot(orig2clust(int2,6),orig2clust(int2,7),'ok','markerfacecolor','k')
    axis([0 1 0 1])
end
%% !!extract the needed clusters
clear
close all
clust=15;
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';

load(fullfile(path1,'clust_Kmeans_PC_F1new'))
ID=ID(:,clust-1);
save(fullfile(path1,'clust_Kmeans_PC_F1new_sel'),'ID','tID','keepOrig','orig2clust')


load(fullfile(path1,'clust_Cmeans_PC_F1new'))
ID=ID(:,clust-1);
save(fullfile(path1,'clust_Cmeans_PC_F1new_sel'),'ID','tID','keepOrig','orig2clust')


load(fullfile(path1,'clust_MixGau_PC_F1new'))
ID=ID(:,clust-1);
save(fullfile(path1,'clust_MixGau_PC_F1new_sel'),'ID','tID','keepOrig','orig2clust')

load(fullfile(path1,'clust_fuzzyGK_PC_F1new'))
ID=ID(:,15);
save(fullfile(path1,'clust_fuzzyGK_PC_F1new_sel'),'ID','tID','keepOrig','orig2clust')
%% !!Rand index
clear 
close all
path1='F:\all_species_temp\newAttempt4\DG_PCanalysis';
load(fullfile(path1,'clust_Kmeans_PC_F1new_sel'))
data1=ID;

load(fullfile(path1,'clust_Cmeans_PC_F1new_sel'))
data2=ID;

load(fullfile(path1,'clust_MixGau_PC_F1new_sel'))
data3=ID;

load(fullfile(path1,'clust_fuzzyGK_PC_F1new_sel'))
data4=ID;

mockup=randi(14,635,1);

addpath(genpath('C:\Program Files\MATLAB\R2014a\toolbox\NCestimation_V2'))
clear AR Ri MI HI
[AR,RI,MI,HI]=RandIndex(data4,mockup);
%% !!Rand index across clustering attempts
clear 
close all
path1='F:\all_species_temp\newAttempt3\DG_PCanalysis';
load(fullfile(path1,'clust_Kmeans_PC_F1new_sel'))
data1=ID;

load(fullfile(path1,'clust_Cmeans_PC_F1new_sel'))
data2=ID;

load(fullfile(path1,'clust_MixGau_PC_F1new_sel'))
data3=ID;

load(fullfile(path1,'clust_fuzzyGK_PC_F1new_sel'))
data4=ID;



path1='F:\all_species_temp\newAttempt4\DG_PCanalysis';
load(fullfile(path1,'clust_Kmeans_PC_F1new_sel'))
data11=ID;

load(fullfile(path1,'clust_Cmeans_PC_F1new_sel'))
data22=ID;

load(fullfile(path1,'clust_MixGau_PC_F1new_sel'))
data33=ID;

load(fullfile(path1,'clust_fuzzyGK_PC_F1new_sel'))
data44=ID;



mockup=randi(14,635,1);

addpath(genpath('C:\Program Files\MATLAB\R2014a\toolbox\NCestimation_V2'))
clear AR Ri MI HI
[AR,RI,MI,HI]=RandIndex(data4,data44);
%% ....just for testing
% inital kmeans step used to initialize EM
K = 2;               % number of mixtures/clusters
cInd = kmeans(data, K, 'EmptyAction','singleton');
% fit a GMM model
gmm = fitgmdist(data, K, 'Options',statset('MaxIter',1000));
% means, covariances, and mixing-weights
mu = gmm.mu;
sigma = gmm.Sigma;
p = gmm.PComponents;
% cluster and posterior probablity of each instance
% note that: [~,clustIdx] = max(p,[],2)
[clustInd,~,p] = cluster(gmm, data);
tabulate(clustInd)
% plot data, clustering of the entire domain, and the GMM contours
clrLite = [1 0.6 0.6 ; 0.6 1 0.6 ; 0.6 0.6 1];
clrDark = [0.7 0 0 ; 0 0.7 0 ; 0 0 0.7];
[X,Y] = meshgrid(linspace(mn(1),mx(1),100), linspace(mn(2),mx(2),100));
% C = cluster(gmm, [X(:) Y(:)]);
% image(X(:), Y(:), reshape(C,size(X))), hold on
% gscatter(data(:,1), data(:,2), species, clrDark)
gscatter(data(:,1), data(:,2))
h = ezcontour(@(x1,x2)pdf(gmm,[x1 x2]), [mn(1) mx(1) mn(2) mx(2)]);
%% ----HERE STARTS THE PER EXPERIMENT ANALYSIS PART ----
%% !!speed tuning from DG
clear
close all
superpath='f:\all_species_temp';
cols='bgr';
Dtype='wph';
spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
speed=spat.*freq/1000; speed=speed';
% uspeed=unique(speed);
% uspeed=[0.5 1 2 4 8 16 32];
uspeed=[1 2 4 8 16];
% uspeed=[1 2 4 8];

ids=[];
for u=1:length(uspeed)
    tmp=find(speed==uspeed(u));
    ids=[ids;tmp];
end

speed=speed(ids);
figure
for d=1:3
    
    load(fullfile(superpath,'newAttempt2','DG',[Dtype(d),'_FFTpeaks']))
    data=reshape(FFTpeaks(:,ids,2),size(FFTpeaks,1),length(ids));
    data=data./repmat(max(data,[],2),1,length(ids));
    
    
    
    newdata=[];
    for s=1:length(uspeed)
        tmp=find(speed==uspeed(s));
        newdata=[newdata max(data(:,tmp),[],2)];
    end
    
    
    DGspeed=[];
    for c=1:size(newdata,1)
        su=sum(newdata(c,:));
        cumdata=cumsum(newdata(c,:));
        tmp1=find(cumdata<su/2);
        if isempty(tmp1)
            optDG=uspeed(1);
        else
            tmp1=tmp1(end);
            tmp2=find(cumdata>su/2);
            if isempty(tmp2)
                optDG=uspeed(end);
            else
                tmp2=tmp2(1);
                slope=(cumdata(tmp2)-cumdata(tmp1))/(uspeed(tmp2)-uspeed(tmp1));
                optDG=(su/2-cumdata(tmp1))/slope+uspeed(tmp1);
            end
        end
        DGspeed=[DGspeed;optDG];
    end
    if d==1
        mall=DGspeed;
    elseif d==2
        pall=DGspeed;
    else
        hall=DGspeed;
    end
    save(fullfile(superpath,[Dtype(d),'_DGspeed']),'DGspeed','togo')

    DGspeed=[DGspeed;-1;max(uspeed)+1];
   
    figure(1)
    subplot(3,2,d*2-1)
    hist(DGspeed,38)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cols(d),'EdgeColor',cols(d))
    med=median(DGspeed(1:end-2));
    st=std(DGspeed(1:end-2));
    conf = 1.96*st/sqrt(length(DGspeed)-2);
    title([Dtype(d), ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(DGspeed)-2),')'])
    hold on
    plot([med-st med+st],[0.8 0.8],'-k','linewidth',2)
    plot(med,0.8,'d','markerfacecolor','k','markersize',7,'markeredgecolor','k')
    yy=get(gca,'ylim');
    plot([med-conf med-conf],[0 yy(2)+1],'--k')
    plot([med+conf med+conf],[0 yy(2)+1],'--k')
    %     set(gca,'xtick',speeds)
    % set(gca,'xticklabel',speeds)
    xlabel('speed (mm/s)')
    axis([0 uspeed(end) 0 yy(2)+1])
    
    subplot(3,2,[2 4])
    c=cdfplot(DGspeed(1:end-2))
    set(c,'color',cols(d),'linestyle','-','linewidth',2)
    hold on
    xlabel('speed (mm/s)')
    ylabel('probability')
    title('DG speed tuning')
    
    figure(2)
    collme=[]; collst=[]; collci=[];
    for s=1:length(uspeed)
        tmp=find(speed==uspeed(s));
        %         tmp=find(speed(13:16)==uspeed(s)); tmp=tmp+12;
        currdata=max(data(:,tmp),[],2);
        %         currdata=mean(data(:,tmp),2);
        me=mean(currdata);
        st=std(currdata);
        ci=1.96*st/sqrt(size(currdata,1));
        collme=[collme;me];
        collst=[collst;st];
        collci=[collci;ci];
        subplot(2,2,d)
        plot([uspeed(s) uspeed(s)],[me-st me+st],'-','color',cols(d))
        hold on
        plot([uspeed(s) uspeed(s)],[me-ci me+ci],'-','color',cols(d),'linewidth',4)
        subplot(2,2,4)
        plot([uspeed(s) uspeed(s)],[me-st me+st],'-','color',cols(d))
        hold on
        plot([uspeed(s) uspeed(s)],[me-ci me+ci],'-','color',cols(d),'linewidth',4)
    end
    
    subplot(2,2,d)
    plot(uspeed,collme,'o-','color',cols(d),'markerfacecolor',cols(d),'markeredgecolor',cols(d),'linewidth',2)
    subplot(2,2,4)
    plot(uspeed,collme,'o-','color',cols(d),'markerfacecolor',cols(d),'markeredgecolor',cols(d),'linewidth',2)
end
p1=ranksum(mall, pall);
p2=ranksum(mall, hall);
p3=ranksum(pall, hall);
%% !!replot speed tuninng 6vel for max(black,white)
clear
close all
superpath='f:\all_species_temp';
savepath=fullfile(superpath,'overview_pics');


names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[0.2; 1; 2; 4; 8; 16];
types='wph';
cols='bgr';
titles={'mouse','pig','human'};

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
    
    allval=zeros(size(info,1),2);
    black=cell2mat(info(goodB,24));
    white=cell2mat(info(goodW,25));
    allval(goodB,1)=black;
    allval(goodW,2)=white;
    allval=max(allval,[],2);
    tmp=find(allval>0);
    allval=allval(tmp);
    if t==1
        mall=allval;
    elseif t==2
        pall=allval;
    else
        hall=allval;
    end
    
    allval=[allval;-1;17];
    
    
    
    subplot(3,2,t*2-1)
    hist(allval,38)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
    med=median(allval(1:end-2));
    st=std(allval(1:end-2));
    conf = 1.96*st/sqrt(length(allval)-2);
    title([titles{t}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(allval)-2),')'])
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
    c=cdfplot(allval(1:end-2))
    set(c,'color',cols(t),'linestyle','-','linewidth',2)
    hold on
    xlabel('speed (mm/s)')
    ylabel('probability')
    title('better bar')
    
    
end

set(gcf,'position',[1 31 1600 794])
saveas(gcf,fullfile(savepath,'speed_all_noSlow.bmp'))
saveas(gcf,fullfile(savepath,'speed_all_noSlow.fig'))
% close

p1=ranksum(mall, pall);
p2=ranksum(mall, hall);
p3=ranksum(pall, hall);


%% 6vel per exp
clear
close all
superpath='f:\all_species_temp';
savepath=fullfile(superpath,'overview_pics');

allmeans=[2.0015 3.1624 3.7842];

names={'0.2 mm/s';'1 mm/s';'2 mm/s';'4 mm/s';'8 mm/s';'16 mm/s';};
speeds=[0.2; 1; 2; 4; 8; 16];
types='wph';
cols='bgr';
titles={'mouse','pig','human'};

collMean=zeros(20,3); NofE=zeros(3,1);
alldata=cell(3,20);
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
    
    allval=zeros(size(info,1),2);
    black=cell2mat(info(goodB,24));
    white=cell2mat(info(goodW,25));
    allval(goodB,1)=black;
    allval(goodW,2)=white;
    allvalorig=max(allval,[],2);
    
    dates=info(:,1);
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    start=[start;length(allvalorig)];
    udates=udates(ia);
    NofE(t)=length(udates);
    
    su=length(udates)/4;
    if mod(su,1)==0
        su=su+1;
    else
        su=ceil(su);
    end
    collmeds=[];
    for u=1:length(udates)
        allval=allvalorig(start(u):start(u+1)-1);
        tmp=find(allval>0);
        allval=allval(tmp);
        alldata{t,u}=allval;
        allval=[allval;-1;17];
        
        
        if length(allval)>2
            figure(t)
            subplot(4,su,u)
            hist(allval,38)
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor',cols(t),'EdgeColor',cols(t))
            med=median(allval(1:end-2));
            collMean(u,t)=med;
            st=std(allval(1:end-2));
            conf = 1.96*st/sqrt(length(allval)-2);
            title([udates{u}, ': ',num2str(med),' +/- ',num2str(st),' (total: ',int2str(length(allval)-2),')'])
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
            
            figure(t)
            subplot(4,su,su*4)
            c=cdfplot(allval(1:end-2))
            set(c,'color',cols(t),'linestyle','-','linewidth',2)
            hold on
            xlabel('speed (mm/s)')
            ylabel('probability')
            title('better bar')
            figure(4)
            subplot(1,2,1)
            c=cdfplot(allval(1:end-2))
            set(c,'color',cols(t),'linestyle','-','linewidth',2)
            hold on
            xlabel('speed (mm/s)')
            ylabel('probability')
            title('better bar')
            figure(4)
            subplot(1,2,2)
            plot(t,med,'o','markerfacecolor',cols(t),'markeredgecolor',cols(t))
            hold on
            %             plot([t t],[med-st med+st],'-','color',cols(t))
            %             plot([t t],[med-conf med+conf],'-','color',cols(t),'linewidth',2)
            axis([0 4 0 inf])
            collmeds=[collmeds;med];
        end
    end
    plot([t-0.2 t+0.2],[mean(collmeds) mean(collmeds)],'-','color',cols(t),'linewidth',2) %mean of means
    plot([t t],[mean(collmeds)-std(collmeds) mean(collmeds)+std(collmeds)],'-','color',cols(t)) %mean of means
    plot([t-0.2 t+0.2],[allmeans(t) allmeans(t)],'-','color',cols(t)) %previously calculated population mean
end

p1=ranksum(collMean(1:NofE(1),1),collMean(1:NofE(2),2)); %means for M vs P
p2=ranksum(collMean(1:NofE(1),1),collMean(1:NofE(3),3));  %means for M vs H
p3=ranksum(collMean(1:NofE(2),2),collMean(1:NofE(3),3));  %means for P vs H
% set(gcf,'position',[1 31 1600 794])
% saveas(gcf,fullfile(savepath,'speed_all_noSlow.bmp'))
% saveas(gcf,fullfile(savepath,'speed_all_noSlow.fig'))
% close

% within species
goodExps=cell(3,1); pall=cell(3,1);

for t=1:3
    goodExp=[];
    for e=1:20
        if ~isempty(alldata{t,e})
            goodExp=[goodExp;e];
        end
    end
    goodExps{t}=goodExp;
    pvals=[];
    for g=1:length(goodExp)-1
        for g2=g+1:length(goodExp)
            p=ranksum(alldata{t,goodExp(g)},alldata{t,goodExp(g2)});
            pvals=[pvals;p];
        end
    end
    pall{t}=pvals;
    matP=zeros(length(goodExp),length(goodExp));
    numOfComp=[0 length(goodExp)-1:-1:1];
    for g=1:length(goodExp)-1
matP(g,g)=1;        
            currID1=sum(numOfComp(1:g))+1;
            currID2=sum(numOfComp(1:g+1));
            matP(g,g+1:end)=pvals(currID1:currID2);
matP(g+1:end,g)=pvals(currID1:currID2);
    end
    matP(g+1,g+1)=1;
    matP=1-matP;
    tmp1=find(matP>0.95);
    tmp2=find(matP>0.99); 
    tmp3=find(matP<0.95); matP(tmp3)=0;
    matP(tmp1)=0.8; matP(tmp2)=1;
figure
imagesc(matP,[0 1])
colormap('gray')
title('p-values across mouse experiments')
end
figure
subplot(1,2,1)
for t=1:3
    plot(t,pall{t},'o','color',cols(t))
    hold on
    plot([t-0.2 t+0.2],[median(pall{t}) median(pall{t})],'-','color',cols(t))
end

plot([0 4],[0.05 0.05],'-k')
plot([0 4],[0.01 0.01],'-k')
axis([0 4 0 1])
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'m','p','h'})
title('ranksum p-values between all experiments within each species')



% across species
pcross=cell(3,1); cnt=0;
for t=1:2
    for tt=t+1:3
        pvals=[]; cnt=cnt+1;
        for g=1:length(goodExps{t})
            for g2=1:length(goodExps{tt})
                p=ranksum(alldata{t,goodExps{t}(g)},alldata{tt,goodExps{tt}(g2)});
                pvals=[pvals;p];
            end
        end
        pcross{cnt}=pvals;
        malength=max(length(goodExps{t}),length(goodExps{tt}));
         matP=zeros(length(goodExps{t}),length(goodExps{tt}));
    numOfComp=[0 length(goodExp)-1:-1:1];
    
    end
end


subplot(1,2,2)
for t=1:3
    plot(t,pcross{t},'o','color',cols(t))
    hold on
    plot([t-0.2 t+0.2],[median(pcross{t}) median(pcross{t})],'-','color',cols(t))
end

plot([0 4],[0.05 0.05],'-k')
plot([0 4],[0.01 0.01],'-k')
axis([0 4 0 1])
set(gca,'xtick',1:3)
set(gca,'xticklabel',{'m-p','m-h','p-h'})
title('ranksum p-values between all experiments across species')
%% spatial / temporal tuning per exp
clear
close all
superpath='f:\all_species_temp';
version=1;

col1={'bbcc','ggkk','rrmm'};
col2={'bbbccc','gggkkk','rrrmmm'};
style1={'-','--','-','--'};
style2={'-','--','-.','-','--','-.'};
cols='bgr';
Dtype='wph';
spatials=[100 200 500 1000 2000 4000];
freqs=[1 2 4 8];
harmlist=[0.5 1 2 3];

collNorm=cell(3,1);

meanPeakT=zeros(3,20,6,4);
meanPeakS=zeros(3,20,4,6);
dataExpT=cell(3,20,6);
dataExpS=cell(3,20,4);
NofE=zeros(3,1);
udat=cell(3,1);
alldata=cell(3,20);
for d=1:3
    
    load(fullfile(superpath,'newAttempt','DG',[Dtype(d),'_normPeaks']))
    collNorm{d,1}=reshape(normPeaks(:,:,2),size(normPeaks,1),24);
    load(fullfile(superpath,[Dtype(d),'_info']))
    dates=info(togo,1);
    [udates,ia] =unique(dates);
    [start,ia]=sort(ia);
    start=[start;length(normPeaks)];
    NofE(d)=length(udates);
    udat{d}=udates;
    for u=1:length(udates)
        
        alldata{d,u}=normPeaks(start(u):start(u+1)-1,:,2);
        
        
        scnt=0;
        for s=1:4:24 %temporal tuning
            scnt=scnt+1;
            
            figure(scnt)
            meanPeak=mean(normPeaks(start(u):start(u+1)-1,s:s+3,2));
            st=std(normPeaks(start(u):start(u+1)-1,s:s+3,2));
            ci=1.96*st/sqrt(length(start(u):start(u+1)-1));
            for g=1:4
                NofCT(d,scnt,g)=size(normPeaks,1);
                goodT{d,scnt,g}=1:size(normPeaks,1);
            end
            
            plot((freqs)+(d-1)*0.1,meanPeak,'color',cols(d),'linewidth',1)
            hold on
            
            
            xlabel('frequency [Hz]')
            axis([0 9 -inf inf])
            title([int2str(spatials(scnt)),'um / f',num2str(harmlist(2))])
            
            dataExpT{d,u,scnt}=normPeaks(start(u):start(u+1)-1,s:s+3,2);
            meanPeakT(d,u,scnt,:)=meanPeak;
            
            
        end
        scnt=0;
        for s=1:4
            scnt=scnt+1;
            figure(scnt+6)
            
            meanPeak=mean(normPeaks(start(u):start(u+1)-1,s:4:24,2));
            st=std(normPeaks(start(u):start(u+1)-1,s:4:24,2));
            ci=1.96*st/sqrt(length(start(u):start(u+1)-1));
            
            
            semilogx(([100 200 500 1000 2000 4000])+spatials*0.1,meanPeak,'color',cols(d),'linewidth',1)
            hold on
            
            xlabel('spatial period [um]')
            axis([50 5000 -inf inf])
            title([int2str(freqs(s)),'Hz / f',num2str(harmlist(2))])
            
            dataExpS{d,u,scnt}=normPeaks(start(u):start(u+1)-1,s:4:24,2);
            meanPeakS(d,u,scnt,:)=meanPeak;
        end
        
    end
    
end

for ss=1:6
    figure(ss)
    for d=1:3
        currmean=mean(reshape(meanPeakT(d,1:NofE(d),ss,:),NofE(d),4));
        plot((freqs)+(d-1)*0.1,currmean,'color',cols(d),'linewidth',4)
    end
end
for ss=1:4
    figure(ss+6)
    for d=1:3
        currmean=mean(reshape(meanPeakS(d,1:NofE(d),ss,:),NofE(d),6));
        semilogx(([100 200 500 1000 2000 4000])+spatials*0.1,currmean,'color',cols(d),'linewidth',4)
    end
end


for d=1:3
    figure(100+d)
    for u=1:NofE(d)
        
        subplot(4,4,u)
        for fcnt=1:6
            plot((freqs)+(fcnt-7)*0.1,reshape(meanPeakT(d,u,fcnt,:),4,1),'color',col2{d}(fcnt),'linewidth',2,'linestyle',style2{fcnt})
            hold on
        end
        
        xlabel('frequency [Hz]')
        axis([0 9 -inf inf])
        title(['temporal tuning (',udat{d}{u},')'])
    end
    
    
    
    figure(1000+d)
    for u=1:NofE(d)
        
        subplot(4,4,u)
        for fcnt=1:4
            semilogx((spatials)+spatials*0.05*(fcnt-3),reshape(meanPeakS(d,u,fcnt,:),6,1),'color',col1{d}(fcnt),'linewidth',2,'linestyle',style1{fcnt})
            hold on
        end
        
        
        xlabel('spatial period [um]')
        axis([50 5000 -inf inf])
        title(['spatial tuning (',udat{d}{u},')'])
    end
end


for d=1:3
    for u=1:NofE(d)
        tdata=alldata{d,u};
        
        meandata=[];
        for m=1:24
            curr=tdata(:,m);
            meandata=[meandata mean(curr)];
        end
        
        normthis=meandata/max(meandata);
        zi=reshape(normthis,length(unique(freqs)),length(unique(spatials)));
        zi=zi';
        zi=flipud(zi);
        figure(10000+d)
        subplot(4,4,u)
        set(gcf,'position',[1 31 1600 794])
        imagesc(1:length(unique(freqs)),1:length(unique(spatials)),zi)
        colormap(gray)
        set(gca,'YTick',1:length(unique(spatials)))
        set(gca,'XTick',1:length(unique(freqs)))
        set(gca,'XTickLabel',unique(freqs))
        xlabel('Hz')
        set(gca,'YTickLabel',fliplr(unique(spatials)))
        ylabel('period (um)')
        title([udat{d}{u},' (',int2str(size(tdata,1)),'cells)'])
    end
end


% ranksum temporal within
pall=cell(3,1);
for d=1:3
    pvals=zeros(200,6,4);cnt=0;
    for u=1:NofE(d)-1
        for uu=u+1:NofE(d)
            cnt=cnt+1;
            for s=1:6
                for t=1:4
                    p=ranksum(dataExpT{d,u,s}(:,t),dataExpT{d,uu,s}(:,t));
                    pvals(cnt,s,t)=p;
                end
            end
        end
    end
    pvals=pvals(1:cnt,:,:);
    pall{d}=pvals;
end
figure(50)
for d=1:3
    for s=1:6
        subplot(2,3,s)
        plot((freqs)+(d-1)*0.1, reshape(pall{d}(:,s,:),size(pall{d},1),4),'o','color',cols(d))
        hold on
        meanP=[];
        for i=1:4
            tmp=pall{d}(:,s,i);
            tst=tmp(~isnan(tmp));
            meanP=[meanP median(tst)];
            if median(tst)<0.05
                plot([(freqs(i))-0.4 (freqs(i))+0.4], [median(tst) median(tst)],'-','color',cols(d),'linewidth',3)
            else
                plot([(freqs(i))-0.4 (freqs(i))+0.4], [median(tst) median(tst)],'-','color',cols(d))
            end
        end
        
        plot([0 9],[0.05 0.05],'-k')
        plot([0 9],[0.01 0.01],'-k')
        title([int2str(spatials(s)),'um '])
    end
end


% ranksum spatial within
pall=cell(3,1);
for d=1:3
    pvals=zeros(200,4,6);cnt=0;
    for u=1:NofE(d)-1
        for uu=u+1:NofE(d)
            cnt=cnt+1;
            for t=1:4
                for s=1:6
                    p=ranksum(dataExpS{d,u,t}(:,s),dataExpS{d,uu,t}(:,s));
                    pvals(cnt,t,s)=p;
                end
            end
        end
    end
    pvals=pvals(1:cnt,:,:);
    pall{d}=pvals;
end
figure(51)
for d=1:3
    for t=1:4
        subplot(2,2,t)
        semilogx(([100 200 500 1000 2000 4000])+spatials*(d-1)*0.1,reshape(pall{d}(:,t,:),size(pall{d},1),6),'o','color',cols(d))
        hold on
        meanP=[];
        for i=1:6
            tmp=pall{d}(:,t,i);
            tst=tmp(~isnan(tmp));
            meanP=[meanP median(tst)];
            if median(tst)<0.05
                plot([(spatials(i))-spatials(i)*0.1 (spatials(i))+spatials(i)*0.1], [median(tst) median(tst)],'-','color',cols(d),'linewidth',3)
            else
                plot([(spatials(i))-spatials(i)*0.1 (spatials(i))+spatials(i)*0.1], [median(tst) median(tst)],'-','color',cols(d))
            end
        end
        
        plot([50 5100],[0.05 0.05],'-k')
        plot([50 5100],[0.01 0.01],'-k')
        title([int2str(freqs(t)),'Hz '])
        axis([50 5100 0 1])
    end
end


% ranksum temporal between
pall=cell(3,1); dcnt=0;
for d=1:2
    for dd=d+1:3
        dcnt=dcnt+1;
        pvals=zeros(400,6,4);cnt=0;
        for u=1:NofE(d)-1
            for uu=1:NofE(dd)
                cnt=cnt+1;
                for s=1:6
                    for t=1:4
                        p=ranksum(dataExpT{d,u,s}(:,t),dataExpT{dd,uu,s}(:,t));
                        pvals(cnt,s,t)=p;
                    end
                end
            end
        end
        pvals=pvals(1:cnt,:,:);
        pall{dcnt}=pvals;
    end
end
figure(52)
for d=1:3
    for s=1:6
        subplot(2,3,s)
        plot((freqs)+(d-1)*0.1, reshape(pall{d}(:,s,:),size(pall{d},1),4),'o','color',cols(d))
        hold on
        meanP=[];
        for i=1:4
            tmp=pall{d}(:,s,i);
            tst=tmp(~isnan(tmp));
            meanP=[meanP median(tst)];
            if median(tst)<0.05
                plot([(freqs(i))-0.4 (freqs(i))+0.4], [median(tst) median(tst)],'-','color',cols(d),'linewidth',3)
            else
                plot([(freqs(i))-0.4 (freqs(i))+0.4], [median(tst) median(tst)],'-','color',cols(d))
            end
        end
        
        plot([0 9],[0.05 0.05],'-k')
        plot([0 9],[0.01 0.01],'-k')
        title([int2str(spatials(s)),'um '])
    end
end


% ranksum spatial between
pall=cell(3,1);dcnt=0;
for d=1:2
    for dd=d+1:3
        dcnt=dcnt+1;
        pvals=zeros(200,4,6);cnt=0;
        for u=1:NofE(d)-1
            for uu=1:NofE(dd)
                cnt=cnt+1;
                for t=1:4
                    for s=1:6
                        p=ranksum(dataExpS{d,u,t}(:,s),dataExpS{dd,uu,t}(:,s));
                        pvals(cnt,t,s)=p;
                    end
                end
            end
        end
        pvals=pvals(1:cnt,:,:);
        pall{dcnt}=pvals;
    end
end
figure(53)
for d=1:3
    for t=1:4
        subplot(2,2,t)
        semilogx(([100 200 500 1000 2000 4000])+spatials*(d-1)*0.1,reshape(pall{d}(:,t,:),size(pall{d},1),6),'o','color',cols(d))
        hold on
        meanP=[];
        for i=1:6
            tmp=pall{d}(:,t,i);
            tst=tmp(~isnan(tmp));
            meanP=[meanP median(tst)];
            if median(tst)<0.05
                plot([(spatials(i))-spatials(i)*0.1 (spatials(i))+spatials(i)*0.1], [median(tst) median(tst)],'-','color',cols(d),'linewidth',3)
            else
                plot([(spatials(i))-spatials(i)*0.1 (spatials(i))+spatials(i)*0.1], [median(tst) median(tst)],'-','color',cols(d))
            end
        end
        
        plot([50 5100],[0.05 0.05],'-k')
        plot([50 5100],[0.01 0.01],'-k')
        title([int2str(freqs(t)),'Hz '])
        axis([50 5100 0 1])
    end
end


%
% for f=1:6
%     figure(f)
%     set(gcf,'position',[1 31 1600 794])
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.fig']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.bmp']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_',int2str(spatials(f)),'um_normPeaks_version',int2str(version),'.emf']))
%
% end
% for f=1:4
%     figure(f+6)
%     set(gcf,'position',[1 31 1600 794])
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.fig']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.bmp']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_',int2str(freqs(f)),'Hz_normPeaks_version',int2str(version),'.emf']))
%
% end
%
% for d=1:3
%     figure(100*d)
%     set(gcf,'position',[1 31 1600 794])
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_temporal_',Dtype(d),'_normPeaks_version',int2str(version),'.emf']))
%
%     figure(1000*d)
%     set(gcf,'position',[1 31 1600 794])
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_spatial_',Dtype(d),'_normPeaks_version',int2str(version),'.emf']))
% end
% for d=1:3
%     figure(10000*d)
%     set(gcf,'position',[1 31 1600 794])
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_heatmap_',Dtype(d),'_normPeaks_version',int2str(version),'.fig']))
%     saveas(gcf,fullfile('E:\all_species_temp\newAttempt\all_stat',['stat_F1_heatmap_',Dtype(d),'_normPeaks_version',int2str(version),'.bmp']))
% end
%% ---- here starts real bootstrapping ---
%% bootstrap chirp => what makes sense??

clear
close all
nfft=6;
nlist=[2000 4000 6000 8000 10000 15000];
superpath='f:\all_species_temp';
load(fullfile(superpath,'newAttempt','chirp',['chirpRatios_',int2str(nlist(nfft))]))

step=5;
data=colldata2;
fr=fr2;
dataset='smooth'; %data2
% dataset='noSmooth_noNorm'; %data1a
% dataset='noSmooth_Norm'; %data1b

% within species
within=zeros(3,1000);
for t=1:3
    datanow=data{t};
    for it=1:1000
        topick1=randi(size(datanow,1),50,1);
        d1=datanow(topick1,:);
        topick2=randi(size(datanow,1),50,1);
        d2=datanow(topick2,:);
        
    p=ranksum(mean(d1,2),mean(d2,2));
   within(t,it)=p;
    end
end


% across species
across=zeros(3,1000);
cnt=0;
for t=1:2
    for tt=t+1:3
        cnt=cnt+1;
    datanow1=data{t};
    datanow2=data{tt};
    for it=1:1000
        topick1=randi(size(datanow1,1),50,1);
        d1=datanow1(topick1,:);
        topick2=randi(size(datanow2,1),50,1);
        d2=datanow2(topick2,:);
        
    p=ranksum(mean(d1,2),mean(d2,2));
   across(cnt,it)=p;
    end
    end
end
%% bootstrap speed pref.bar
% result: m-p (0.038), m-h (0.003), p-h (0.381)
% result: m (0.507), p (0.502), h (0.511)
clear
close all
superpath='f:\all_species_temp';

types='wph';
cols='bgr';
titles={'mouse','pig','human'};

speedpref=[]; realIDs=[]; tID=[]; dates=[];
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
    
    date=info(:,1);
    allval=zeros(size(info,1),2);
    black=cell2mat(info(goodB,24));
    white=cell2mat(info(goodW,25));
    allval(goodB,1)=black;
    allval(goodW,2)=white;
    allvalorig=max(allval,[],2);
    tmp=find(allvalorig>0);
    date=date(tmp);
    dates=[dates;date];
    allvalorig=allvalorig(tmp);
    realIDs=[realIDs;tmp];
    tID=[tID;repmat(t,length(tmp),1)];
    speedpref=[speedpref;allvalorig];
end
    
save(fullfile(superpath,'newAttempt','speed_bar','speedBar'),'speedpref','tID','realIDs')

dlength=[];
for t=1:3
    tpm=find(tID==t);
    
    dlength=[dlength;length(tpm)];
end
    
% samplesize=round(length(tID)/length(unique(dates)));
samplesize=25;

% within species
within=zeros(3,1000);
for t=1:3
    datanow=find(tID==t);
    datanow=speedpref(datanow);
    for it=1:1000
        topick1=randi(size(datanow,1),samplesize,1);
        d1=datanow(topick1,:);
        topick2=randi(size(datanow,1),samplesize,1);
        d2=datanow(topick2,:);
        
    p=ranksum(mean(d1,2),mean(d2,2));
   within(t,it)=p;
    end
end


% across species
across=zeros(3,1000);
cnt=0;
for t=1:2
    for tt=t+1:3
        cnt=cnt+1;
     datanow1=find(tID==t);
    datanow1=speedpref(datanow1);
     datanow2=find(tID==tt);
    datanow2=speedpref(datanow2);
    for it=1:1000
        topick1=randi(size(datanow1,1),samplesize,1);
        d1=datanow1(topick1,:);
        topick2=randi(size(datanow2,1),samplesize,1);
        d2=datanow2(topick2,:);
        
    p=ranksum(mean(d1,2),mean(d2,2));
   across(cnt,it)=p;
    end
    end
end

clc
median(across,2)
median(within,2)
