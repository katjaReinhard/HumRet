%%
clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
load(fullfile(path1,'h_responding'))

%%
load(fullfile(path1,['h_info']))
load(fullfile(path1,['h_normPeaks']))
rDates=info(togo,1);
rUnits=info(togo,2);

cd('R:\Users\Katja\human_retina\_scripts')

dglist={'0100um%dur_12000_speed_1';'0100um%dur_12000_speed_2';'0100um%dur_12000_speed_4';'0100um%dur_12000_speed_8';'0200um%dur_12000_speed_1';'0200um%dur_12000_speed_2';'0200um%dur_12000_speed_4';'0200um%dur_12000_speed_8';'0500um%dur_12000_speed_1';'0500um%dur_12000_speed_2';'0500um%dur_12000_speed_4';'0500um%dur_12000_speed_8';'1000um%dur_12000_speed_1%';'1000um%dur_12000_speed_2';'1000um%dur_12000_speed_4';'1000um%dur_12000_speed_8';'2000um%dur_12000_speed_1';'2000um%dur_12000_speed_2';'2000um%dur_12000_speed_4';'2000um%dur_12000_speed_8';'4000um%dur_12000_speed_1';'4000um%dur_12000_speed_2';'4000um%dur_12000_speed_4';'4000um%dur_12000_speed_8'};
spats=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
frs=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];

harmonics=zeros(length(togo),20);
for g=1:length(togo)
    
    currUnit=rUnits{g};
    currDate=rDates{g};
    
    %     startF=info{resp(g),4};
    %     stopF=info{resp(g),5};
    startF=info{togo(g),4};
    stopF=info{togo(g),5};
    
    
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
    
    
    [dg id]=max(normPeaks(g,:,2));
    currdg=dglist{id};
    
    heka=unit{1,2}(:,1);
    file=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},currdg))
            file=[file;h];
        end
    end
    %     tmp=find(file>=startF);
    %     file=file(tmp);
    %     tmp=find(file<=stopF);
    %     file=file(tmp);
    
    
    allrate=[];
    for r=1:length(file)
        spikes=unit{1,2}{file(r),2};
        spikes=round(spikes);
        rate=zeros(1,16500);
        rate(spikes)=1;
        rate(spikes(diff(spikes)==0))=2;
        allrate=[allrate;rate(2200:14200)];
    end
    
    if mod(r,2)==0
        endr=r;
    else
        endr=r-1;
    end
    
    cnt=0;
    for j=1:2:endr
        cnt=cnt+1;
        if j==endr
            meanrate=mean(allrate(j:end,:),1);
        else
            meanrate=mean(allrate(j:j+1,:),1);
        end
        
        Y=fft(meanrate,15000)/numel(meanrate);
        freq=1000/2*linspace(0,1,15000/2+1);
        data=2*abs(Y(1:15000/2+1));
        data=data';
        
        tmp1=find(freq>=0.35);
        tmp2=find(freq<=30);
        data=data(tmp1(1):tmp2(end));
        freq=freq(tmp1(1):tmp2(end));
        
%         figure
%         plot(freq,data)
        
        tempData=data;
        excluded=[];
        for f=[0.5 1 2 3];
            tmp=find(freq==frs(id)*f);
            excluded=[excluded tmp-2:tmp+2];
        end
        included=1:length(tempData);
        start=find(freq<=frs(id)*0.5);
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
            tmp=find(freq==frs(id)*testf(f));
            if isempty(tmp)
                tst=find(freq<frs(id)*testf(f));
                tst=tst(end);
                flist=tst:tst+1;
                [~,totake]=max(data(flist));
            else
                flist=tmp-1:tmp+1;
                [~,totake]=max(data(tmp-1:tmp+1));
            end
            tmp=flist(totake);
            aa=find(freq<=frs(id)*testf(f)-frs(id)/2);
            if isempty(aa)
                aa=1;
            end
            aa=aa(end);
            bb=find(freq>=frs(id)*testf(f)+frs(id)/2);
            bb=bb(1);
            
            maxs=[maxs max(data(aa:bb))];
            peaks=[peaks data(tmp)];
            ids=[ids tmp];
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
        
        %         tmp=find(freq==frs(id));
        %         if isempty(tmp)
        %             tst=find(freq<frs(id));
        %             tst=tst(end);
        %             flist=tst:tst+1;
        %             [~,totake]=max(data(flist));
        %         else
        %             flist=tmp-1:tmp+1;
        %             [~,totake]=max(data(tmp-1:tmp+1));
        %         end
        %         tmp=flist(totake);
        %         aa=find(freq<=frs(id)-frs(id)/2);
        %         if isempty(aa)
        %             aa=1;
        %         end
        %         aa=aa(end);
        %         bb=find(freq>=frs(id)+frs(id)/2);
        %         bb=bb(1);
        %
        %         maxs=max(data(aa:bb));
        %         peak=data(tmp);
        %
        %         if peak<maxs
        %             peak=0;
        %         end
        %         if peak<=thr
        %             peak=0;
        %         end
        %
        %
        if peaks(2)>0
            harmonics(g,cnt)=peaks(2);
        else
            harmonics(g,cnt)=max(peaks);
        end
        
    end
    harmonics(g,cnt+1:end)=-1;
end

%%
close all
% figure

normFFT=zeros(length(harmonics),10);
stoppingCells=0;

for h=1:length(harmonics)
    curr=harmonics(h,:);
    tmp=find(curr>-1);
    dat=curr(tmp);
%     plot(dat)
%     hold on
    tmp2=find(curr>0);
    dat2=curr(tmp2);
    dat2=dat2/max(dat2);
    normFFT(h,tmp2)=dat2;
    
    first=0; foundzero=0;
    for j=1:length(tmp)
        if dat(j)>0 && first==0
            first=j;
        end
        if dat(j)==0 && foundzero<j;
            foundzero=j;
        end
    end
    if foundzero>first
        stoppingCells=stoppingCells+1;
    end
    
    
end


figure
for i=1:7
    curr=normFFT(:,i);
    tmp=find(curr>0);
    plot(i,mean(curr(tmp)),'o')
    hold on
    plot([i i],[mean(curr(tmp))-std(curr(tmp)) mean(curr(tmp))+std(curr(tmp))],'-')
end
    
    
    



