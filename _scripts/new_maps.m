%% electrode maps

clear
close all

path1='R:\Users\Katja\human_retina';
types='wph';

FLASH=zeros(3,1);
LF=zeros(3,1);
DS=zeros(3,1);
DG=zeros(3,1);
totCells=zeros(3,1);

clc
for t=3
    
    if t==1
        currpath='S:\data\Katja\WT';
    elseif t==2
        currpath='S:\data\Katja\Pig';
        
    else
        currpath='R:\Users\Katja\human_retina\Human';
    end
    
    
    load(fullfile(path1,[types(t),'_info']))
    goodones=cell2mat(info(:,3));
    goodones=find(goodones>0);
    
    if t==2
        badIDs=202:250;
        bads=[];
        for b=1:length(badIDs)
            tmp=find(goodones==badIDs(b));
            bads=[bads;tmp];
        end
        goodones=setdiff(goodones,bads);
    end
    
    new=info(goodones,:);
    dates=new(:,1);
    udates=unique(dates);
    
    
    
    %flash
    flash=cell2mat(new(:,26));
    ID=find(flash==1);
    goodIDs=ID;
    
    
    %LF
    ID=cell2mat(new(:,9));
    ID=find(ID~=0);
    goodIDs=[goodIDs;ID];
    
    
    
    %DG
    load(fullfile(path1,'newAttempt2','DG',[types(t),'_normPeaks']))
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    % chirp
    load(fullfile(path1,'newAttempt3','chirp',['chirpRatios_',int2str(8000)]))
    togo=togos{t};
    ID=[];
    for to=1:length(togo)
        tmp=find(goodones==togo(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    
    %6vel
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
    
    ID=[];
    for to=1:length(goodB)
        tmp=find(goodones==goodB(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    ID=[];
    for to=1:length(goodW)
        tmp=find(goodones==goodW(to));
        ID=[ID;tmp];
    end
    goodIDs=[goodIDs;ID];
    
    
    goodIDs=unique(goodIDs);
    
    
    
    % now go through all dates
    dates=info(goodones,1);
    units=info(goodones,2);
    udates=unique(dates);
    [udates,ia] =unique(dates);
    start=[1;ia+1];
%     [start,ia]=sort(ia);
%     start=[start;length(dates)+1];
%     udates=udates(ia);
    
    CHmap=[0 24 26 29 32 35 37 0; 21 22 25 30 31 36 39 40; 19 20 23 28 33 38 41 42; 16 17 18 27 34 43 44 45; 15 14 13 4 57 48 47 46;12 11 8 3 58 53 50 49;10 9 6 1 60 55 52 51;0 7 5 2 59 56 54 0];
    CHmap=flipud(CHmap);
    for u=1:length(udates)
        figure
        un=units(start(u):start(u+1)-1);
        upath=fullfile(currpath,udates{u},'units');
        ulist=dir(fullfile(upath,'*mat'));
        CH=[];
        for uu=1:length(un)
            m=0; found=0;
            while found==0
                m=m+1;
                if ~isempty(regexp(ulist(m).name,un{uu}))
                    found=1;
                end
            end
            name=ulist(m).name;
            tmp=regexp(name,'CH');
            tmp2=regexp(name(tmp:end),'_');
            num=str2num(name(tmp+2:tmp2(1)+tmp-2));
            CH=[CH;num];
        end
        
        
        goodies=[];
        for s=start(u):start(u+1)-1
            if ~isempty(find(goodIDs==s))
                goodies=[goodies;1];
            else
                goodies=[goodies;0]; %responding cells
            end
        end
        
        for xx=1:8
            for yy=1:8
                plot(xx,yy,'o','color',[0.6 0.6 0.6],'markerfacecolor',[0.6 0.6 0.6],'markersize',10)
                hold on
            end
        end
        plot(1,1,'ow','markerfacecolor','w','markersize',10)
        plot(8,1,'ow','markerfacecolor','w','markersize',10)
        plot(1,8,'ow','markerfacecolor','w','markersize',10)
        plot(8,8,'ow','markerfacecolor','w','markersize',10)
        
        bad=find(goodies==0);
%         for b=1:length(bad)
%             [xx,yy]=find(CHmap==CH(bad(b)));
%             plot(xx,yy,'o','color','r','markerfacecolor','r','markersize',45)
%         end
        nice=find(goodies==1);
        for n=1:length(nice)
            [yy,xx]=find(CHmap==CH(nice(n)));
            plot(xx,yy,'o','color','g','markerfacecolor','g','markersize',28)
        end
        axis([0 9 0 9])
%         title([types(t),' / ', udates{u}])
%         set(gcf,'position',[1 31 1600 794])
%         saveas(gcf,fullfile('F:\all_species_temp\CHmaps',[types(t),'_',udates{u},'.emf']))
        
    end
    
end