

%% this is the one

clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
path1='R:\Users\Katja\human_retina\newAttemptOnlyHum\DG_PCanalysis';

file='clust_Kmeans_PC_F1';
load(fullfile(path1,file))
sel=1; clust=12; %actual clust
gau=1;
st=0;
fact=2;

if sel==1
    data=ID(:,clust-1);
end

superpath='R:\Users\Katja\human_retina';
Dtype='wph';
cols='bgry';
normDG=[];
for ty=3
    load(fullfile(superpath,'newAttempt','DG',[Dtype(ty),'_normPeaks']))
    normDG=[normDG;normPeaks];
end

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0.2 0.2 0.2; 0 1 1; 1 1 0;0.5 0 0];

clu=unique(data);
for clus=1:length(clu)
    currID=find(data==clu(clus));
    
    if gau==1
        
        for t=1:4
            figure(t)
            subplot(4,4,clus)
            if t==4
                tst=1:length(currID);
            else
                tst=find(tID(currID)==t);
            end
            
            if ~isempty(tst)
                curr=normDG(currID(tst),:,2);
                curr=mean(curr,1);
                curr=flipud(reshape(curr,4,6));
                [x,y] = meshgrid(1:6,1:4);
                imagesc(x(:)',y(:)',curr);
                hold on
                %                 colormap('gray')
                title([int2str(clus),' (mean; ',int2str(length(tst)),'cells)'])
                
                
                if t<4
                    curr=normDG(currID(tst),:,2);
                    curr=mean(curr,1);
                    dat=reshape(curr,4,6);
                    dat=flipud(fliplr(dat));
                    %                 dat=fliplr(dat);
                    
                    minData = mean(mean(dat))+st*std(reshape(dat,24,1));
                    p = (dat - minData);
                    tmp=find(p<0);
                    p(tmp)=0;
                    sumData = sum(sum(p));
                    p=p/sumData;
                    
                    
                    [x,y] = meshgrid(6:-1:1,1:4);
                    %                 [x,y] = meshgrid(1:6,1:4);
                    xx=reshape(x,24,1);
                    yy=reshape(y,24,1);
                    pp=reshape(p,24,1);
                    mx=dot(xx,pp); %centerpoint
                    my=dot(yy,pp); %centerpoint
                    
                    c=dot((xx-mx).^2,pp);
                    b=dot((xx-mx).*(yy-my),pp);
                    a=dot((yy-my).^2,pp);
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
                    %                 AngD=rad2deg(Ang); %angle in degrees
                    AngD=Ang*(180/pi);
                    
                    ymin=get(gca,'ylim');
                    ymin=ymin(1);
                    
                    centerX=mx;%x-axis center;1=left
                    centerY=my;%y-axis center; 1=bottom
                    
                    
                    
                    centerX2=mx;%x-axis center;1=left
                    centerY2=7-my;%y-axis center; 1=bottom
                    Ang2=-Ang;
                    %                     [Y,X]=meshgrid(0.5:0.01:4.5,0.5:0.01:6.5);
                    [Y,X]=meshgrid(0:0.01:5,0:0.01:7);
                    
                    theta = pi/2-Ang;
                    %
                    if ~isempty(find(D==0))
                        tmp=find(D==0)
                        D(tmp)=0.1;
                    end
                    
                    
                    aa = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
                    bb = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
                    %                 bb = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
                    %                 aa = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
                    cc = sin(theta)^2/2/(D(1)/1.5)^2 + cos(theta)^2/2/(D(2)/1.5)^2;
                    allClust = exp( - (aa*(Y-centerY).^2 + 2*bb*(Y-centerY).*(X-centerX2) + cc*(X-centerX2).^2)) ;
                    
                    contour(X,Y,allClust(:,:),1,'color','r')
                    hold on
                    plot(centerX,centerY,'or','markerfacecolor','r','markersize',2)
                end
                
            end
            if t==1
                figure(100)
                subplot(4,4,clus)
                [x,y] = meshgrid(0:7,0:5);
                imagesc(x(:)',y(:)',zeros(6,4));
                colormap('gray')
                hold on
            end
            
            
            if ~isempty(tst)
                figure(100)
                subplot(4,4,clus)
                contour(X,Y,allClust(:,:),1,'color','r')
                hold on
                plot(centerX,centerY,'or','markerfacecolor','r','markersize',2)
                plot([0.5 0.5],[0.5 4.5],'-w')
                plot([0.5 6.5],[0.5 0.5],'-w')
                plot([6.5 6.5],[0.5 4.5],'-w')
                plot([0.5 6.5],[4.5 4.5],'-w')
            end
        end
    end
    
    
end
figure
curr=normDG(:,:,2);
curr=mean(curr,1);
curr=flipud(reshape(curr,4,6));
[x,y] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',curr);
%% this is the one: log

clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
path1='R:\Users\Katja\human_retina\newAttemptOnlyHum\DG_PCanalysis';

file='clust_Kmeans_PC_F1';
load(fullfile(path1,file))
sel=1; clust=12; %actual clust
gau=1;
st=0;
fact=2;

if sel==1
    data=ID(:,clust-1);
end

superpath='R:\Users\Katja\human_retina';
Dtype='wph';
cols='bgry';
normDG=[];
for ty=3
    load(fullfile(superpath,'newAttempt','DG',[Dtype(ty),'_normPeaks']))
    normDG=[normDG;normPeaks];
end

colors=[1 0 0; 0 1 0; 0 0 1;0.5 0 0.5;0 0.5 0.5; 0.5 0.5 0.5; 1 0.5 0; 0 1 0.5; 1 0 1; 0.2 0.2 0.2; 0 1 1; 1 1 0;0.5 0 0];

clu=unique(data);
for clus=1:length(clu)
    currID=find(data==clu(clus));
    
    if gau==1
        
        for t=1:4
            figure(t)
            subplot(4,4,clus)
            if t==4
                tst=1:length(currID);
            else
                tst=find(tID(currID)==t);
            end
            
            if ~isempty(tst)
                curr=normDG(currID(tst),:,2);
                tmp=find(curr==0);
                curr=log10(curr*1000);
                curr(tmp)=0;
                curr=mean(curr,1);
                curr=curr/max(curr);
                curr=flipud(reshape(curr,4,6));
                [x,y] = meshgrid(1:6,1:4);
                imagesc(x(:)',y(:)',fliplr(curr));
                hold on
                %                 colormap('gray')
                title([int2str(clus),' (mean; ',int2str(length(tst)),'cells)'])
                
                
                if t<4
                    curr=normDG(currID(tst),:,2);
                    ctmp=find(curr==0);
                    curr=log10(curr*1000);
                    curr(tmp)=0;
                    curr=mean(curr,1);
                    curr=curr/max(curr);
                    dat=reshape(curr,4,6);
                    dat=flipud(fliplr(dat));
                    %                 dat=fliplr(dat);
                    
                    minData = mean(mean(dat))+st*std(reshape(dat,24,1));
                    p = (dat - minData);
                    tmp=find(p<0);
                    p(tmp)=0;
                    sumData = sum(sum(p));
                    p=p/sumData;
                    
                    
                    [x,y] = meshgrid(6:-1:1,1:4);
                    %                 [x,y] = meshgrid(1:6,1:4);
                    xx=reshape(x,24,1);
                    yy=reshape(y,24,1);
                    pp=reshape(p,24,1);
                    mx=dot(xx,pp); %centerpoint
                    my=dot(yy,pp); %centerpoint
                    
                    c=dot((xx-mx).^2,pp);
                    b=dot((xx-mx).*(yy-my),pp);
                    a=dot((yy-my).^2,pp);
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
                    %                 AngD=rad2deg(Ang); %angle in degrees
                    AngD=Ang*(180/pi);
                    
                    ymin=get(gca,'ylim');
                    ymin=ymin(1);
                    
                    centerX=mx;%x-axis center;1=left
                    centerY=my;%y-axis center; 1=bottom
                    
                    
                    
                    centerX2=mx;%x-axis center;1=left
                    centerY2=7-my;%y-axis center; 1=bottom
                    Ang2=-Ang;
                    %                     [Y,X]=meshgrid(0.5:0.01:4.5,0.5:0.01:6.5);
                    [Y,X]=meshgrid(0:0.01:5,0:0.01:7);
                    
                    theta = pi/2-Ang;
                    %
                    if ~isempty(find(D==0))
                        tmp=find(D==0)
                        D(tmp)=0.1;
                    end
                    
                    
                    aa = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
                    bb = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
                    %                 bb = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
                    %                 aa = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
                    cc = sin(theta)^2/2/(D(1)/1.5)^2 + cos(theta)^2/2/(D(2)/1.5)^2;
                    allClust = exp( - (aa*(Y-centerY).^2 + 2*bb*(Y-centerY).*(X-centerX2) + cc*(X-centerX2).^2)) ;
                    
                    contour(X,Y,flipud(allClust(:,:)),1,'color','r')
                    hold on
                    plot(7-centerX,centerY,'or','markerfacecolor','r','markersize',2)
                end
                
            end
            if t==1
                figure(100)
                subplot(4,4,clus)
                [x,y] = meshgrid(0:7,0:5);
                imagesc(x(:)',y(:)',zeros(6,4));
                colormap('gray')
                hold on
            end
            
            
            if ~isempty(tst)
                figure(100)
                subplot(4,4,clus)
                contour(X,Y,flipud(allClust(:,:)),1,'color','r')
                hold on
                plot(7-centerX,centerY,'or','markerfacecolor','r','markersize',2)
                plot([0.5 0.5],[0.5 4.5],'-w')
                plot([0.5 6.5],[0.5 0.5],'-w')
                plot([6.5 6.5],[0.5 4.5],'-w')
                plot([0.5 6.5],[4.5 4.5],'-w')
            end
        end
    end
    
    
end
figure
curr=normDG(:,:,2);
curr=mean(curr,1);
tmp=find(curr==0);
curr=log10(curr*1000);
curr(tmp)=0;
curr=curr/max(curr);
curr=flipud(reshape(curr,4,6));
[x,y] = meshgrid(1:6,1:4);
imagesc(x(:)',y(:)',fliplr(curr));
%% !!plot all heatmaps per cluster
clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
% file='clust_MixGau_GaussianNew2';
path1='R:\Users\Katja\human_retina\newAttemptOnlyHum\DG_PCanalysis';
file='clust_Kmeans_amp_F1';
load(fullfile(path1,file))
sel=1; clust=10;


if sel==1
    data=ID(:,clust-1);
end


superpath='R:\Users\Katja\human_retina';
Dtype='wph';
normDG=[];
for ty=3
    load(fullfile(superpath,'newAttempt','DG',[Dtype(ty),'_normPeaks']))
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
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
        
        curr=flipud(reshape(curr,4,6));
        [y,x] = meshgrid(1:6,1:4);
        imagesc(x(:)',y(:)',curr,[0 1]);
        colormap('jet')
        title(['clust ',int2str(c),' / ',int2str(currID(i))])
    end
end

%% plot mean heatmaps for cluster result
clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
path1='R:\Users\Katja\human_retina\newAttemptOnlyHum\DG_PCanalysis';
file='clust_Kmeans_amp_F1';
load(fullfile(path1,file))
sel=1; clust=12; %actual clust
gau=0;

if sel==1
    data=ID(:,clust-1);
end

superpath='R:\Users\Katja\human_retina';
Dtype='wph';
normDG=[];
for ty=3
    load(fullfile(superpath,'newAttempt','DG',[Dtype(ty),'_normPeaks']))
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
%% !!new silhouette attempt
clear
close all
plotting=1;
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
path1='R:\Users\Katja\human_retina\newAttemptOnlyHum\DG_PCanalysis'; %newAttempt3: corrected for outlier
% type='GaussianNew2';
type='amp_F1';
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
%% spatial tuning

% for 2Hz
id=2:4:24;
data=normPeaks(:,id,2);

figure
coll1=[];
coll2=[];
cnt=0;
for i=6:-1:1
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 2Hz')
end
subplot(2,1,1)
loglog([4000 2000 1000 500 200 100],coll1)
hold on
title('mean 2Hz')
axis([100 10000 0.1 1])

subplot(2,1,2)
loglog([4000 2000 1000 500 200 100],coll2)
hold on
title('median 2Hz')
axis([100 10000 0.1 1])

% for 4Hz
id=3:4:24;
data=normPeaks(:,id,2);


coll1=[];
coll2=[];
cnt=0;
for i=1:6
    cnt=cnt+1;
    %     subplot(2,2,3)
    %     semilogx(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 4Hz')
    
    % subplot(2,2,4)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 4Hz')
end
subplot(2,1,1)
loglog([100 200 500 1000 2000 4000],coll1)
hold on
title('mean 4Hz')
%  axis([100 10000 0.1 1])
subplot(2,1,2)
loglog([100 200 500 1000 2000 4000],coll2)
hold on
title('median 4Hz')
%  axis([100 10000 0.1 1])
%% spatial tuning: max of all and mean of all
% for 2Hz
id=2:4:24;
data=normPeaks(:,id,2);

figure
coll1=[];
coll2=[];
cnt=0;
for i=6:-1:1
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 2Hz')
end
subplot(2,2,1)
loglog([4000 2000 1000 500 200 100],coll1)
hold on
title('mean 2Hz')
axis([100 10000 0.1 1])

subplot(2,2,3)
loglog([4000 2000 1000 500 200 100],coll1)
hold on
title('mean 2Hz')
axis([100 10000 0.1 1])


% for 4Hz
id=3:4:24;
data=normPeaks(:,id,2);


coll1=[];
coll2=[];
cnt=0;
for i=1:6
    cnt=cnt+1;
    %     subplot(2,2,3)
    %     semilogx(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 4Hz')
    
    % subplot(2,2,4)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 4Hz')
end
subplot(2,2,2)
loglog([100 200 500 1000 2000 4000],coll1)
hold on
title('mean 4Hz')
%  axis([100 10000 0.1 1])
subplot(2,2,4)
loglog([100 200 500 1000 2000 4000],coll2)
hold on
title('median 4Hz')
%  axis([100 10000 0.1 1])


data=normPeaks(:,:,2);

% figure
coll1=[];
coll2=[];
cnt=0;
for i=6:-1:1
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    curr=data(:,i*4-3:i*4);
    ma=max(curr,[],2);
    coll1=[coll1 mean(ma)];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    me=mean(curr,2);
    coll2=[coll2 mean(me)];
    % title('median 2Hz')
end
subplot(2,2,1)
loglog([4000 2000 1000 500 200 100],coll1,'-r')
title('max')
axis([100 10000 0.1 1])
subplot(2,2,2)
loglog([4000 2000 1000 500 200 100],coll1,'-r')
title('max')
axis([100 10000 0.1 1])

subplot(2,2,4)
loglog([4000 2000 1000 500 200 100],coll2,'-r')
title('median')
axis([100 10000 0.1 1])
subplot(2,2,3)
loglog([4000 2000 1000 500 200 100],coll2,'-r')
title('median')
axis([100 10000 0.1 1])

%% spatial LOG TO TAKE
data=normPeaks(:,:,2);

% figure
coll1=[];
coll2=[];
coll3=[];
cnt=0;
for i=6:-1:1
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    curr=data(:,i*4-3:i*4);
    ma=max(curr,[],2);
    coll1=[coll1 median(ma)];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    me=median(curr,2);
    coll2=[coll2 median(me)];
    % title('median 2Hz')
    tmp=find(ma>0);
    coll3=[coll3 median(ma(tmp))];
end
figure
loglog([0.075 0.15 0.3 0.6 1.5 3],coll1,'-r')
title('median')
axis([0.01 1 0.1 1])

figure
loglog([0.075 0.15 0.3 0.6 1.5 3],coll3,'-r')
title('median non-zero')
axis([0.01 1 0.1 1])
%% temporal tuning

% for 2000um
id=17:20;
data=normPeaks(:,id,2);

figure
coll1=[];
coll2=[];
cnt=0;
for i=1:4
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 2Hz')
end
subplot(2,2,1)
loglog([1 2 4 8],coll1)
title('mean 2000um')
axis([1 8 0.01 1])

subplot(2,2,2)
loglog([1 2 4 8],coll2)
title('median 2000um')
axis([1 8 0.01 1])


% for 4000um
id=21:24;
data=normPeaks(:,id,2);


coll1=[];
coll2=[];
cnt=0;
for i=1:4
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 2Hz')
end
subplot(2,2,3)
loglog([1 2 4 8],coll1)
title('mean 4000um')
axis([1 8 0.01 1])

subplot(2,2,4)
loglog([1 2 4 8],coll2)
title('median 4000um')
axis([1 8 0.01 1])
%% temporal tuning: mean and max

id=17:20;
data=normPeaks(:,id,2);

figure
coll1=[];
coll2=[];
cnt=0;
for i=1:4
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 2Hz')
end
subplot(2,2,1)
loglog([1 2 4 8],coll1)
hold on
title('mean 2000um')
axis([1 8 0.01 1])

subplot(2,2,3)
loglog([1 2 4 8],coll1)
hold on
title('mean 2000um')
axis([1 8 0.01 1])


% for 4000um
id=21:24;
data=normPeaks(:,id,2);


coll1=[];
coll2=[];
cnt=0;
for i=1:4
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 2Hz')
end
subplot(2,2,2)
loglog([1 2 4 8],coll1)
hold on
title('mean 4000um')
axis([1 8 0.01 1])

subplot(2,2,4)
loglog([1 2 4 8],coll1)
hold on
title('mean 4000um')
axis([1 8 0.01 1])



data=normPeaks(:,:,2);

coll1=[];
coll2=[];
cnt=0;
for i=1:4
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    curr=data(:,i:4:24);
    ma=max(curr,[],2);
    coll1=[coll1 mean(ma)];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    me=mean(curr,2);
    coll2=[coll2 mean(me)];
    % title('median 2Hz')
end
subplot(2,2,1)
loglog([1 2 4 8],coll1,'-r')
title('max')
axis([1 8 0.01 1])

subplot(2,2,2)
loglog([1 2 4 8],coll1,'-r')
title('max')
axis([1 8 0.01 1])

subplot(2,2,3)
loglog([1 2 4 8],coll2,'-r')
title('median')
axis([1 8 0.01 1])
subplot(2,2,4)
loglog([1 2 4 8],coll2,'-r')
title('median')
axis([1 8 0.01 1])
%% temporal LOG TO TAKE MAX
data=normPeaks(:,:,2);

coll1=[];
coll2=[];
coll3=[];
cnt=0;
for i=1:4
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    curr=data(:,i:4:24);
    ma=max(curr,[],2);
    coll1=[coll1 median(ma)];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    me=mean(curr,2);
    coll2=[coll2 median(me)];
    % title('median 2Hz')
    
    tmp=find(ma>0);
    coll3=[coll3 median(ma(tmp))];
end
figure
loglog([1 2 4 8],coll1,'-r')
title('median')
axis([1 10 0.1 1])

figure
loglog([1 2 4 8],coll3,'-r')
title('median non-zero')
axis([1 10 0.1 1])

%% temporal tuning2

% for 2000um
id=17:20;
data=normPeaks(:,id,2);

figure
coll1=[];
coll2=[];
cnt=0;
for i=1:4
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 2Hz')
end
subplot(2,2,1)
semilogx([1 2 4 8],coll1)
title('mean 2000um')
axis([1 8 0 1])

subplot(2,2,2)
semilogx([1 2 4 8],coll2)
title('median 2000um')
axis([1 8 0 1])


% for 4000um
id=21:24;
data=normPeaks(:,id,2);


coll1=[];
coll2=[];
cnt=0;
for i=1:4
    cnt=cnt+1;
    %     semilogx(2,2,1)
    %     plot(cnt,mean(data(:,i)),'o')
    %     hold on
    coll1=[coll1 mean(data(:,i))];
    %     title('mean 2Hz')
    
    % subplot(2,2,2)
    %     semilogx(cnt,median(data(:,i)),'o')
    %     hold on
    coll2=[coll2 median(data(:,i))];
    % title('median 2Hz')
end
subplot(2,2,3)
semilogx([1 2 4 8],coll1)
title('mean 4000um')
axis([1 8 0 1])

subplot(2,2,4)
semilogx([1 2 4 8],coll2)
title('median 4000um')
axis([1 8 0 1])

%% plot scale bar heatmaps
figure
imagesc(0:0.01:1)
colormap gray
%% CORR: spatial and temporalLOG TO TAKE
data=normPeaks(:,:,2);

% figure
collT=[];
collS=[];

cnt=0;

close all
for i=1:length(togo)
    curr=data(i,:);
    curr2=reshape(curr,4,6);
    suSpat=sum(curr2,1);
    suTemp=sum(curr2,2);
    [~,maSpat]=max(suSpat);%1=100um
    [~,maTemp]=max(suTemp);%1=8Hz
    
    figure(1)
    loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(curr2(maTemp,:)/max(curr2(maTemp,:))),'-b')
    hold on
    
    figure(2)
    loglog([1 2 4 8],curr2(:,maSpat)'/max(curr2(:,maSpat)),'-b')
    hold on
    
    collS=[collS;curr2(maTemp,:)];
    collT=[collT;curr2(:,maSpat)'];
end

figure(1)
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(mean(collS)/max(mean(collS))),'-r')
figure(2)
loglog([1 2 4 8],mean(collT)/max(mean(collT)),'-r')

figure
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(median(collS)),'-r')
title('median spatial')
axis([0.01 4 0.1 1])

figure
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(mean(collS)),'-r')
title('mean spatial')
axis([0.01 4 0.1 1])

med=[];
for j=1:size(collS,2)
    tmp=find(collS(:,j)>0);
    med=[med median(collS(tmp,j))];
end
figure
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(med),'-r')
title('median non-zero spatial')
axis([0.01 1 0.1 1])

figure
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(median(collS)),'-k')
hold on
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(mean(collS)),'-b')
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(med),'--k')
title('spatial')
axis([0.01 1 0.1 1])



figure
loglog([1 2 4 8],median(collT),'-r')
title('median temporal')
axis([1 10 0.1 1])

figure
loglog([1 2 4 8],mean(collT),'-r')
title('mean temporal')
axis([1 10 0.1 1])

med=[];
for j=1:size(collT,2)
    tmp=find(collT(:,j)>0);
    med=[med median(collT(tmp,j))];
end
figure
loglog([1 2 4 8],med,'-r')
title('median non-zero temporal')
axis([1 10 0.1 1])



figure
loglog([1 2 4 8],median(collT),'-k')
hold on
loglog([1 2 4 8],mean(collT),'-b')
loglog([1 2 4 8],med,'--k')
title('temporal')
axis([1 10 0.1 1])
%% POP HEATMAP

data=normPeaks(:,:,2);


% heatMean=mean(data)/max(mean(data));
heatMean=mean(data);
heatMeanR=reshape(heatMean,4,6);
heatMedian=median(data)/max(median(data));
heatMedianR=reshape(heatMedian,4,6);

med=[];
for j=1:size(data,2)
    tmp=find(data(:,j)>0);
    med=[med median(data(tmp,j))];
end
med=med/max(med);
heatMedianNZ=reshape(med',4,6);
figure('color','w')
imagesc(flipud(fliplr(heatMeanR)))
colormap gray
set(gca,'box','off','ycolor','w','xcolor','w','plotboxaspectratio',[1 1/6*4 1])
% figure
% imagesc(flipud(fliplr(heatMedianNZ)))
colorbar
%% HISTOGRAM TEMP, SPAT, SPEED
% clear
% close all
data=normPeaks(:,:,2);
temp = [1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
spat = [0.075 0.075 0.075 0.075 0.15 0.15 0.15 0.15 0.3  0.3 0.3 0.3 0.6 0.6 0.6 0.6 1.5 1.5 1.5 1.5  3 3 3 3 ];
spatum = [4000 4000 4000 4000 2000 2000 2000 2000 1000 1000 1000 1000 500 500 500 500 200 200 200 200 100 100 100 100];

spat = fliplr(spat);
 uspat = [0.075 0.15 0.3 0.6 1.5 3];
 utemp = [1 2 4 8];
 speeds = spatum/1000.*temp;
 uspeeds = unique(speeds);
 
[m1,m2] = max(data');

mspat = spat(m2);
hmspat = [];
for u=1:length(uspat)
   tmp = find(mspat == uspat(u));
   hmspat = [hmspat;length(tmp)]; 
end

figure('color','w')
bar(hmspat)
set(gca,'xtick',1:6,'xticklabel',[0.075 0.15 0.3 0.6 1.5 3])
set(gca,'box','off','fontsize',10)
xtickangle(30)
%-----------
mtemp = temp(m2);
hmtemp = [];
for u=1:length(utemp)
   tmp = find(mtemp == utemp(u));
   hmtemp = [hmtemp;length(tmp)]; 
end

figure('color','w')
bar(hmtemp)
set(gca,'xtick',1:4,'xticklabel',[1 2 4 8])
set(gca,'box','off','fontsize',10)

%------------
mspeeds = speeds(m2);
hmspeeds = [];
for u=1:length(uspeeds)
   tmp = find(mspeeds == uspeeds(u));
   hmspeeds = [hmspeeds;length(tmp)]; 
end

figure('color','w')
bar(hmspeeds)
set(gca,'xtick',1:length(uspeeds),'xticklabel',uspeeds)
% xtickformat('%.1f')

set(gca,'box','off','fontsize',10)
xtickangle(30)

%% spatial and temporal based on average heatmap

data=normPeaks(:,:,2);


heatMean=mean(data);
heatMeanR=reshape(heatMean,4,6);
heatMedian=median(data)/max(median(data));
heatMedianR=reshape(heatMedian,4,6);

med=[];
for j=1:size(data,2)
    tmp=find(data(:,j)>0);
    med=[med median(data(tmp,j))];
end
med=med/max(med);
heatMedianNZ=reshape(med',4,6);
figure
imagesc(flipud(fliplr(heatMeanR)))
colorbar
% saveas(gcf,fullfile('R:\Users\Katja\human_retina\Figs','heatmap_corrMean.emf'))

curr2=heatMeanR;
suSpat=sum(curr2,1);
suTemp=sum(curr2,2);
[~,maSpat]=max(suSpat);%1=100um
[~,maTemp]=max(suTemp);%1=8Hz

%version 1
% spat=curr2(maTemp,:);
% temp=curr2(:,maSpat)';

%version 2
spat=max(curr2);
temp=max(curr2');

close all
figure(1)
loglog([0.07 0.13 0.27 0.53 1.33 2.66],fliplr(spat),'-b')
hold on

figure(2)
loglog([1 2 4 8],temp,'-b')
hold on

curr2=heatMedianR;
suSpat=sum(curr2,1);
suTemp=sum(curr2,2);
[~,maSpat]=max(suSpat);%1=100um
[~,maTemp]=max(suTemp);%1=8Hz

spat=curr2(maTemp,:);
temp=curr2(:,maSpat)';


figure(1)
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(spat),'-k')
figure(2)
loglog([1 2 4 8],temp,'-k')

curr2=heatMedianNZ;
suSpat=sum(curr2,1);
suTemp=sum(curr2,2);
[~,maSpat]=max(suSpat);%1=100um
[~,maTemp]=max(suTemp);%1=8Hz

spat=curr2(maTemp,:);
temp=curr2(:,maSpat)';

figure(1)
loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(spat),'--k')
title('spatial')
loglog([0.075 3],[0.1 0.1],'-r')
axis([0.01 5 0.06 1])
% saveas(gcf,fullfile('R:\Users\Katja\human_retina\Figs','spatial_newCorr.emf'))
figure(2)
loglog([1 2 4 8],temp,'--k')
title('temporal')
axis([1 10 0.1 1])
% saveas(gcf,fullfile('R:\Users\Katja\human_retina\Figs','temporal_newCorr.emf'))



