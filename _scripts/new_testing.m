clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
for i=1:30
    curr=normPeaks(i,:,2);
    
    curr2=flipud(reshape(curr,4,6));
    
    
    [X,Y]=meshgrid(0:0.02:5,0:0.02:7);
    theta = -deg2rad(mean(keepOrig(currID,5)));
    aa = cos(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + sin(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
    bb = -sin(2*theta)/4/(mean(keepOrig(currID,3))/1.5)^2 + sin(2*theta)/4/(mean(keepOrig(currID,4))/1.5)^2 ;
    cc = sin(theta)^2/2/(mean(keepOrig(currID,3))/1.5)^2 + cos(theta)^2/2/(mean(keepOrig(currID,4))/1.5)^2;
    allClust = exp( - (aa*(X-mean(keepOrig(currID,1))).^2 + 2*bb*(X-mean(keepOrig(currID,1))).*(Y-mean(keepOrig(currID,2))) + cc*(Y-mean(keepOrig(currID,2))).^2)) ;
    
    contour(X,Y,allClust(:,:),1,'color',[0.9 0.9 0.9])
    
    
    
    %     [x,y]=meshgrid([100 200 500 1000 2000 4000],[ 1 2 4 8]);
    %
    %     curr2=curr2(:);
    %     x=x(:);
    %     y=y(:);
    %
    %     gaussmodel='a1*exp(-((x-x0).^2/(2*sigmax^2)+(y-y0).^2/(2*sigmay^2)))';
    %     fitobject = fit([x,y],curr2,gaussmodel,'start',[1 1500 500 3 1.5]);
    %     a1=fitobject.a1;
    %     sigmax=fitobject.sigmax;
    %     sigmay=fitobject.sigmay;
    %     x0=fitobject.x0;
    %     y0=fitobject.y0;
    %
    %
    %     formula=a1*exp(-((x-x0).^2/(2*sigmax^2)+(y-y0).^2/(2*sigmay^2)));
    %
    %
    %     [x,y]=meshgrid([100 200 500 1000 2000 4000],[ 1 2 4 8]);
    %
    %
    %     figure
    %     imagesc(reshape(curr2,4,6))
    %     %       colormap('gray')
    %     hold on
    %     contour(reshape(formula,4,6))
end
%% output
a1=1.254;
sigmax=2.352;
sigmay=0.9872;
x0=7.32;
y0=2.844;

formula=a1*exp(-((x-x0).^2/(2*sigmax^2)+(y-y0).^2/(2*sigmay^2)));

figure
imagesc(reshape(formula,4,6))

figure
contour(reshape(formula,4,6))

figure
imagesc(reshape(curr2,4,6))
colormap('gray')
hold on
contour(reshape(formula,4,6))

%%
close all
clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
charact=zeros(length(normPeaks),11);

spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];
speeds{1}=[17 22];
speeds{2}=[13 18 23];
speeds{3}=[14 19 24];
speeds{4}=[15 20];

for i=1:length(normPeaks)
    curr=normPeaks(i,:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    
    %     curr2=flipud(reshape(curr,4,6));
    curr2=reshape(curr,4,6);
    %
    %     figure
    %     imagesc(curr2)
    
    %spatial
    %     sums=sum(curr2,2);
    %     [ma id]=max(sums); %this is the temporal freq for which we check
    for id=1:4
        cumu=cumsum(curr2(id,:));
        cumu50=find(cumu>=cumu(end)/2);
        cumu50=cumu50(1);
        charact(i,id)=cumu50;
    end
    
    %     sums(id)=0;
    %     [ma id]=max(sums); %this is the temporal freq for which we check (2nd highest)
    %     cumu=cumsum(curr2(id,:));
    %     cumu50=find(cumu>=cumu(end)/2);
    %     cumu50=cumu50(1);
    %     charact(i,2)=cumu50;
    
    %temporal
    %     sums=sum(curr2,1);
    %     [ma id]=max(sums); %this is the spatial freq for which we check
    for id=1:6
        cumu=cumsum(curr2(:,id));
        cumu50=find(cumu>=cumu(end)/2);
        cumu50=cumu50(1);
        charact(i,id+4)=cumu50;
    end
    
    %     sums(id)=0;
    %     [ma id]=max(sums); %this is the temporal freq for which we check (2nd highest)
    %     cumu=cumsum(curr2(:,id));
    %
    %
    %     cumu50=find(cumu>=cumu(end)/2);
    %     cumu50=cumu50(1);
    %     charact(i,4)=cumu50;
    
    %iso-speed
    sums=[];
    for s=1:length(speeds)
        su=sum(curr2(speeds{s}));
        sums=[sums su];
    end
    st=std(sums);
    charact(i,11)=st;
    
    
    %     title(['spatial: ',int2str(spatial(charact(i,1))),' / ',int2str(spatial(charact(i,2))) '; temporal: ',int2str(temporal(charact(i,3))),' / ',int2str(temporal(charact(i,4))),' / std(iso-speed): ',num2str(st)])
    
end

charact2=zeros(length(normPeaks),3);
for c=1:length(charact)
    curr=charact(c,1:4);
    tmp=find(curr>1);
    m=median(curr(tmp));
    charact2(c,1)=m;
    
    curr=charact(c,5:10);
    m=median(curr);
    charact2(c,2)=m;
    
    charact2(c,3)=charact(c,11);
end
%% pick cells
close all
tocheck=find(charact2(:,2)>=3);
for t=1:length(tocheck)
    curr=normPeaks(tocheck(t),:,2);
    
    %log
    %     tmp=find(curr==0);
    %     curr=log10(curr*1000);
    %     curr(tmp)=0;
    %     curr=mean(curr,1);
    %     curr=curr/max(curr);
    
    
    %     curr2=flipud(reshape(curr,4,6));
    curr2=reshape(curr,4,6);
    %
    figure
    imagesc(curr2)
    title([int2str(charact2(tocheck(t),2)),' / ',num2str(charact2(tocheck(t),3))])
end

%% new Gauss (Xu)

close all
clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
charact=zeros(length(normPeaks),4);

spatial=[0.1 0.2 0.5 1 2 4];
temporal=[1 2 4 8];
sf=spatial;
tf=temporal;




for i=1:30
    curr=normPeaks(i,:,2);
    curr2=reshape(curr,4,6);
    
    
    
    
    f0=curr2;
    f0(f0<=0)=0.0001; % remove negative values
    
    [sfmat,tfmat] = meshgrid(sf,tf);
    [mat_fit,xhat]=fitsftf(f0,sf,tf);
    
    %     pars =  [1,.3,.1,4,2,-.5];
    % %     pars=[1,0.2,0.7,3,1.5,-.3];
    %     A = pars(1); % amplitude
    %     mu_sf = pars(2); % center of SF
    %     sigma_sf = pars(3);  % sigma of SF tuning
    %     mu_tf = pars(4);  % center of TF
    %     sigma_tf = pars(5);   % sigma of TF tuning
    %     Q = pars(6); % -1 to 1, the tan of x/y
    %
    %
    %     tfp = 2.^((Q+1)*( log2(sfmat) - log2(mu_sf) ) + log2(mu_tf));
    %     t1 = exp( -(log2(sfmat)-log2(mu_sf)).^2 / log2(sigma_sf).^2 );
    %     t2 = exp( -(log2(tfmat)-log2(tfp)).^2 / log2(sigma_tf).^2 );
    %
    %     f = A*t1.*t2;
    %
    %      x0 =  [1,.2,.1,4,2,-.5];
    %     lb =  [0,min(sf),.01,min(tf),.1,-1]; % lower boundary
    % ub =  [1000,max(sf),1,max(tf),20,1]; % upper boundary
    % opt = optimset('display','on');
    % xhat = lsqcurvefit(f,x0,sfmat,f0,lb,ub,opt);
    %
    %
    %     [mat_fit,xhat]=fitsftf(curr2,sf,tf);
    %
    formula=xhat(1)*exp(-((tfmat-xhat(4)).^2/(2*xhat(5)^2)+(sfmat-xhat(2)).^2/(2*xhat(3)^2)));
    
    
    
end
%%
close all
clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
charact=zeros(length(normPeaks),3);
vars=zeros(length(normPeaks),11);

spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];
speeds{1}=[17 22];
speeds{2}=[13 18 23];
speeds{3}=[14 19 24];
speeds{4}=[15 20];

for i=1:length(normPeaks)
    curr=normPeaks(i,:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    
    curr2=flipud(reshape(curr,4,6));
    %
    %         figure
    %         imagesc(curr2)
    %
    %spatial
    sums=sum(curr2,2);
    [ma id]=max(sums); %this is the temporal freq for which we check
    % for id=1:4
    cumu=cumsum(curr2(id,:));
    cumu50=find(cumu>=cumu(end)/2);
    cumu50=cumu50(1);
    charact(i,1)=cumu50;
    vars(i,1:6)=max(curr2,[],1);
    % end
    
    %     sums(id)=0;
    %     [ma id]=max(sums); %this is the temporal freq for which we check (2nd highest)
    %     cumu=cumsum(curr2(id,:));
    %     cumu50=find(cumu>=cumu(end)/2);
    %     cumu50=cumu50(1);
    %     charact(i,2)=cumu50;
    
    %temporal
    sums=sum(curr2,1);
    [ma id]=max(sums); %this is the spatial freq for which we check
    % for id=1:6
    cumu=cumsum(curr2(:,id));
    cumu50=find(cumu>=cumu(end)/2);
    cumu50=cumu50(1);
    charact(i,2)=cumu50;
    vars(i,7:10)=max(curr2,[],2);
    % end
    
    %     sums(id)=0;
    %     [ma id]=max(sums); %this is the temporal freq for which we check (2nd highest)
    %     cumu=cumsum(curr2(:,id));
    %
    %
    %     cumu50=find(cumu>=cumu(end)/2);
    %     cumu50=cumu50(1);
    %     charact(i,4)=cumu50;
    
    %iso-speed
    sums=[];
    for s=1:length(speeds)
        su=sum(curr2(speeds{s}));
        sums=[sums su];
    end
    st=std(sums);
    charact(i,3)=st;
    vars(i,11)=st;
    
    %     title(['spatial: ',int2str(spatial(charact(i,1))),' / ',int2str(spatial(charact(i,2))) '; temporal: ',int2str(temporal(charact(i,3))),' / ',int2str(temporal(charact(i,4))),' / std(iso-speed): ',num2str(st)])
    
end

%% pick cells
close all
tocheck=find(charact(:,1)==3);
for t=1:length(tocheck)
    curr=normPeaks(tocheck(t),:,2);
    
    %log
    %     tmp=find(curr==0);
    %     curr=log10(curr*1000);
    %     curr(tmp)=0;
    %     curr=mean(curr,1);
    %     curr=curr/max(curr);
    
    
    %     curr2=flipud(reshape(curr,4,6));
    curr2=reshape(curr,4,6);
    %
    figure
    imagesc(curr2)
    title([int2str(charact(tocheck(t),2)),' / ',num2str(charact(tocheck(t),3))])
end
%% clustering
opts=statset('maxiter',1000);
ID=[];
for k=2:15
    [idx c]=kmeans(vars,k,'emptyaction','drop','start','cluster','replicates',10000,'options',opts);
    k
    ID=[ID idx];
end
%% !!new silhouette attempt


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
            curr2=vars(curr,:);
            
            if u==6 || u==10
                tmp=find(curr2(:,5)<1.5);
                curr2(tmp,5)=pi-curr2(tmp,5);
                vars(curr,5)=curr2(:,5);
            end
            
        end
    end
    
    
    
    
    ais=[]; bis=[];
    for c=1:size(vars,1)
        cpoint=vars(c,:);
        others=1:size(vars,1);
        others=setdiff(others,c)';
        osame=find(data(others)==data(c));
        osame=others(osame);
        odiff=find(data(others)~=data(c));
        odiff=others(odiff);
        
        
        disame=sqrt(sum((repmat(cpoint,length(osame),1)-vars(osame,:)).^2));
        disame=mean(disame);
        didiff=sqrt(sum((repmat(cpoint,length(odiff),1)-vars(odiff,:)).^2));
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

%% !!plot all heatmaps per cluster

close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
% file='clust_MixGau_GaussianNew2';

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
% clear
close all
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% file='clust_Kmeans_GaussianNew2_corr';
% path1='R:\Users\Katja\human_retina\newAttemptOnlyHum\DG_PCanalysis';
% file='clust_Kmeans_amp_F1';
% load(fullfile(path1,file))
sel=1; clust=10; %actual clust
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
%% find parasol cells

close all
% clear
% superpath='R:\Users\Katja\human_retina';
superpath = 'C:\Users\KatjaR\OneDrive - imec\HumRet';
% load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
load(fullfile(superpath,['h','_normPeaks']))

load(fullfile(superpath,'h_info'))
dates=info(togo,1);
charact=zeros(length(normPeaks),3);
vars=zeros(length(normPeaks),11);

spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];

plotit=0;

candidates=[];
for i=1:length(normPeaks)
    curr=normPeaks(i,:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    
    curr2=flipud(reshape(curr,4,6));
    
    spat=max(curr2,[],1);
    spat4=curr2(3,:);
    spat8=curr2(4,:);
    temp=max(curr2,[],2);
    temp_alt=sum(curr2,2)/4;
    temp4000=curr2(:,6);
    
    
    %     if sum(spat(2:end)>0.5)==5 && temp(1)<0.8
    %         if sum(spat4(2:end)>0.5)==5 || sum(spat8(2:end)>0.5)==5
    %             %             if sum(temp(1:4)>0.8)==3 && temp_alt(4)>0.8
    %             if sum(temp(2:4)>0.8)==3 && temp_alt(1)<0.8
    %                 candidates=[candidates; i 1];
    % %                 disp('yes')
    %                 if plotit==1
    %                     figure
    %                     imagesc(curr2)
    %                 end
    %             elseif temp_alt(1)<temp_alt(4) && sum(temp_alt(2:4)>0.8)==3
    %                 candidates=[candidates; i 1];
    % %                 disp('yes')
    %             end
    %         else
    %             if sum(temp(2:4)>0.8)==3
    %                 %                 candidates=[candidates;i 0.5];
    %                 if plotit==1
    %                     figure
    %                     imagesc(curr2)
    %                     title([num2str(candidates(end,2)),' (bcs of spatial)'])
    %                 end
    %             end
    %         end
    
    if sum(spat4(2:end)>0.5)==5 || sum(spat8(2:end)>0.5)==5
        %          if sum(temp(2:4)>0.8)==3 && temp_alt(4)>0.8
        if sum(temp(3:4)>0.8)==2 && temp(1)<0.8
            candidates=[candidates; i 1];
            if plotit==1
                figure
                imagesc(curr2)
                title([num2str(candidates(end,2))])
            end
        elseif sum(temp4000(3:4)>0.8)==2 && temp4000(1)<0.8
            candidates=[candidates; i 1];
            if plotit==1
                figure
                imagesc(curr2)
                title([num2str(candidates(end,2))])
            end
            %                 disp('yes')
        elseif temp_alt(1)<temp_alt(4) && sum(temp_alt(3:4)>0.8)==2
            candidates=[candidates;i 0.5];
            if plotit==1
                figure
                imagesc(curr2)
                title([num2str(candidates(end,2)),', ',num2str(temp_alt(1))])
            end
        elseif temp4000(1)<temp4000(4) && sum(temp4000(3:4)>0.8)==2
            candidates=[candidates;i 0.5];
            if plotit==1
                figure
                imagesc(curr2)
                title([num2str(candidates(end,2)),', ',num2str(temp_alt(1))])
            end
        end
    end
    
end
good=find(candidates(:,2)==1);
% clc
display(['good candidates: ',int2str(length(good)),'; others: ',int2str(length(candidates)-length(good))])

coll=[];
for i=1:length(candidates)
    curr=normPeaks(candidates(i,1),:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    
    curr2=flipud(reshape(curr,4,6));
    spat=max(curr2,[],1);
    coll=[coll;spat];
end
figure
loglog([0.075 0.15 0.3 0.6 1.5 3],coll)
hold on
loglog([0.075 0.15 0.3 0.6 1.5 3],mean(coll),'linewidth',4)
axis([0.01 10 0.1 1])
parasol=candidates(:,1);

%% find midget cells

% close all
% clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
charact=zeros(length(normPeaks),3);
vars=zeros(length(normPeaks),11);

spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];

plotit=0;

candidates=[];
for i=1:length(normPeaks)
    curr=normPeaks(i,:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    
    curr2=flipud(reshape(curr,4,6));
    
    spat=max(curr2,[],1);
    spat4=curr2(3,:);
    spat8=curr2(4,:);
    temp=max(curr2,[],2);
    temp_alt=sum(curr2,2)/4;
    
    
    %     if spat(2)>0.5 && spat(5)<spat(3)
    if spat8(2)>0.5 && spat8(5)<spat8(3) % strong for 200, stronger for 500 than 4000
        %         if sum(spat8(2:end)>0.5)==5
        candidates=[candidates;i 1];
        if plotit==1
            figure
            imagesc(curr2)
            title([num2str(candidates(end,2))])
        end
        %         end
        %         elseif sum(spat8(2:end)>0.5)==5
        %             candidates=[candidates;i 0.5];
        %             if plotit==1
        %                 figure
        %                 imagesc(curr2)
        %                 title([num2str(candidates(end,2))])
        %             end
        %         end
        %     elseif spat8(2)>0.5 && spat8(5)<spat8(3)
        %         if sum(spat8(2:end)>0.5)==5
        %             candidates=[candidates;i 1];
        %             if plotit==1
        %                 figure
        %                 imagesc(curr2)
        %                 title([num2str(candidates(end,2))])
        %             end
        %         end
        %     end
        %
    end
end
good=find(candidates(:,2)==1);
clc
display(['good candidates: ',int2str(length(good)),'; others: ',int2str(length(candidates)-length(good))])


coll=[];
for i=1:length(candidates)
    curr=normPeaks(candidates(i,1),:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    
    curr2=flipud(reshape(curr,4,6));
    spat=max(curr2,[],1);
    coll=[coll;spat];
end
figure
loglog([0.075 0.15 0.3 0.6 1.5 3],coll)
hold on
loglog([0.075 0.15 0.3 0.6 1.5 3],mean(coll),'linewidth',4)
axis([0.01 10 0.1 1])
midget=candidates(:,1);

%% find blue-yellow cells

% close all
% clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))


spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];

plotit=0;

candidates=[];
for i=1:length(normPeaks)
    curr=normPeaks(i,:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    
    curr2=flipud(reshape(curr,4,6));
    
    spat=max(curr2,[],1);
    spat2=curr2(2,:);
    spat4=curr2(3,:);
    spat8=curr2(4,:);
    temp=max(curr2,[],2);
    temp_alt=sum(curr2,2)/4;
    spat_alt=sum(curr2,1)/6;
    temp2000=curr2(:,5);
    temp4000=curr2(:,6);
    
    %     if spat(2)>0.2 && sum(spat(3:6)>0.8)==4
    if spat4(2)<0.2 && sum(spat4(3:6)>0.8)==4
        if sum(temp(2:end)>0.8)==3
            candidates=[candidates;i 1];
            if plotit==1
                figure
                imagesc(curr2)
                title([num2str(candidates(end,2))])
            end
            
            %         elseif sum(temp_alt>0.8)==4
        end
    elseif spat2(2)<0.2 && sum(spat2(3:6)>0.8)==4
        if sum(temp(2:end)>0.8)==3
            candidates=[candidates;i 0.5];
            if plotit==1
                figure
                imagesc(curr2)
                title([num2str(candidates(end,2))])
            end
        end
    end
    
    
end
good=find(candidates(:,2)==1);
clc
display(['good candidates: ',int2str(length(good)),'; others: ',int2str(length(candidates)-length(good))])

coll=[];
for i=1:length(candidates)
    curr=normPeaks(candidates(i,1),:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    
    curr2=flipud(reshape(curr,4,6));
    spat=max(curr2,[],1);
    coll=[coll;spat];
end
figure
loglog([0.075 0.15 0.3 0.6 1.5 3],coll)
hold on
loglog([0.075 0.15 0.3 0.6 1.5 3],mean(coll),'linewidth',4)
axis([0.01 10 0.1 1])
blueyellow=candidates(:,1);
%% find most parasol-like cell per experiment
close all
clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
load(fullfile(superpath,'h_info'))
dates=info(togo,1);
undates=unique(dates);


spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];


plotit=1;

candidates=[];
for i=1:length(undates)
    currD=undates{i};
    
    
    m=0; found=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(dates{m},currD))
            found=1;
        end
    end
    start=m;
    found=0;
    while found==0 && m<length(dates)
        m=m+1;
        if isempty(regexp(dates{m},currD))
            found=1;
        end
    end
    list=start:m-1;
    
    cand=[];
    for j=1:length(list)
        curr=normPeaks(list(j),:,2);
        
        %log
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
        
        
        curr2=flipud(reshape(curr,4,6));
        
        spat=max(curr2,[],1);
        if sum(spat(2:end)>0.5)==5 && sum(spat(4:end)>0.85)==3
            if spat(1)<0.4
                cand=[cand;list(j)];
            end
        end
    end
    
    if plotit==1
        figure
        subs=ceil(length(cand)/2);
        for c=1:length(cand)
            subplot(2,subs,c)
            curr=normPeaks(cand(c),:,2);
            
            %log
            tmp=find(curr==0);
            curr=log10(curr*1000);
            curr(tmp)=0;
            curr=mean(curr,1);
            curr=curr/max(curr);
            
            mas=[];
            for ii=6:-1:1
                dat=curr(:,ii*4-3:ii*4);
                ma=max(dat,[],2);
                mas=[mas;ma];
                
            end
            tmp=find(mas==0);
            mas(tmp)=0.0001;
            %             curr2=flipud(reshape(curr,4,6));
            curr=flipud(reshape(curr,4,6));
            %             imagesc(curr)
            loglog([0.075 0.15 0.3 0.6 1.5 3],mas)
            axis([0.01 10 0.1 1])
            title([currD,' / ',int2str(cand(c))])
        end
        
        if isempty(cand)
            subs=ceil(length(list)/2);
            for c=1:length(list)
                subplot(2,subs,c)
                curr=normPeaks(list(c),:,2);
                
                %log
                tmp=find(curr==0);
                curr=log10(curr*1000);
                curr(tmp)=0;
                curr=mean(curr,1);
                curr=curr/max(curr);
                
                mas=[];
                for ii=6:-1:1
                    dat=curr(:,ii*4-3:ii*4);
                    ma=max(dat,[],2);
                    mas=[mas;ma];
                    
                end
                tmp=find(mas==0);
                mas(tmp)=0.0001;
                %             curr2=flipud(reshape(curr,4,6));
                curr=flipud(reshape(curr,4,6));
                %             imagesc(curr)
                loglog([0.075 0.15 0.3 0.6 1.5 3],mas)
                axis([0.01 10 0.1 1])
                title([currD,' / ',int2str(list(c)),' (all)'])
            end
        end
    else
        figure
        subs=ceil(length(gofor)/2);
        for c=1:length(cand)
            subplot(2,subs,c)
            curr=normPeaks(cand(c),:,2);
            
            %log
            tmp=find(curr==0);
            curr=log10(curr*1000);
            curr(tmp)=0;
            curr=mean(curr,1);
            curr=curr/max(curr);
            
            mas=[];
            for ii=6:-1:1
                dat=curr(:,ii*4-3:ii*4);
                ma=max(dat,[],2);
                mas=[mas;ma];
                
            end
            tmp=find(mas==0);
            mas(tmp)=0.0001;
            %             curr2=flipud(reshape(curr,4,6));
            curr=flipud(reshape(curr,4,6));
            %             imagesc(curr)
            loglog([0.075 0.15 0.3 0.6 1.5 3],mas)
            axis([0.01 10 0.1 1])
            title([currD,' / ',int2str(cand(c))])
        end
        
    end
    
    
    
end
%% ... continued
close all
clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
load(fullfile(superpath,'h_info'))
dates=info(togo,1);
undates=unique(dates);

gofor=[13 33 51 67 84 93 118 136 161 189 220 233 250 256 269];

spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];
cols=colormap;
cols=cols(1:3:end,:);

figure
coll=[];
for c=1:length(gofor)
    curr=normPeaks(gofor(c),:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
    mas=[];
    for ii=6:-1:1
        dat=curr(:,ii*4-3:ii*4);
        ma=max(dat,[],2);
        mas=[mas;ma];
        
    end
    tmp=find(mas==0);
        mas(tmp)=0.0001;
    coll=[coll;mas'];
    %             curr2=flipud(reshape(curr,4,6));
    loglog([0.075 0.15 0.3 0.6 1.5 3],mas,'-','color',cols(c,:))
    hold on
    axis([0.01 10 0.01 1])
    
    %     title([dates{gofor(c)},' / ',int2str(gofor(c))])
end
%  loglog([0.075 0.15 0.3 0.6 1.5 3],coll)
%     axis([0.01 10 0.1 1])
legend(undates)
%% automatically find parasol, midget, blueyellow, others => futher compare
% with temporal properties

close all
clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
load(fullfile(superpath,'h_info'))

parasol=[0.95 1 1 1 0.65 0.3];
bluey=[1 0.9 0.9 0.75 0.2 0.01];
midget=[0.6 0.7 0.85 1 0.7 0.1];
plotit=0;

cellType=zeros(length(togo),1);
for c=1:length(togo)
    curr=normPeaks( c,:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    curr2=flipud(reshape(curr,4,6));
    spat=max(curr2,[],1);
    spat=fliplr(spat);
    spat4=fliplr(curr2(3,:));
    spat8=fliplr(curr2(4,:));
    
    
    
    pDiff=parasol-spat;
    mDiff=midget-spat;
    bDiff=bluey-spat;
    
    Sum(1)=sum(abs(pDiff(1:4)));
    Sum(2)=sum(abs(mDiff(1:4)));
    Sum(3)=sum(abs(bDiff(1:4)));
    Max(1)=max(abs(pDiff(1:4)));
    Max(2)=max(abs(mDiff(1:4)));
    Max(3)=max(abs(bDiff(1:4)));
    
    if spat(5)<0.5
        possib=[0 2 3];
    elseif spat(6)>0.5
        possib=[0 4];
    else
        possib=[0 1 2];
    end
    
    possib2=find(Max<0.2);
    if isempty(possib2)
        possib2=0;
    else
        possib2=[possib2 4];
    end
    
    possib3=find(Sum<0.5);
    if isempty(possib3)
        possib3=0;
    else
        possib3=[possib3 4];
    end
    
    theOne=intersect(possib2,possib3);
    theOne=intersect(possib,theOne);
    
    if isempty(theOne)
        theOne=0;
    elseif length(theOne)>1
        [~,id]=min(Sum);
        if ~isempty(find(theOne==id))
            theOne=id;
        else
            Sum(id)=1;
            [~,theOne]=min(Sum);
        end
    end
    
    cellType(c)=theOne;
    
    if plotit==1
        figure
        
        loglog([0.075 0.15 0.3 0.6 1.5 3],parasol,'-k')
        hold on
        loglog([0.075 0.15 0.3 0.6 1.5 3],midget,'-m')
        loglog([0.075 0.15 0.3 0.6 1.5 3],bluey,'-b')
        loglog([0.075 0.15 0.3 0.6 1.5 3],spat,'-g','linewidth',4)
        loglog([0.075 0.15 0.3 0.6 1.5 3],spat8/max(spat8),'-c','linewidth',2)
        loglog([0.075 0.15 0.3 0.6 1.5 3],spat4/max(spat4),'-g','linewidth',2)
        axis([0.01 10 0.1 1])
        saveas(gcf,fullfile('R:\Users\Katja\human_retina\spatial_curves',[int2str(c),'.bmp']))
        close
    end
    
    
end

%% adjust cell types by considering temporal aspects
cellType(:,2)=cellType(:,1);
for c=1:length(cellType)
    if cellType(c,1)==1 || cellType(c,1)==3
        curr=normPeaks( c,:,2);
        
        %log
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
        curr2=flipud(reshape(curr,4,6));
        temp=max(curr2,[],2);
        if temp(1)<0.5 && cellType(c,1)==3
            cellType(c,2)=0;
        end
        if temp(1)<0.2 && cellType(c,1)==1
            cellType(c,2)=0;
        end
        if temp(2)<0.5
            cellType(c,2)=0;
        end
        if temp(2)<0.5
            cellType(c,2)=0;
        end
        if sum(temp(3:4)>0.7)<2
            cellType(c,2)=0;
        end
        
        
        
    end
end

%% plot temporal, spatial and heatmap
close all
list={'parasol (midget)','blue-yellow','type 4','others','possible midget'};

cT=[1 3 4 0 1];
for c=1:length(cT)
    if c<4
    cand=find(cellType(:,2)==cT(c));
    elseif c==4
        cand=find(cellType(:,1)==0);
    elseif c==5
        tmp=find(cellType(:,1)==1);
        cand=find(cellType(tmp,2)==0);
        cand=tmp(cand);
    end
    
    coll=[]; coll2=[]; coll3=[];
    for i=1:length(cand)
        curr=normPeaks(cand(i),:,2);
        
        %log
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
        coll2=[coll2;curr];
        
        curr2=flipud(reshape(curr,4,6));
        spat=max(curr2,[],1);
        temp=max(curr2,[],2);
        coll=[coll;spat];
        coll3=[coll3;temp'];
    end
    figure
    subplot(2,2,1)
    loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(coll))
    hold on
    loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(median(coll)),'linewidth',4)
    axis([0.01 10 0.1 1])
    title([list{c}, ' (spatial) ',int2str(length(cand)),' cells'])
    
    subplot(2,2,2)
    loglog([1 2 4 8],fliplr(coll3))
    title('max')
    hold on
    loglog([1 2 4 8],fliplr(median(coll3)),'linewidth',4)
     axis([1 10 0.1 1])
     title('temporal')
    
    subplot(2,2,3)
    curr2=reshape(mean(coll2),4,6);
    %
    
    imagesc(flipud(fliplr(curr2)))
    set(gcf,'position',[1 41 1600 784])
    saveas(gcf,fullfile('R:\Users\Katja\human_retina\Figs',['cellType ',list{c},'.emf']))
end
%% plot all 0 cells
close all
    cand=find(cellType(:,1)==0);
    
figure
cnt=0;
    for i=1:length(cand)
        cnt=cnt+1;
        curr=normPeaks(cand(i),:,2);
        
        %log
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
%         coll2=[coll2;curr];
    
        if cnt>16
            cnt=1;
            figure
        end
    subplot(4,4,cnt)
    curr2=reshape(curr,4,6);
    imagesc(flipud(fliplr(curr2)))
    title(int2str(cand(i)))
    
    end


%% CORR: find most parasol-like cell per experiment
close all
clear
% superpath='R:\Users\Katja\human_retina';
superpath = 'C:\Users\KatjaR\OneDrive - imec\HumRet'
% load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
load(fullfile(superpath,['h','_normPeaks']))
load(fullfile(superpath,'h_info'))
dates=info(togo,1);
undates=unique(dates);


spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];


plotit=1;

candidates=[];
for i=1:length(undates)
    currD=undates{i};
    
    
    m=0; found=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(dates{m},currD))
            found=1;
        end
    end
    start=m;
    found=0;
    while found==0 && m<length(dates)
        m=m+1;
        if isempty(regexp(dates{m},currD))
            found=1;
        end
    end
    list=start:m-1;
    
    cand=[];
    for j=1:length(list)
        curr=normPeaks(list(j),:,2);
        
        %log
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
        
        
        curr2=flipud(reshape(curr,4,6));
        
        suSpat=sum(curr2,1);
        suTemp=sum(curr2,2);
        [~,maSpat]=max(suSpat);%1=100um
        [~,maTemp]=max(suTemp);%1=8Hz
        spat=curr2(maTemp,:);
        temp=curr2(:,maSpat)';

        
        if sum(spat(2:end)>0.5)==5 && sum(spat(4:end)>0.85)==3
            if spat(1)<0.4
                cand=[cand;list(j)];
            end
        end
    end
    
    if plotit==1
        figure
        subs=ceil(length(cand)/2);
        for c=1:length(cand)
            subplot(2,subs,c)
            curr=normPeaks(cand(c),:,2);
            
            %log
            tmp=find(curr==0);
            curr=log10(curr*1000);
            curr(tmp)=0;
            curr=mean(curr,1);
            curr=curr/max(curr);
              curr2=flipud(reshape(curr,4,6));
              
            suSpat=sum(curr2,1);
        suTemp=sum(curr2,2);
        [~,maSpat]=max(suSpat);%1=100um
        [~,maTemp]=max(suTemp);%1=8Hz
        spat=curr2(maTemp,:);
        temp=curr2(:,maSpat)';
            tmp=find(spat==0);
            spat(tmp)=0.0001;
            %             curr2=flipud(reshape(curr,4,6));
            curr=flipud(reshape(curr,4,6));
            %             imagesc(curr)
            loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(spat))
            axis([0.01 10 0.1 1])
            title([currD,' / ',int2str(cand(c))])
        end
        
        if isempty(cand)
            subs=ceil(length(list)/2);
            for c=1:length(list)
                subplot(2,subs,c)
                curr=normPeaks(list(c),:,2);
                
                %log
                tmp=find(curr==0);
                curr=log10(curr*1000);
                curr(tmp)=0;
                curr=mean(curr,1);
                curr=curr/max(curr);
                
                     curr2=flipud(reshape(curr,4,6));
              
            suSpat=sum(curr2,1);
        suTemp=sum(curr2,2);
        [~,maSpat]=max(suSpat);%1=100um
        [~,maTemp]=max(suTemp);%1=8Hz
        spat=curr2(maTemp,:);
        temp=curr2(:,maSpat)';
            tmp=find(spat==0);
            spat(tmp)=0.0001;
                %             curr2=flipud(reshape(curr,4,6));
                curr=flipud(reshape(curr,4,6));
                %             imagesc(curr)
                loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(spat))
                axis([0.01 10 0.1 1])
                title([currD,' / ',int2str(list(c)),' (all)'])
            end
        end
    else
        figure
        subs=ceil(length(gofor)/2);
        for c=1:length(cand)
            subplot(2,subs,c)
            curr=normPeaks(cand(c),:,2);
            
            %log
            tmp=find(curr==0);
            curr=log10(curr*1000);
            curr(tmp)=0;
            curr=mean(curr,1);
            curr=curr/max(curr);
            
           curr2=flipud(reshape(curr,4,6));
              
            suSpat=sum(curr2,1);
        suTemp=sum(curr2,2);
        [~,maSpat]=max(suSpat);%1=100um
        [~,maTemp]=max(suTemp);%1=8Hz
        spat=curr2(maTemp,:);
        temp=curr2(:,maSpat)';
            tmp=find(spat==0);
            spat(tmp)=0.0001;
            %             curr2=flipud(reshape(curr,4,6));
            curr=flipud(reshape(curr,4,6));
            %             imagesc(curr)
            loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(spat))
            axis([0.01 10 0.1 1])
            title([currD,' / ',int2str(cand(c))])
        end
        
    end
    
    
    
end
%% ... continued
close all
clear
superpath='R:\Users\Katja\human_retina';
load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
load(fullfile(superpath,'h_info'))
dates=info(togo,1);
undates=unique(dates);

gofor=[13 34 48 67 84 100 118 131 161 189 220 233 249 256 270];

spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];
cols=colormap;
cols=cols(1:3:end,:);

figure
coll=[];
for c=1:length(gofor)
    curr=normPeaks(gofor(c),:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    
     curr2=flipud(reshape(curr,4,6));
              
            suSpat=sum(curr2,1);
        suTemp=sum(curr2,2);
        [~,maSpat]=max(suSpat);%1=100um
        [~,maTemp]=max(suTemp);%1=8Hz
        spat=curr2(maTemp,:);
        temp=curr2(:,maSpat)';
            tmp=find(spat==0);
            spat(tmp)=0.01;
    %             curr2=flipud(reshape(curr,4,6));
    loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(spat),'-r')
    hold on
    axis([0.01 10 0.01 1])
    
    %     title([dates{gofor(c)},' / ',int2str(gofor(c))])
end
%  loglog([0.075 0.15 0.3 0.6 1.5 3],coll)
%     axis([0.01 10 0.1 1])
% legend(undates)    
 loglog([0.075 3],[0.1 0.1],'-k')
 saveas(gcf,fullfile('R:\Users\Katja\human_retina\Figs','parasol_like_corr.emf'))

%% CORR: automatically find parasol, midget, blueyellow, others => futher compare
% with temporal properties

close all
clear
% superpath='R:\Users\Katja\human_retina';
superpath='C:\Users\KatjaR\OneDrive - imec\HumRet';
% load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
load(fullfile(superpath,['h','_normPeaks']))
load(fullfile(superpath,'h_info'))

parasol=[0.95 1 1 1 0.65 0.3];
bluey=[1 0.9 0.9 0.75 0.2 0.01];
midget=[0.6 0.7 0.85 1 0.7 0.1];
plotit=0;

cellType=zeros(length(togo),1);
for c=1:length(togo)
    curr=normPeaks( c,:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    curr2=flipud(reshape(curr,4,6));
    spat=max(curr2,[],1);
    spat=fliplr(spat);
    spat4=fliplr(curr2(3,:));
    spat8=fliplr(curr2(4,:));
    
    
    
    pDiff=parasol-spat;
    mDiff=midget-spat;
    bDiff=bluey-spat;
    
    Sum(1)=sum(abs(pDiff(1:4)));
    Sum(2)=sum(abs(mDiff(1:4)));
    Sum(3)=sum(abs(bDiff(1:4)));
    Max(1)=max(abs(pDiff(1:4)));
    Max(2)=max(abs(mDiff(1:4)));
    Max(3)=max(abs(bDiff(1:4)));
    
    if spat(5)<0.5
        possib=[0 2 3];
    elseif spat(6)>0.5
        possib=[0 4];
    else
        possib=[0 1 2];
    end
    
    
    possib2=find(Max<0.2);
    if isempty(possib2)
        possib2=0;
    else
        possib2=[possib2 4];
    end
    
    possib3=find(Sum<0.5);
    if isempty(possib3)
        possib3=0;
    else
        possib3=[possib3 4];
    end
    
    theOne=intersect(possib2,possib3);
    theOne=intersect(possib,theOne);
    
    if isempty(theOne)
        theOne=0;
    elseif length(theOne)>1
        [~,id]=min(Sum);
        if ~isempty(find(theOne==id))
            theOne=id;
        else
            Sum(id)=1;
            [~,theOne]=min(Sum);
        end
    end
    
    cellType(c)=theOne;
    
    if plotit==1
        figure
        
        loglog([0.075 0.15 0.3 0.6 1.5 3],parasol,'-k')
        hold on
        loglog([0.075 0.15 0.3 0.6 1.5 3],midget,'-m')
        loglog([0.075 0.15 0.3 0.6 1.5 3],bluey,'-b')
        loglog([0.075 0.15 0.3 0.6 1.5 3],spat,'-g','linewidth',4)
        loglog([0.075 0.15 0.3 0.6 1.5 3],spat8/max(spat8),'-c','linewidth',2)
        loglog([0.075 0.15 0.3 0.6 1.5 3],spat4/max(spat4),'-g','linewidth',2)
        axis([0.01 10 0.1 1])
%         saveas(gcf,fullfile('R:\Users\Katja\human_retina\spatial_curves',[int2str(c),'.bmp']))
%         close
    end
    
    
end

%% adjust cell types by considering temporal aspects
cellType(:,2)=cellType(:,1);

for c=1:length(cellType)
    if cellType(c,1)==1 || cellType(c,1)==3
        curr=normPeaks( c,:,2);
        
        %log
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
        curr2=flipud(reshape(curr,4,6));
        temp=max(curr2,[],2);
        if temp(1)<0.5 && cellType(c,1)==3
            cellType(c,2)=0;
        end
        if temp(1)<0.2 && cellType(c,1)==1
            cellType(c,2)=0;
        end
        if temp(2)<0.5
            cellType(c,2)=0;
        end
        if temp(2)<0.5
            cellType(c,2)=0;
        end
        if sum(temp(3:4)>0.7)<2
            cellType(c,2)=0;
        end
        
        
        
    end
end

%% plot temporal, spatial and heatmap
close all
list={'parasol (midget)','blue-yellow','type 4','others','possible midget'};

cT=[1 3 4 0 1];
for c=1:length(cT)
    if c<4
    cand=find(cellType(:,2)==cT(c));
    elseif c==4
        cand=find(cellType(:,1)==0);
    elseif c==5
        tmp=find(cellType(:,1)==1);
        cand=find(cellType(tmp,2)==0);
        cand=tmp(cand);
    end
    
    coll=[]; coll2=[]; coll3=[];
    for i=1:length(cand)
        curr=normPeaks(cand(i),:,2);
        
        %log
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
        coll2=[coll2;curr];
        
        curr2=flipud(reshape(curr,4,6));
        spat=max(curr2,[],1);
        temp=max(curr2,[],2);
        coll=[coll;spat];
        coll3=[coll3;temp'];
    end
    figure
    med = median(coll);
    tmp = find(med==0); med(tmp)=0.001;
    for cc = 1:size(coll,1)
        tmp = find(coll(cc,:)>0);
        if tmp(1)>1
           coll(cc,tmp(1)-1)=0.001; 
        end
        if tmp(end)<6
            coll(cc,tmp(end)+1:end)=0.001;
        end
        
         tmp = find(coll3(cc,:)>0);
        if tmp(1)>1
           coll3(cc,tmp(1)-1)=0.001; 
        end
        if tmp(end)<4
            coll3(cc,tmp(end)+1:end)=0.001;
        end
%          tmp = find(coll(cc,:)==0);
%          if ~isempty(tmp) 
%              if tmp(end)==size(coll,2) || isempty(find(coll(cc,tmp(end)+1:end)>0))
%          coll(cc,tmp(end))=0.001;
%              end
%          end
    end
    subplot(2,2,1)
    loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(coll))
    hold on
    loglog([0.075 0.15 0.3 0.6 1.5 3],fliplr(med),'linewidth',4)
    axis([0.01 10 0.1 1])
    title([list{c}, ' (spatial) ',int2str(length(cand)),' cells'])
    
    subplot(2,2,2)
    med = median(coll3);
    tmp = find(med==0); med(tmp)=0.001;
    loglog([1 2 4 8],fliplr(coll3))
    title('max')
    hold on
    loglog([1 2 4 8],fliplr(med),'linewidth',4)
     axis([1 10 0.1 1])
     title('temporal')
    
    subplot(2,2,3)
    curr2=reshape(mean(coll2),4,6);
    if c==2
        curr2(1:end,1)=0;
    end
    %
    
    imagesc(flipud(fliplr(curr2)))
    set(gcf,'position',[1 41 1600 784])
    colormap gray
    saveas(gcf,fullfile('C:\Users\KatjaR\OneDrive - imec\HumRet\Figs',['final_cellType ',list{c},'.emf']))
end
%% plot all 0 cells
close all
    cand=find(cellType(:,1)==0);
    
figure
cnt=0;
    for i=1:length(cand)
        cnt=cnt+1;
        curr=normPeaks(cand(i),:,2);
        
        %log
        tmp=find(curr==0);
        curr=log10(curr*1000);
        curr(tmp)=0;
        curr=mean(curr,1);
        curr=curr/max(curr);
%         coll2=[coll2;curr];
    
        if cnt>16
            cnt=1;
            figure
        end
    subplot(4,4,cnt)
    curr2=reshape(curr,4,6);
    imagesc(flipud(fliplr(curr2)))
    title(int2str(cand(i)))
    
    end


    %% test
    allSpat =[];
    dat = zeros(length(togo),24);
for c=1:length(togo)
    curr=normPeaks( c,:,2);
    
    %log
    tmp=find(curr==0);
    curr=log10(curr*1000);
    curr(tmp)=0;
    curr=mean(curr,1);
    curr=curr/max(curr);
    dat(c,:) = curr;
    curr2=flipud(reshape(curr,4,6));
    spat=max(curr2,[],1);
    
    allSpat = [allSpat;spat];
end
heatMean=mean(data)/max(mean(data));
heatMeanR=reshape(heatMean,4,6);
curr2=heatMeanR;
suSpat=sum(curr2,1);
suTemp=sum(curr2,2);
[~,maSpat]=max(suSpat);%1=100um
[~,maTemp]=max(suTemp);%1=8Hz

spat=curr2(maTemp,:);

figure
semilogx([0.07 0.13 0.27 0.53 1.33 2.66],fliplr(median(allSpat)))
hold on
figure
semilogx([0.07 0.13 0.27 0.53 1.33 2.66],fliplr(spat))