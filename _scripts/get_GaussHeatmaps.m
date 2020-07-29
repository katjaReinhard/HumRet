function [allCenterX,allCenterY,allContX,allContY,dimX,dimY] = get_GaussHeatmaps(normDG)


st = 0; fact=2;

allCenterX = []; allCenterY = []; allContX = cell(1,1); allContY = cell(1,1);
dimX = []; dimY = [];
for nn = 1:length(normDG)
    %%
    curr=normDG(nn,:,2);
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
        tmp=find(D==0);
        D(tmp)=0.1;
    end
    
    
    aa = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
    bb = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
    %                 bb = cos(theta)^2/2/(D(1)/1.5)^2 + sin(theta)^2/2/(D(2)/1.5)^2;
    %                 aa = -sin(2*theta)/4/(D(1)/1.5)^2 + sin(2*theta)/4/(D(2)/1.5)^2 ;
    cc = sin(theta)^2/2/(D(1)/1.5)^2 + cos(theta)^2/2/(D(2)/1.5)^2;
    allClust = exp( - (aa*(Y-centerY).^2 + 2*bb*(Y-centerY).*(X-centerX2) + cc*(X-centerX2).^2)) ;
    
    allCenterX = [allCenterX;centerX]; allCenterY = [allCenterY;centerY];
     dimX = [dimX;D(1)]; dimY = [dimY;D(1)];
    figure
    c = contour(X,Y,allClust(:,:),[0.5 0.5]);
    allContX{nn} = c(1,2:end);
    allContY{nn} = c(2,2:end);
    close
    
end
% contour(X,Y,allClust(:,:),1,'color','r')
% hold on
% plot(centerX,centerY,'or','markerfacecolor','r','markersize',2)