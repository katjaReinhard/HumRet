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
charact=zeros(length(normPeaks),4);

spatial=[100 200 500 1000 2000 4000];
temporal=[1 2 4 8];
for i=1:30
    curr=normPeaks(i,:,2);
    
    %     curr2=flipud(reshape(curr,4,6));
    curr2=reshape(curr,4,6);
    
    figure
    imagesc(curr2)
    
    %spatial
    sums=sum(curr2,2);
    [ma id]=max(sums); %this is the temporal freq for which we check
    cumu=cumsum(curr2(id,:));
    cumu50=find(cumu>=cumu(end)/2);
    cumu50=cumu50(1);
    charact(i,1)=cumu50;
    
    sums(id)=0;
    [ma id]=max(sums); %this is the temporal freq for which we check (2nd highest)
    cumu=cumsum(curr2(id,:));
    cumu50=find(cumu>=cumu(end)/2);
    cumu50=cumu50(1);
    charact(i,2)=cumu50;
    
    %temporal
    sums=sum(curr2,1);
    [ma id]=max(sums); %this is the spatial freq for which we check
    cumu=cumsum(curr2(:,id));
    cumu50=find(cumu>=cumu(end)/2);
    cumu50=cumu50(1);
    charact(i,3)=cumu50;
    
    
    sums(id)=0;
    [ma id]=max(sums); %this is the temporal freq for which we check (2nd highest)
    cumu=cumsum(curr2(:,id));
    cumu50=find(cumu>=cumu(end)/2);
    cumu50=cumu50(1);
    charact(i,4)=cumu50;
    
    title(['spatial: ',int2str(spatial(charact(i,1))),' / ',int2str(spatial(charact(i,2))) '; temporal: ',int2str(temporal(charact(i,3))),' / ',int2str(temporal(charact(i,4)))])
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

% x0 =  [1,.2,.1,4,2,-.5]; 
pars=[1,1.5,0.7,3,1.5,-.5];
A = pars(1); % amplitude
mu_sf = pars(2); % center of SF
sigma_sf = pars(3);  % sigma of SF tuning
mu_tf = pars(4);  % center of TF
sigma_tf = pars(5);   % sigma of TF tuning 
Q = pars(6); % -1 to 1, the tan of x/y

% computer terms
tfp = 2.^((Q+1)*( log2(sfmat) - log2(mu_sf) ) + log2(mu_tf));
t1 = exp( -(log2(sfmat)-log2(mu_sf)).^2 / log2(sigma_sf).^2 );
t2 = exp( -(log2(tfmat)-log2(tfp)).^2 / log2(sigma_tf).^2 );

R = A*t1.*t2;





[mat_fit,xhat]=fitsftf(curr2,sf,tf);

formula=xhat(1)*exp(-((tfmat-xhat(4)).^2/(2*xhat(5)^2)+(sfmat-xhat(2)).^2/(2*xhat(3)^2)));

end
