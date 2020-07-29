function [mat_fit,xhat]=fitsftf(mat_raw,sf,tf)
% fit of 2d gaussian in sftf space

tf=tf(end:-1:1);

f0=mat_raw;
f0(f0<=0)=0.0001; % remove negative values

[sfmat,tfmat] = meshgrid(sf,tf);

f = @(par,sfmat) mymodel_sftf(par,sfmat,tfmat);
x0 =  [1,1.5,1,3,1,-.5];   % Amp,mu_sf,width_sf,mu_tf,width_tf,Q(tan)
% x0=[1,1500,1000,3,1,-.5];
lb =  [0,min(sf),.01,min(tf),.1,-1]; % lower boundary
ub =  [1000,max(sf),1,max(tf),20,1]; % upper boundary
opt = optimset('display','on');
xhat = lsqcurvefit(f,x0,sfmat,f0,lb,ub,opt);
mat_fit = mymodel_sftf(xhat,sfmat,tfmat);


