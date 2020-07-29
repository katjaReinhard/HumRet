%% 0) phase calculations of simple sine
fs2=2000;
T=4;
% phase=1;
phase=pi/2;
frq1=1;
t=1/fs2:1/fs2:T;
resp1=sin(2*pi*frq1*t+phase);
resp2=sin(2*pi*frq1*2*t+phase);
% resp2=0;
x=resp1+resp2;
figure;subplot(2,2,1);plot(t,resp1,'-k'); hold on;plot(t,resp2,'-c');plot(t,x,'linewidth',3)
title(['f1: ',int2str(frq1),' / f2: ',int2str(frq1*2)])

n2=5*fs2;
y2=fft(x,n2)/numel(x);
f2=fs2/2*linspace(0,1,n2/2+1);
y2=y2(1:n2/2+1);
c2=2*abs(y2);
subplot(2,2,3);plot(f2,c2)
axis([0 20 -inf inf])

y02=fftshift(y2);
 f02 = (-n2/4:n2/4)*(fs2/n2);
 % f02=linspace(0,2,n2/2+1)*(fs2/2);
% f02=f02-max(f02)/2;
% f02=round(f02*1000)/1000;
phases2 = angle(y02);
subplot(2,2,4);plot(f02,phases2)
axis([0 20 -inf inf])
title(['phase1: ',num2str(phase),' / phase2: ',num2str(phase*2)])
%%
tmp1=find(f02==frq1);
tmp2=find(f02==frq1*2);
phs=[phases2(tmp1) phases2(tmp2)];
ratios=(phs(2)-2*pi*(-2:2))/phs(1)


