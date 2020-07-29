figure('position',[987 1311 560 99.3333])
rectangle('position',[-1 1 4 1885],'facecolor','w','edgecolor','none')
hold on
t = -3*pi:0.01:3*pi; t= t';
imagesc(sin(t))
colormap gray
axis([-1 3 -inf inf])
axis off
box off
exportgraphics(gcf,fullfile('C:\Users\katja\OneDrive - imec\HumRet\PlosOne\stim','Stim5.png'),'resolution',600)


figure('position',[987 1311 560 99.3333])
rectangle('position',[-1 1 4 1885],'facecolor','w','edgecolor','none')
hold on
t = -3*pi:0.01:3*pi; t= t';
imagesc(sin(t))
colormap gray
hold on
text(-0.1,900,'best','fontname','Arial','fontsize',20)
axis([-1 3 -inf inf])
axis off
box off
exportgraphics(gcf,fullfile('C:\Users\katja\OneDrive - imec\HumRet\PlosOne\stim','Stim6.png'),'resolution',600)