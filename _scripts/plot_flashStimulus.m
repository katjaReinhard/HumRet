 figure('position',[987 1311 560 99.3333])
        rectangle('position',[0 0 1 0.1],'facecolor',[0.5 0.5 0.5],'edgecolor','none')
        hold on
        rectangle('position',[1 0 1 0.1],'facecolor',[0 0 0],'edgecolor','none')
        rectangle('position',[2 0 1 0.1],'facecolor',[0.5 0.5 0.5],'edgecolor','none')
        axis([0 3 -0.2 0.3])
        axis off
        box off
        saveas(gcf,fullfile('C:\Users\katja\OneDrive - imec\HumRet\PlosOne\stim','Stim1.png'))