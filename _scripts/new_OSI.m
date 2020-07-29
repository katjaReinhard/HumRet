clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';

load(fullfile(path1,['h_info']))
rDates=info(:,1);
rUnits=info(:,2); 

test=[313 314 316 322 324 337 341 346 352 362 376 388 417 424 431 435 447 471 481 491 500 501 505 520 530 587 628 650 664 669 709 720 741 743 752 808 811 826 1042];
for g=1:length(test) 

currUnit=rUnits{test(g)};
    currDate=rDates{test(g)};
    startF=info{test(g),4};
    stopF=info{test(g),5};

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
     heka=unit{1,2}(:,1);
    
black=[];
    for f=startF:stopF
        if ~isempty(regexp(heka{f},'DS_alldirections_black'))
            black=[black;f];
        end
    end
 
    
    first=black(1);
    prot=read_header_field_heka(fullfile(dpath,'HEKA'),heka{first},'Stimulus Protocol');
    cons=prot(1:end,2);
    changes=find(cons==-1);
    startpoints=round(prot(changes(1:8),1));
    endpoints=round(prot(changes(2:9)-1,1));
    cols=[0.1 0.6 0.9 0.6 0.9 0.6 0.9 0.6 0.9 0.1];
    changes=[0; startpoints; endpoints(end); round(prot(end,1)+1000)];
    tt=round(prot(end,1)+1000);
    
    figure
    allrateB=[];
    for f=1:length(black)
        subplot(2,2,1)
        spikes=unit{1,2}{black(f),2};
        for ss=1:length(spikes)
           plot([spikes(ss) spikes(ss)],[f-0.4 f+0.4],'-')
           hold on
        end
        rate=MEA_spikerates(spikes,60,tt);
        allrateB=[allrateB;rate];
    end
    for s=1:length(startpoints)
        plot([startpoints(s) startpoints(s)],[0 f+1],'-k')
    end
    meanrateB=mean(allrateB);
    subplot(2,2,3)
    plot(meanrateB)
    hold on
     for s=1:length(startpoints)
        plot([startpoints(s) startpoints(s)],[0 max(meanrateB)],'-k')
    end
    
   
    
%     figure
%     set(gcf,'position',[1601 1 1600 824])
%     subplot(2,2,1)
%     for c=1:length(changes)-1
%         rectangle('position',[changes(c),0,changes(c+1)-changes(c),max(meanrateB)+2],'facecolor',[cols(c) cols(c) cols(c)],'edgecolor',[cols(c) cols(c) cols(c)])
%         hold on
%     end
%     plot(meanrateB)
%     axis tight
%     
%     subplot(2,2,2)
%     for c=1:length(changes)-1
%         rectangle('position',[changes(c),0,changes(c+1)-changes(c),max(meanrateW)+2],'facecolor',[cols(c) cols(c) cols(c)],'edgecolor',[cols(c) cols(c) cols(c)])
%         hold on
%     end
%     plot(meanrateW)
%     axis tight
%     
    backB=mean(meanrateB(1:startpoints(1)));
    newrateB=meanrateB-backB;

    
    startpoints=[startpoints; endpoints(end)];
    vals=[];
    for fli=1:length(startpoints)-1
        curr=max(newrateB(startpoints(fli):startpoints(fli+1)));
        vals=[vals;curr];
    end
    tmp=find(vals<0);
    vals(tmp)=0;
    normv=vals/max(vals);
    xcoord=(normv(1)-normv(2)+normv(8)*cos(45)+normv(5)*cos(45)-normv(6)*cos(45)-normv(7)*cos(45))/(sum(normv));
    ycoord=(normv(4)-normv(3)+normv(6)*cos(45)+normv(8)*cos(45)-normv(5)*cos(45)-normv(7)*cos(45))/(sum(normv));
    DSindex=sqrt(xcoord^2+ycoord^2);
    
    indices=[];
    pref=normv(1)+normv(2);
    orth=normv(3)+normv(4);
    diff=pref-orth;
    indices=[indices; abs(diff)];
    pref=normv(5)+normv(6);
    orth=normv(7)+normv(8);
    diff=pref-orth;
    indices=[indices; abs(diff)];
  
%     %
    subplot(2,2,2)
    plot([0 normv(1)],[0 0],'-k')
    hold on
    plot([0 -normv(2)],[0 0],'-k')
    plot([0 0],[0 -normv(3)],'-k')
    plot([0 0],[0 normv(4)],'-k')
    plot([0 normv(5)*cosd(45)],[0 -normv(5)*cosd(45)],'-k')
    plot([0 -normv(6)*cosd(45)],[0 normv(6)*cosd(45)],'-k')
    plot([0 -normv(7)*cosd(45)],[0 -normv(7)*cosd(45)],'-k')
    plot([0 normv(8)*cosd(45)],[0 normv(8)*cosd(45)],'-k')
    plot([0 xcoord],[0 ycoord],'-r','linewidth',3)
    plot([normv(1) normv(5)*cosd(45)],[0 -normv(5)*cosd(45)],'b-')
     plot([normv(5)*cosd(45) 0],[-normv(5)*cosd(45) -normv(3)],'b-')
     plot([0 -normv(7)*cosd(45)],[-normv(3) -normv(7)*cosd(45)],'b-')
     plot([-normv(7)*cosd(45) -normv(2)],[-normv(7)*cosd(45) 0],'b-')
     plot([-normv(2) -normv(6)*cosd(45)],[0 normv(6)*cosd(45)],'b-')
     plot([-normv(6)*cosd(45) 0],[normv(6)*cosd(45) normv(4)],'b-')
     plot([0 normv(8)*cosd(45)],[normv(4) normv(8)*cosd(45)],'b-')
        plot([normv(8)*cosd(45) normv(1)],[normv(8)*cosd(45) 0],'b-')
        title(['DSI: ',num2str(DSindex),' || OSI: ',num2str(indices(1)),' / ',num2str(indices(2))])
    axis([-1 1 -1 1])
    rectangle('position',[-0.2,-0.2,0.4,0.4],'curvature',[1,1],'edgecolor','g')
%     rectangle('position',[-0.35,-0.35,0.7,0.7],'curvature',[1,1],'edgecolor','g')
%     rectangle('position',[-0.5,-0.5,1,1],'curvature',[1,1],'edgecolor','g')
    rectangle('position',[-1,-1,2,2],'curvature',[1,1],'edgecolor','g')
    axis([-1 1 -1 1])
    set(gcf,'position',[1601 -179 1920 1004])
    saveas(gcf,fullfile(path1,'OSI_DSI',[int2str(test(g)),'_',currDate,'_',currUnit,'.bmp']))
    close
end