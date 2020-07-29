
%% CORR: automatically find parasol, midget, blueyellow, others => futher compare
% with temporal properties

close all
clear
% superpath='R:\Users\Katja\human_retina';
superpath='C:\Users\KatjaR\OneDrive - imec\HumRet';
% load(fullfile(superpath,'newAttempt','DG',['h','_normPeaks']))
load(fullfile(superpath,['w','_normPeaks']))
load(fullfile(superpath,'w_info'))

parasol=[0.95 1 1 1 0.65 0.3];
% parasol = 10.^parasol; parasol = parasol/(max(parasol));
bluey=[1 0.9 0.9 0.75 0.2 0.01];
% bluey = 10.^bluey; bluey = bluey/(max(bluey));
midget=[0.6 0.7 0.85 1 0.7 0.1];
plotit=0;

col = colormap('jet');
col = col(1:3:end,:);
close


cellType=zeros(length(togo),1);
SPAT = [];
TEMP = [];
for c=1:length(togo)
    curr=normPeaks( c,:,2);
    curr2 = reshape(curr,4,6);
    suSpat=sum(curr2,1);
    suTemp=sum(curr2,2);
    [~,maSpat]=max(suSpat);%1=100um
    [~,maTemp]=max(suTemp);%1=8Hz
    spat=fliplr(curr2(maTemp,:));
    
    temp=fliplr(curr2(:,maSpat)');
    %     spat = fliplr(max(curr2));
    %     temp = fliplr(max(curr2'));
    % spat=fliplr(suSpat/max(suSpat));
    % temp=fliplr(suTemp'/max(suTemp'));
    spat = spat/max(spat);
    temp = temp/max(temp);
    SPAT = [SPAT;spat];
    TEMP = [TEMP;temp];
    
    possib = [];
    % 4 = broad
    % 3 = blue-yellow
    % 0 = other
    % 1 = parasol
    % 2 = midget
    if spat(5)<0.2 && spat(6)<0.2
        possib=[0 3];
    elseif spat(5)>0.2 && spat(6)>0.2
        possib=[0 4];
    elseif spat(5)>0.4
        possib=[0 1];
    else
        possib=[0];
    end
    
    possib2=find(spat(1:4)<0.2);
    if ~isempty(possib2)
        possib=0;
    elseif length(possib)>1
        possib=setdiff(possib,0);
    end
    
    cellType(c)=possib;
    
    
end

%% adjust cell types by considering temporal aspects

% 4 = broad
% 3 = blue-yellow
% 0 = other
% 1 = parasol
% 2 = midget

parasol = [0.4 0.6 0.8 0.9];
bluey = [0.55 0.7 0.9 0.95];
parasol = 10.^parasol; parasol = parasol/(max(parasol));
bluey = 10.^bluey; bluey = bluey/(max(bluey));

for c=1:length(cellType)
    if cellType(c,1)==1 || cellType(c,1)==3
        %         if cellType(c,1)==1
        %             break
        %         end
        curr=normPeaks( c,:,2);
        curr2 = reshape(curr,4,6);
        suSpat=sum(curr2,1);
        suTemp=sum(curr2,2);
        [~,maSpat]=max(suSpat);%1=100um
        [~,maTemp]=max(suTemp);%1=8Hz
        spat=fliplr(curr2(maTemp,:));
        % spat = log(spat)
        temp=fliplr(curr2(:,maSpat)');
        %         spat = fliplr(max(curr2));
        %         temp = fliplr(max(curr2'));
        %         spat=fliplr(suSpat/max(suSpat));
        % temp=fliplr(suTemp'/max(suTemp'));
        spat = spat/max(spat);
    temp = temp/max(temp);
        
        if temp(1)<0.2
            cellType(c)=0;
        end
        if temp(2)<0.2
            cellType(c)=0;
        end
        
        if sum(temp(3:4)>0.3)<2
            cellType(c)=0;
        end
        
        
        
    end
end

%% plot temporal, spatial and heatmap
close all
list={'parasol (midget)','blue-yellow','type 4','others','possible midget'};
% 4 = broad
% 3 = blue-yellow
% 0 = other
% 1 = parasol
% 2 = midget


cT=[1 3 4 0];
for c=1:length(cT)
    
    cand=find(cellType==cT(c));
    
    coll2=[];
    for i = 1:length(cand)
        curr=normPeaks(cand(i),:,2);
        coll2=[coll2;curr];
    end
    
    figure
    coll = SPAT(cand,:);
    coll3 = TEMP(cand,:);
    %     med = median(coll);
    med = mean(coll);
    tmp = find(med==0); med(tmp)=0.1;
    for cc = 1:size(coll,1)
        %         tmp = find(coll(cc,:)>0);
        %         if tmp(1)>1
        %             coll(cc,tmp(1)-1)=0.001;
        %         end
        %         if tmp(end)<6
        %             coll(cc,tmp(end)+1:end)=0.001;
        %         end
        %
        %         tmp = find(coll3(cc,:)>0);
        %         if tmp(1)>1
        %             coll3(cc,tmp(1)-1)=0.001;
        %         end
        %         if tmp(end)<4
        %             coll3(cc,tmp(end)+1:end)=0.001;
        %         end
        tmp = find(coll==0);
        coll(tmp)=0.1;
        tmp = find(coll3==0);
        coll3(tmp)=0.1;
    end
    
    subplot(2,2,1)
    ct = 0;
%     for co = 1:length(coll)
%         if mod(co,10)>0 || ct>=size(col,1)
%             loglog([0.07 0.13 0.27 0.53 1.33 2.66],coll(co,:),'-k')
%         else
%             ct = ct+1;
%             loglog([0.07 0.13 0.27 0.53 1.33 2.66],coll(co,:),'-','color',col(ct,:))
%         end
%         hold on
%     end
    ct = 0;
    special = [];
    for co = 1:length(coll)
        if coll(co,6)<0.2
            loglog([0.07 0.13 0.27 0.53 1.33 2.66],coll(co,:),'-k')
        else
            ct = ct+1;
            loglog([0.07 0.13 0.27 0.53 1.33 2.66],coll(co,:),'-','color',col(ct,:))
            special = [special;co];
        end
        hold on
    end
    
    
    hold on
    loglog([0.07 0.13 0.27 0.53 1.33 2.66],med,'linewidth',4)
    axis([0.01 10 0.1 1])
    title([list{c}, ' (spatial) ',int2str(length(cand)),' cells'])
    
    subplot(2,2,2)
    %     med = median(coll3);
    med = mean(coll3);
    tmp = find(med==0); med(tmp)=0.1;
    ct = 0;
    for co = 1:length(coll3)
        if isempty(find(special==co))
            loglog([1 2 4 8],coll3(co,:),'k')
        else
            ct = ct+1;
            loglog([1 2 4 8],coll3(co,:),'-','color',col(ct,:))
        end
        hold on
    end
    title('max')
    hold on
    loglog([1 2 4 8],med,'linewidth',4)
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

figure
loglog([0.07 0.13 0.27 0.53 1.33 2.66],mean(SPAT),'linewidth',4)
axis([0.01 10 max(mean(SPAT))/10 1])
hold on
plot([0.01 10],[0.1 0.1],'-k')


figure
loglog([1 2 4 8],mean(TEMP),'linewidth',4)
axis([1 10 max(mean(TEMP))/10 1])
hold on
plot([1 10],[0.1 0.1],'-k')

%%
figure
loglog([0.07 0.13 0.27 0.53 1.33 2.66],mean(SPAT),'linewidth',4)
axis([0.01 10 max(mean(SPAT))/10 1])
hold on
plot([0.01 10],[0.1 0.1],'-k')



