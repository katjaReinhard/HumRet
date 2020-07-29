clear
close all


Dtype='h'; % h = human, p = pig, w = wildtype mouse
% superpath='S:\data\Katja\allspecies';
% superpath='/gpfs01/muench/data/Katja/allspecies';
% superpath='C:\Users\Katja\Desktop\all_species_temp';
% superpath='S:\data\Katja\all_species';
superpath='R:\Users\Katja\human_retina';
%% data



switch Dtype
    
    
    case 'h'
        % Human data
        dates={'20120210a','20120224','20120320c','20120608b','20130124a','20130211a',...,
            '20130211c','20130211d','20130510a','20130510d','20130510e','20130718b','20130718c','20131210a','20131217f'};
        mainpath='S:\data\Katja\Human';
        starts=[255 245 8 91 1 97 192 97 1 1 1 173 1 1 88]; %index as in unit files
        stops=[374 365 277 280 191 190 449 285 246 164 246 429 172 86 260];
        NDs=[4 4 4 4 4 3 4 4 4 4 4 3 4 4 3];
        dnum=3;
    case 'p'
        % Pig data
        dates={'20121127','20121211a','20130613a','20130613b','20130613c','20140723d','20120213','20120319a','20120319b'};
        mainpath='S:\data\Katja\Pig';
        starts=[1 1 1 1 1 156 383 1 361];
        stops=[190 190 164 164 164 285 502 270 630];
        NDs=[4 4 4 4 4 4 3 3 3];
        dnum=2;
        
    case 'm'
        % Mouse data OPA
        dates={'20120126','20120127','20120221','20120222'};
        mainpath='S:\data\Katja\OPA';
        starts=[363 363 363 363];
        stops=[543 543 543 543];
        NDs=[4 4 4 4];
    case 'w'
        % Mouse data C3H
        dates={'20150114b','20150128b','20150203b','20150324b','20150424a','20150507a','20150507b'};
        %         mainpath='S:\data\Katja\WT';
        mainpath='C:\Users\Katja\Desktop\all_species_temp';
        starts=[412 412 412 412 412 412 412];
        stops=[577 577 577 577 577 577 577];
        NDs=[4 4 4 4 4 4 4];
        dnum=1;
end

%% !!clustering kmeans
clear
% path1='F:\all_species_temp\newAttempt\DG_Gaussian';
% type='GaussianNew2';
% path1='F:\all_species_temp\newAttempt\DG_PCanalysis';
type='amp_F1';
% path1='F:\all_species_temp\newAttempt6\DG_PCanalysis';
path1='R:\Users\Katja\human_retina\newAttemptOnlyHum\DG_PCanalysis';
% type='means_F1';

outlier=1;

load(fullfile(path1,type))
savename='clust_Kmeans_';
close all

orig=amps;
gau=0;
% orig=Gauss;
% gau=1;
% orig=amps;
% gau=0;

if gau==1
    
    orig2clust=[]; tID=zeros(sum(datalength),1);
    last=1;
    for i=1:3
        curr=reshape(Gauss(i,:,:),300,5); %centerX, centerY, x-axis, y-axis, angle
        curr=curr(1:datalength(i),:);
        orig2clust=[orig2clust;curr];
        tID(last:last+datalength(i)-1)=repmat(i,datalength(i),1);
        last=last+datalength(i);
    end
    tmp=find(orig2clust(:,5)<0);
    tmp2=abs(-deg2rad(180)-orig2clust(tmp,5));
    orig2clust(tmp,5)=tmp2;
    keepOrig=orig2clust;
else
    orig2clust=orig;
end
if outlier==1
    for p=1:size(orig2clust,2)-1
        curr=orig2clust(:,p);
        norm=(curr-min(curr))/(max(curr)-min(curr));
        orig2clust(:,p)=norm;
    end
    curr=orig2clust(:,end);
    [tmp tmp2]=sort(curr);
    norm=(curr-min(curr))/(tmp(end-1)+0.2-min(curr));
    norm(tmp2(end))=1;
    orig2clust(:,end)=norm;
else
    for p=1:size(orig2clust,2)
        curr=orig2clust(:,p);
        norm=(curr-min(curr))/(max(curr)-min(curr));
        orig2clust(:,p)=norm;
    end
end
opts=statset('maxiter',1000);
ID=[];
for k=2:15
    [idx c]=kmeans(orig2clust,k,'emptyaction','drop','start','cluster','replicates',10000,'options',opts);
    k
    ID=[ID idx];
end

if gau==1
    ang=rad2deg(keepOrig(:,5));
    keepOrig(:,5)=ang;
else
    keepOrig=orig2clust;
end
save(fullfile(path1,[savename,type]),'orig2clust','keepOrig','ID','tID')


