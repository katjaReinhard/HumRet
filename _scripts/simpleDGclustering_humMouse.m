%%
clear
close all
path1 = 'C:\Users\KatjaR\OneDrive - imec\HumRet';
load(fullfile(path1,'h_normPeaks'))
orig2clust = normPeaks(:,:,2);
clear normPeaks
load(fullfile(path1,'w_normPeaks'))
orig2clust = [orig2clust;normPeaks(:,:,2)];

%%
Sim= squareform(pdist(orig2clust));
%%
opts=statset('maxiter',1000);
ID=[];
for k=2:20
    [idx c]=kmeans(orig2clust,k,'emptyaction','drop','start','cluster','replicates',10000,'options',opts);
    k
    ID=[ID idx];
end
%%
h = 293;
for c = 8:19
   curr = ID(:,c);
   figure
   un=unique(curr);
   for cc = 1:un(end)
       subplot(5,4,cc)
       tmp = find(curr==cc);
       dat= orig2clust(tmp,:);
       dat2=mean(dat)';
        dat3=reshape(dat2,4,6);
       imagesc(flipud(fliplr(dat3)))
       hu = length(find(tmp<=h));
       mo=length(tmp)-hu;
       title(['human: ',int2str(hu),' | mouse: ',int2str(mo)])
   end
   
end

