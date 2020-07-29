clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
load(fullfile(path1,'h_responding'))

%% spont background

load(fullfile(path1,'h_spont'))

load(fullfile(path1,['h_info']))
rDates=info(resp,1);
rUnits=info(resp,2);

cd('R:\Users\Katja\human_retina\_scripts')

meanSpont=[];
for g=1:length(spontRates)
    data=spontRates{g};
    m=mean(data);
    meanSpont=[meanSpont;m];
end

non=setdiff(1:length(meanSpont),resp);
tmp=find(isnan(meanSpont(non)));
tmp2=1:length(non);
tmp2=setdiff(tmp2,tmp);
non=non(tmp2);

figure
plot([0.8 1.2],[mean(meanSpont(non)) mean(meanSpont(non))],'-r')
hold on
plot([1.8 2.2],[mean(meanSpont(resp)) mean(meanSpont(resp))],'-g')
plot([0.8 1.2],[median(meanSpont(non)) median(meanSpont(non))],'--r')
plot([1.8 2.2],[median(meanSpont(resp)) median(meanSpont(resp))],'--g')
plot([1 1],[min(meanSpont(non)) max(meanSpont(non))],'-r')
plot([2 2],[min(meanSpont(resp)) max(meanSpont(resp))],'-g')
plot([1 1],[mean(meanSpont(non))-std(meanSpont(non)) mean(meanSpont(non))+std(meanSpont(non))],'-r','linewidth',2)
plot([2 2],[mean(meanSpont(resp))-std(meanSpont(resp)) mean(meanSpont(resp))+std(meanSpont(resp))],'-g','linewidth',2)
ci1=-1.96/sqrt(length(non));
ci2=-1.96/sqrt(length(resp));
% plot([1 1],[mean(meanSpont(non))-ci1 mean(meanSpont(non))+ci1],'-r','linewidth',4)
% plot([2 2],[mean(meanSpont(resp))-ci2 mean(meanSpont(resp))+ci2],'-g','linewidth',4)
axis([0 3 -inf inf])

name=['backgroundFiring_nonResponder_responder.emf'];
% saveas(gcf,fullfile(path1,'Figs',name))



[hnon,x]=hist([meanSpont(non); 60],15);
hnon=hnon/length(non);
hresp=hist([meanSpont(resp); 60],15);
hresp=hresp/length(resp);
figure
bar(x,[hnon' hresp'])
legend('non','resp')
