clear
close all
% path1='/media/NERFFS01/Users/Katja/human_retina';
path1='R:\Users\Katja\human_retina';
load(fullfile(path1,'h_responding'))

%%

for r=1:size(DGspikes,3)
                    spikes=DGspikes{i,dg,r};
                    spikes=round(spikes);
                    rate=zeros(1,14500);
                    rate(spikes)=1;
                    rate(spikes(diff(spikes)==0))=2;
                    allrate=[allrate;rate(1:14500)];
                end
                meanrate1=mean(allrate);
                meanrate=meanrate1(2200:14200);
            end