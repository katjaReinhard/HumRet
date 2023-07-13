clear
close all
clc
% creates list_date_responses, which contains:
% bar, chrp, chrp_tested,dates,DG,flash,good,good2,idxPerDate,R
%%
% path1 = 'C:\Users\Katja\OneDrive - imec\HumRet';
% path2 = 'C:\Users\Katja\OneDrive - imec\HumanRGC_processedData';
path1 = 'D:\_data\all_species';
speciesnow = 'p';
load(fullfile(path1,[speciesnow,'_info']))
%% full-field flash

good=find(goodones==1);
good2=[];
for i=1:length(good)
    if info{good(i),8}~=0
        good2=[good2;good(i)];
    end
end

flash = good2;
%% bar
load(fullfile(path1,[speciesnow,'_vel6_new_Black']))

tmp = find(~cellfun(@isempty,vel6_newB));
bar = tmp;
%% gratings
tmp = find(~cellfun(@isempty,info(:,18)));
DG = tmp;
%% chirp
chrp_tested = find(~cellfun(@isempty,info(:,19)));
tmp = cell2mat(info(chrp_tested,19));
chrp = find(tmp==1);chrp = chrp_tested(chrp);

%% all responding
R = [DG;chrp;bar;flash;]; R = unique(R);

dates = unique(info(:,1));
idxPerDate = cell(length(dates),1);
for d = 1:length(dates)
    tmp = strfind(info(:,1),dates{d});
    tmp = find(~cellfun(@isempty,tmp));
    idxPerDate{d} = tmp;
end
%% get numbers
T1 = []; %grating, chirp, bar, flash,chirp tested
for d = 1:length(dates)
    id = idxPerDate{d};
    tmp1 = intersect(id,DG);
    tmp2 = intersect(id,chrp);
    tmp3 = intersect(id,bar);
    tmp4 = intersect(id,flash);
    tmp5 = intersect(id,chrp_tested);
    tmp6 = intersect(id,R);
    T1 = [T1;length(tmp1) length(tmp2) length(tmp3) length(tmp4) length(tmp5) length(tmp6)];
end

%%
save(fullfile(path1,[speciesnow,'_list_date_responses']),'bar', 'chrp', 'chrp_tested','dates','DG','flash','good','good2','idxPerDate','R')




