clear
close all
clc

%%
path1 = 'C:\Users\Katja\OneDrive - imec\HumanRetina\PlosOne';
path2 = 'C:\Users\Katja\OneDrive - imec\HumRet';
pathcode = 'C:\Users\katja\OneDrive - imec\HumRet\_scripts';
pathSave = 'C:\Users\katja\OneDrive - imec\HumanRetina\PlosOne\resubmission';

addpath(genpath(pathcode))
%% do only once (replace with actual data from Cameron)
% I = imread(fullfile(path1,'Cowan2020.png'));
%
% figure
% imagesc(I)
% coords = cell(5,2);
% for ii = 2:6
%     [x,y] = ginput;
%     coords{ii,1} = x;
%     coords{ii,2} = y;
% end
% [x,y] = ginput;
% twosec = x(2)-x(1);
% save(fullfile(path2,'Cowan2020'),'coords','twosec')
%%
prot=read_header_field_heka(fullfile('C:\Users\Katja\OneDrive - imec\HumRet\data\20130211c\HEKA'),'20130211_C1#0211_HumanPig_FFFlashBW%dur_2000_gap_2000%_FW0ND4.phys','Stimulus Protocol');
prot2=read_header_field_heka(fullfile('C:\Users\Katja\OneDrive - imec\HumRet\data\20130211c\HEKA'),'20130211_C1#0283_HumanPig_chirp%mean_128_amplitude_128_pause_500%_FW0ND4.phys','Stimulus Protocol');

%% get spike rates chrp and flash and transiency
load(fullfile(path2,'h_info'))
load(fullfile(path2,'list_date_responses'))

chrp_spike = [];
for g=1:length(chrp)
    
    currUnit=info{chrp(g),2};
    currDate=info{chrp(g),1};
    
    startF=info{chrp(g),4};
    stopF=info{chrp(g),5};
    
    %pathes%
    dpath=fullfile(path2,'data',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    %load current unit
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    %find files with chirp stimulus
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'chirp'))
            flashfiles=[flashfiles;h];
        end
    end
    tmp=find(flashfiles>startF);
    flashfiles=flashfiles(tmp);
    tmp=find(flashfiles<=stopF);
    flashfiles=flashfiles(tmp);
    
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(round(spikes),50,25000);
        allrate=[allrate;rate];
    end
    meanrate=mean(allrate,1);
    chrp_spike = [chrp_spike;meanrate];
    
end

flash_spike = []; flash_trans = [];
for g=1:length(flash)
    
    currUnit=info{flash(g),2};
    currDate=info{flash(g),1};
    
    startF=info{flash(g),4};
    stopF=info{flash(g),5};
    
    %pathes%
    dpath=fullfile(path2,'data',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    %load current unit
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    %find files with chirp stimulus
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlashBW'))
            flashfiles=[flashfiles;h];
        end
    end
    tmp=find(flashfiles>startF);
    flashfiles=flashfiles(tmp);
    tmp=find(flashfiles<=stopF);
    flashfiles=flashfiles(tmp);
    
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(round(spikes),20,12000);
        allrate=[allrate;rate];
    end
    meanrate=mean(allrate,1);
    flash_spike = [flash_spike;meanrate];
    bck = mean(meanrate(50:1500));
    stback = std( meanrate(50:1500));
    [m1,m2] = max(meanrate);
    tmp = find(meanrate(m2:end)<=bck+2*stback); tmp = tmp(1);
    flash_trans = [flash_trans;tmp];
end

flash_all = [];
for g=1:length(R)
    
    currUnit=info{R(g),2};
    currDate=info{R(g),1};
    
    startF=info{R(g),4};
    stopF=info{R(g),5};
    
    %pathes%
    dpath=fullfile(path2,'data',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    %load current unit
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    %find files with chirp stimulus
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'FFFlashBW'))
            flashfiles=[flashfiles;h];
        end
    end
    tmp=find(flashfiles>startF);
    flashfiles=flashfiles(tmp);
    tmp=find(flashfiles<=stopF);
    flashfiles=flashfiles(tmp);
    
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(round(spikes),20,12000);
        allrate=[allrate;rate];
    end
    meanrate=mean(allrate,1);
    flash_all = [flash_all;meanrate];
    
end


chrp_all = [];
for g=1:length(chrp_tested)
    
    currUnit=info{chrp_tested(g),2};
    currDate=info{chrp_tested(g),1};
    
    startF=info{chrp_tested(g),4};
    stopF=info{chrp_tested(g),5};
    
    %pathes%
    dpath=fullfile(path2,'data',currDate);
    upath=fullfile(dpath,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    %load current unit
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,strcat('unit_',currUnit)))
            found=1;
        end
    end
    
    load(fullfile(upath,ulist(m).name))
    
    %find files with chirp stimulus
    heka=unit{1,2}(:,1);
    flashfiles=[];
    for h=1:length(heka)
        if ~isempty(regexp(heka{h},'chirp'))
            flashfiles=[flashfiles;h];
        end
    end
    tmp=find(flashfiles>startF);
    flashfiles=flashfiles(tmp);
    tmp=find(flashfiles<=stopF);
    flashfiles=flashfiles(tmp);
    
    allrate=[];
    for f=1:length(flashfiles)
        spikes=unit{1,2}{flashfiles(f),2};
        rate=MEA_spikerates(round(spikes),50,25000);
        allrate=[allrate;rate];
    end
    meanrate=mean(allrate,1);
    chrp_all = [chrp_all;meanrate];
    
end
%% collect data for responsive units
int = intersect(chrp,flash);
comb_spike = [];
flash_corr = [];
flash_corr_trans = [];

for i = 1:length(chrp)
    %     tmp = find(flash==int(i));
    %     f = flash_spike(tmp,:);
    %     f = [f(7750:10000) f(4000:6250)]; %250ms before white, white, black, +250ms
    tmp = find(chrp==chrp(i));
    c = chrp_spike(tmp,:);
    c = [c(3020:10970) c(11450:20000)];%end plus 500ms
    c2 = medfilt1(c,200); c2(1:300) = c2(301); c2(end-299) = c2(end-300);
    %     comb_spike = [comb_spike;f c];
    comb_spike = [comb_spike;c];
    
    if ~isempty(find(flash==chrp(i)))
        tmp = find(flash==chrp(i));
        f = flash_spike(tmp,:);
        %     f = [f(7750:9950) f(4000:5950)]; %250ms before white, white, black, +250ms
        flash_corr = [flash_corr;f(1:12100)];
        flash_corr_trans = [flash_corr_trans;flash_trans(tmp)];
    else
        flash_corr = [flash_corr;zeros(1,12100)];
        flash_corr_trans = [flash_corr_trans;-1];
        
    end
end

comb_all = [];
for i = 1:length(chrp_tested)
    %     tmp = find(flash==int(i));
    %     f = flash_spike(tmp,:);
    %     f = [f(7750:10000) f(4000:6250)]; %250ms before white, white, black, +250ms
    tmp = find(chrp_tested==chrp_tested(i));
    c = chrp_all(tmp,:);
    c = [c(3020:10970) c(11450:20000)];%end plus 500ms
    c2 = medfilt1(c,200); c2(1:300) = c2(301); c2(end-299) = c2(end-300);
    %     comb_spike = [comb_spike;f c];
    comb_all = [comb_all;c];
    
    
end

% for i = 1:length(int)
%     tmp = find(flash==int(i));
%     f = flash_spike(tmp,:);
%     f = [f(7750:10000) f(4000:6250)]; %250ms before white, white, black, +250ms
%     tmp = find(chrp==int(i));
%     c = chrp_spike(tmp,:);
%     c = [c(3020:10970) c(11450:20000)];%end plus 500ms
%     c2 = medfilt1(c,200); c2(1:300) = c2(301); c2(end-299) = c2(end-300);
% %     comb_spike = [comb_spike;f c];
%     comb_spike = [comb_spike;c];
% end
%% adjut Cowan 2020 to the same temporal resolution as our data and cut relevant parts
load(fullfile(path2,'Cowan2020'))

coords2 = [];
mi = 100; ma = 0;
for c = 1:6
    mi = min(mi,coords{c,1}(1));
    ma = max(ma,coords{c,1}(end));
end
mi = floor(mi); ma = ceil(ma);
ms = (ma-mi)/twosec*2*1000;
xq = mi:(ma-mi)/ms:ma;
for c = 1:6
    currx = coords{c,1};
    [s1,s2] = sort(currx);
    currx = s1;
    curry = coords{c,2}; curry = curry(s2);
    [u1 u2] = unique(currx); currx = u1;
    curry = curry(u2);
    vq1 = interp1(currx,curry,xq);
    tmp = find(isnan(vq1(1:2000)));
    vq1(tmp) = vq1(tmp(end)+1);
    tmp = find(isnan(vq1(25000:end)));
    vq1(tmp) = vq1(tmp(1)-1);
    
    vq2 = medfilt1(vq1(9700:24370),600);
    vq1(9710:24360) = vq2(11:end-10);
    coords2 = [coords2;300-vq1];
    
end
%
co = [  coords2(:,8446:8446+250) coords2(:,9912:24360+ 250)];
co2 = [];
for c = 1:size(co,1)
    now = co(c,:);
    new = resample(now,16522,length(now));
    co2 = [co2;new(11:16512)];
end

% coords3 = [coords2(:,1322-249:1322+2000) coords2(:,4939:4939+2000)   coords2(:,8503:8503+250) ...
%     co2];

coords3 = [ co2];


%% normalize templates and our data for matching
templ = [];
for c = 1:5
    now = coords3(c,:);
    now(1:220) = now(221);
    now = (now-min(now))./(max(now)-min(now));
    sub = now(8468:9005);
    templ = [templ;now];
end

ours = []; ours2 = [];
for c = 1:size(comb_spike,1)
    now = comb_spike(c,:);
    [yupper,ylower] = envelope(now,200,'peak');
    now2 = movmedian(now,600);
    now2 = (now2-min(now2))./(max(now2)-min(now2));
    sub = now2(8468:9005);
    ours = [ours;now2];
    ours2 = [ours2;yupper];
end
%% find cells with flash + chirp responses and collect polarity & transience
vers = 2;

S = sum(flash_corr');
IDboth = find(S~=0);
flash_type = [];
for ii = 1:length(IDboth)
    curr = flash_corr(IDboth(ii),:);
    tr = flash_corr_trans(IDboth(ii),:);
    bckgr = mean(curr);
    st = std(curr);
    [m(1),~] = max(curr(1950:3950));%b1
    [m(2),~] = max(curr(3950:5950));%w1
    [m(3),~] = max(curr(5950:7950));%w2
    [m(4),~] = max(curr(7950:9950));%b2
    tmp = find(m>bckgr+3*st);
    if length(tmp)>2
        flash_type = [flash_type;2 tr];
    elseif ~isempty(find(tmp==1)) ||~isempty(find(tmp==4))
        if  ~isempty(find(tmp==2)) ||~isempty(find(tmp==3))
            flash_type = [flash_type;2 tr];
        else
            flash_type = [flash_type;-1 tr];
        end
    else
        flash_type = [flash_type;1 tr];
    end
    
end
% --- manual corrections ----
flash_type(11,1) = -1;

% separate into on trans, on sust, off trans, off sust and bistrat based on
% flash
comp = find(flash_type(:,1)==1);
tmp = find(flash_type(comp,2)<400);
ontr = comp(tmp);
comp = find(flash_type(:,1)==1);
tmp = find(flash_type(comp,2)>=400);
onsu = comp(tmp);
comp = find(flash_type(:,1)==-1);
tmp = find(flash_type(comp,2)<400);
oftr = comp(tmp);
comp = find(flash_type(:,1)==-1);
tmp = find(flash_type(comp,2)>=400);
ofsu = comp(tmp);
comp = find(flash_type(:,1)==2);
tmp = find(flash_type(comp,2)<400);
onof = comp(tmp);

matches = [];
for d = 1:size(ours,1)
    d1 = ours(d,:);
    D = templ-repmat(d1,5,1);
    DS = sum(abs(D)')./size(d1,2);
    matches = [matches;DS min(DS)];
end
if vers == 1
    % find similarity to the 5 templates
    
    
    matches2 = matches(IDboth,:);
    
    ONTR = matches2(ontr,2);
    ONSU = matches2(onsu,1);
    OFTR = matches2(oftr,4);
    OFSU = matches2(ofsu,5);
    ONOF = matches2(onof,3);
    
    % only consider cells that have <0.25 average difference
    tmp = find(ONSU<0.25);
    onsu2 = onsu(tmp);
    tmp = find(ONTR<0.25);
    ontr2 = ontr(tmp);
    tmp = find(OFTR<0.25);
    oftr2 = oftr(tmp);
    tmp = find(OFSU<0.25);
    ofsu2 = ofsu(tmp);
    tmp = find(ONOF<0.25);
    onof2 = onof(tmp);
else
    early = median(ours2(IDboth(ontr),500:2000),2);
    late = median(ours2(IDboth(ontr),5000:7000),2);
    di = early-late;
    tmp = find(di<0);
    ontr2 = ontr(tmp);
    early = median(ours2(IDboth(onsu),500:2000),2);
    late = median(ours2(IDboth(onsu),5000:7000),2);
    di = early-late;
    tmp = find(di<0);
    onsu2 = onsu(tmp);
    early = median(ours2(IDboth(onof),500:2000),2);
    late = median(ours2(IDboth(onof),5000:7000),2);
    di = early-late;
    tmp = find(di<0);
    onof2 = onof(tmp);
    early = median(ours2(IDboth(oftr),500:2000),2);
    late = median(ours2(IDboth(oftr),5000:7000),2);
    di = early-late;
    tmp = find(di<0);
    oftr2 = oftr(tmp);
    %     early = median(ours2(ofsu,500:2000),2);
    %     late = median(ours2(ofsu,5000:7000),2);
    %     di = early-late;
    %     tmp = find(di<0);
    %     ofsu2 = ofsu(tmp);
    ofsu2 = ofsu;
end

ALL_CHRP = intersect(R,chrp_tested);
ASSIGN = zeros(length(ALL_CHRP),1);
ASSIGN(chrp(IDboth(onsu2))) = 1;
ASSIGN(chrp(IDboth(ontr2))) = 2;
ASSIGN(chrp(IDboth(onof2))) = 3;
ASSIGN(chrp(IDboth(oftr2))) = 4;
ASSIGN(chrp(IDboth(ofsu2))) = 5;

IDnoflash =  find(S==0);
matchesO = matches(IDnoflash,:);
tmp = find(matchesO(:,6)<0.25);
IDnoflash2 = IDnoflash(tmp);
matchesO2 = matchesO(tmp,:);
IDnoflashNoTop = IDnoflash(find(matchesO(:,6)>=0.25));
IDnoflashNoTop = setdiff(IDnoflashNoTop,78);
%!!!!!----------- !!! this is not incorporated yet, but remaining into additional cluster

IDflashNotAs = setdiff(IDboth,[IDboth(onsu2) IDboth(ontr2) IDboth(onof2) IDboth(oftr2) IDboth(ofsu2)]);
flashNotAs = zeros(length(IDflashNotAs),1);
for ii = 1:length(IDflashNotAs)
    tmp = find(IDboth==IDflashNotAs(ii));
    if flash_type(tmp,1)==1
        ASSIGN(chrp(IDflashNotAs(ii))) = 6;%new ON
        flashNotAs(ii) = 6;
    elseif flash_type(tmp,1)==-1
        ASSIGN(chrp(IDflashNotAs(ii))) = 7; % new OFF
        flashNotAs(ii) = 7;
    elseif flash_type(tmp,1)==2
        ASSIGN(chrp(IDflashNotAs(ii))) = 8; % new ON-OFF
        flashNotAs(ii) = 8;
    end
end



% for ii = 1:length(IDflashNotAs)
%    figure('position',[235.4000 407.4000 1.0496e+03 200.8000])
%    subplot(1,2,1)
%    plot(flash_corr(IDflashNotAs(ii),:))
%    subplot(1,2,2)
%    plot(comb_spike(IDflashNotAs(ii),:))
%     w = waitforbuttonpress
%     close
% end
% 278 cells were tested for chrp and responded to any stimulus
% intersect(R,chrp_tested) gives the list of chrp tested & responsive in
% general
%% plot results
load(fullfile(path2,'h_normPeaks'))
tot = 278;
alph = 0.35;
[m1,m2] = min(matchesO2');

figure
subplot(6,3,2)
data = comb_spike(IDboth(ontr2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(2,:),'-','color',[0 0 1 alph])
hold on
plot(mean(data2,1),'-k')
axis tight
box off
set(gca,'xtick',[])

subplot(6,3,1)
data = flash_corr(IDboth(ontr2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(mean(data2(:,200:5900),1),'-k')
axis tight
box off
set(gca,'xtick',[])
title(['putative ON parasol n = ',int2str(length(ontr2)),' / ',num2str(round(100/tot*length(ontr2)*10)/10),'%'])


subplot(6,3,3)
idnow = find(m2==2);
ASSIGN(chrp(IDnoflash(idnow)))=-2;
data = comb_spike(IDnoflash(idnow),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(2,:),'-','color',[0 0 1 alph])
hold on
plot(mean(data2,1),'-k')
axis tight
box off
set(gca,'xtick',[])
% title(['putative ON midget (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])


subplot(6,3,5)
data = comb_spike(IDboth(oftr2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(4,:),'-','color',[0 1 0.3 alph])
hold on
plot(mean(data2),'-k')
axis tight
box off
set(gca,'xtick',[])

subplot(6,3,4)
data = flash_corr(IDboth(oftr2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(mean(data2(:,200:5900)),'-k')
axis tight
box off
set(gca,'xtick',[])
title(['putative OFF parasol n = ',int2str(length(oftr2)),' / ',num2str(round(100/tot*length(oftr2)*10)/10),'%'])


subplot(6,3,6)
idnow = find(m2==4);
ASSIGN(chrp(IDnoflash(idnow)))=-4;
data = comb_spike(IDnoflash(idnow),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(4,:),'-','color',[0 1 0.3 alph])
hold on
plot(mean(data2),'-k')
axis tight
box off
set(gca,'xtick',[])
% title(['putative ON parasol (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])


subplot(6,3,8)
data = comb_spike(IDboth(onsu2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(1,:),'-','color',[1 0.2 0.2 alph])
hold on
plot(mean(data2,1),'-k')
axis tight
box off
set(gca,'xtick',[])

subplot(6,3,7)
data = flash_corr(IDboth(onsu2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(mean(data2(:,200:5900)),'-k')
axis tight
box off
set(gca,'xtick',[])
title(['putative ON midget n = ',int2str(length(onsu2)),' / ',num2str(round(100/tot*length(onsu2)*10)/10),'%'])


subplot(6,3,9)
idnow = find(m2==1);
ASSIGN(chrp(IDnoflash(idnow)))=-1;
data = comb_spike(IDnoflash(idnow),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(1,:),'-','color',[1 0.2 0.2 alph])
hold on
plot(mean(data2),'-k')
axis tight
box off
set(gca,'xtick',[])
% title(['putative ON parasol (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])


subplot(6,3,11)
data = comb_spike(IDboth(ofsu2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(5,:),'-','color',[1 0 1 alph])
hold on
plot(mean(data2,1),'-k')
axis tight
box off
set(gca,'xtick',[])

subplot(6,3,10)
data = flash_corr(IDboth(ofsu2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(mean(data2(:,200:5900),1),'-k')
axis tight
box off
set(gca,'xtick',[])
title(['putative OFF midget n = ',int2str(length(ofsu2)),' / ',num2str(round(100/tot*length(ofsu2)*10)/10),'%'])


subplot(6,3,12)
idnow = find(m2==5);
ASSIGN(chrp(IDnoflash(idnow)))=-5;
data = comb_spike(IDnoflash(idnow),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(5,:),'-','color',[1 0 1 alph])
hold on
plot(mean(data2),'-k')
axis tight
box off
set(gca,'xtick',[])
% title(['putative OFF parasol (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])



subplot(6,3,14)
data = comb_spike(IDboth(onof2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(4,:),'-','color',[1 0.4 0 alph])
hold on
plot(mean(data2),'-k')
axis tight
box off
set(gca,'xtick',0:2000:16000,'xticklabel',0:2:16)
xlabel('sec')

subplot(6,3,13)
data = flash_corr(IDboth(onof2),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(mean(data2(:,200:5900)),'-k')
set(gca,'xtick',0:1000:5000,'xticklabel',0:1:5)
xlabel('sec')
axis tight
box off
title(['putative ON-OFF bistratified n = ',int2str(length(onof2)),' / ',num2str(round(100/tot*length(onof2)*10)/10),'%'])


subplot(6,3,15)
idnow = find(m2==4);
ASSIGN(chrp(IDnoflash(idnow)))=-4;
data = comb_spike(IDnoflash(idnow),:);
data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
plot(templ(4,:),'-','color',[1 0.4 0 alph])
hold on
plot(mean(data2,1),'-k')
axis tight
set(gca,'xtick',[])
title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
box off



subplot(6,3,16)
plot([1 1984-200],[128 128],'-','color',[0.4 0.4 0.4],'linewidth',2)
hold on
plot([1984-200 1984-200],[128 0],'-','color',[0.4 0.4 0.4],'linewidth',2)
plot([1984-200 3968-200],[0 0],'-','color',[0.4 0.4 0.4],'linewidth',2)
plot([3968-200 3968-200],[0 128],'-','color',[0.4 0.4 0.4],'linewidth',2)
plot([3968-200 5952-200],[128 128],'-','color',[0.4 0.4 0.4],'linewidth',2)
axis([ 1 5952-200 0 255])
box off
set(gca,'xtick',[])
ylabel('RGB')

subplot(6,3,17)
plot([prot2(3:481,4); prot2(485:970,4)],'-','color',[0.4 0.4 0.4])
axis tight
box off
axis([ 1 965 0 255])
set(gca,'xtick',[])
ylabel('RGB')

subplot(6,3,18)
plot([prot2(3:481,4); prot2(485:970,4)],'-','color',[0.4 0.4 0.4])
axis tight
box off
axis([ 1 965 0 255])
set(gca,'xtick',[])
ylabel('RGB')


% ---- add plotting spatio-temporal maps ---
%% how about cells that do not respond to chirp
load(fullfile(path2,'EXAMPLE_INFO'),'SUPER_COLL','CLUSTER_INFO')

IDnoChirp = setdiff(chrp_tested,chrp);
IDnoChirp = intersect(R,IDnoChirp);

% %%
% ii = 3;
% tmp2 = find(flash==tmp(ii));
% figure
% plot(flash_spike(tmp2,:))

%
% --- 125 cells respond NOT to chirp but to DG ---
normP = []; ids = [];
for n = 1:length(IDnoChirp)
    tmp = find(togo==IDnoChirp(n));
    ids = [ids;tmp];
    now = normPeaks(tmp,:,2);
    normP = [normP;now(:)' ];
end

% CLUST = [];
% EVAL = [];
% for k = 2:30
%     [idx,C] = kmeans(normP,k,'Replicates',1000);
%     eva1 = evalclusters(normP,idx,'CalinskiHarabasz');
%     eva2 = evalclusters(normP,idx,'DaviesBouldin');
%     CLUST = [CLUST idx];
%     evals = [eva1.CriterionValues; eva2.CriterionValues];
%     EVAL = [EVAL evals];
% end
% save(fullfile(path2,'cluster_DG_noMPB'),'CLUST','EVAL','normP','ids')
load(fullfile(path2,'cluster_DG_noMPB'))
CC = 6;
CL = CLUST(:,CC-1);
figure
HEATMAP = []; NoC = [];
for c = 1:CC
    tmp  =find(CL==c);
    NoC = [NoC;length(tmp)];
    htmp = 1-normP(tmp,:); htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    
    HEATMAP = [HEATMAP htmp];
    if c < CC
        HEATMAP = [HEATMAP ones(4,1)];
    end
end

imagesc(HEATMAP)
colormap(gca,'gray')
axis off
box off
%
load(fullfile(path2,'abs_ID_exampleCells'))
load(fullfile(path2,'EXAMPLE_INFO'))
% fig5_clust = []; clust_new = [];
% for ii = 1:length(COLL_ID)
%     inow = COLL_ID(ii);
%
%            clust_new = [clust_new;ASSIGN(inow)];
%
% %     pos = find(SUPER_COLL(:,1)==inow);
% %     cl = SUPER_COLL(pos,2);
% %     fig5_clust = [fig5_clust;cl];
% end

%% main plot
% yh1 = 0.025; yh2 = 0.028; yh3 = 0.023;
% col = [0 0 1 0.35;0 0.2 0.8 0.35;1 0.2 0 0.35;1 0 1 0.35;1 0.4 0 0.35];
% figure
% for ax = 1:5
%     if ax == 1
%         IDX = ontr2;
%         tm = 2;
%     elseif ax == 2
%         IDX = oftr2;
%         tm = 4;
%     elseif ax ==3
%         IDX = onsu2;
%         tm = 1;
%     elseif ax == 4
%         IDX = ofsu2;
%         tm = 5;
%     elseif ax ==5
%         IDX = onof2;
%         tm = 3;
%     end
%     axes('position',[0.08 0.95-yh2*(2*(ax))+yh1 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     data = flash_corr(IDboth(IDX),:);
%     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
%     plot(mean(data2(:,200:5900),1),'-k')
%     axis tight
%     box off
%     axis off
%     set(gca,'xtick',[])
%     
%     
%     axes('position',[0.08 0.95-yh2*(2*(ax)) 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     imagesc(1-data2(:,200:5900))
%     colormap gray
%     axis off
%     
%     axes('position',[0.2 0.95-yh2*(2*(ax))+yh1 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     data = comb_spike(IDboth(IDX),:);
%     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
%     if ~isempty(IDX)
%         plot(templ(tm,:),'-','color',col(ax,:))
%         hold on
%         plot(mean(data2,1),'-k')
%     end
%     axis tight
%     axis off
%     box off
%     set(gca,'xtick',[])
%     
%     axes('position',[0.2 0.95-yh2*(2*(ax)) 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     imagesc(1-data2)
%     colormap gray
%     axis off
%     
%     
%     axes('position',[0.32 0.95-yh2*(2*(ax))+yh1 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     idnow = find(m2==tm);
%     ASSIGN(chrp(IDnoflash(idnow)))=-tm;
%     data = comb_spike(IDnoflash(idnow),:);
%     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
%     if ~isempty(idnow)
%         plot(templ(tm,:),'-','color',col(ax,:))
%         hold on
%         plot(mean(data2,1),'-k')
%     end
%     axis tight
%     box off
%     axis off
%     set(gca,'xtick',[])
%     % title(['putative ON midget (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
%     %     title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'],...
%     %         'fontname','Arial','fontsize',9)
%     
%     
%     axes('position',[0.32 0.95-yh2*(2*(ax)) 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     imagesc(1-data2)
%     colormap gray
%     axis off
%     
%     axes('position',[0.02 0.95-yh2*(2*(ax)) 0.035 yh3*2])%[0.08 0.9 0.8 0.02])
%     if ax == 1
%         text(0.05,0.9,'cluster 1','fontname','Arial','fontsize',7,'fontweight','bold')
%         text(0.05,0.7,'put. ON parasol','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['put. ON parasol, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%     elseif ax == 2
%         text(0.05,0.9,'cluster 2','fontname','Arial','fontsize',7,'fontweight','bold')
%         text(0.05,0.7,'put. OFF parasol','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['put. OFF parasol, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%     elseif ax ==3
%         text(0.05,0.9,'cluster 3','fontname','Arial','fontsize',7,'fontweight','bold')
%         text(0.05,0.7,'put. ON midget','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['putative ON midget, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%         
%     elseif ax ==4
%         text(0.05,0.9,'cluster 4','fontname','Arial','fontsize',7,'fontweight','bold')
%         
%         text(0.05,0.7,'put. OFF midget','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['putative OFF midget, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%         
%     elseif ax ==5
%         text(0.05,0.9,'cluster 5','fontname','Arial','fontsize',7,'fontweight','bold')
%         
%         text(0.05,0.65,'put. ON-OFF bistrat.','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['putative ON-OFF bistratified, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%         
%     end
%     text(0.05,0.4,['n = ',int2str(length(IDX)),...
%         ' / ',int2str(length(idnow))],'fontname','Arial','fontsize',7)
%     %         text(0.05,0.3,['n_n_f = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'],...
%     %         'fontname','Arial','fontsize',7)
%     text(0.05,0.15,['total n = ',int2str(length(idnow)+length(IDX)),' / ',num2str(round(100/tot*(length(idnow)+length(IDX))*10)/10),'%'],...
%         'fontname','Arial','fontsize',7)
%     box off
%     axis off
%     
%     axes('position',[0.44 0.95-yh2*(2*(ax)) 0.06 yh3*2])%[0.08 0.9 0.8 0.02])
%     normP = [];
%     for n = 1:length(IDX)
%         tmp = find(togo == chrp(IDboth(IDX(n))));
%         now = normPeaks(tmp,:,2);
%         normP = [normP;now(:)' ];
%     end
%     for n = 1:length(idnow)
%         tmp = find(togo == chrp(IDnoflash(idnow(n))));
%         now = normPeaks(tmp,:,2);
%         normP = [normP;now(:)' ];
%     end
%     htmp = 1-normP; htmp = mean(htmp,1);
%     htmp=reshape(htmp,4,6);
%     htmp = flipud(fliplr(htmp));
%     
%     imagesc(htmp)
%     colormap gray
%     axis tight
%     axis off
%     set(gca,'plotboxaspectratio',[1 1/6*4 1])
% end
% 
% % --- cells with flash+chrp but not assigned to the top 5
% for ax = 6:8
%     IDX = IDflashNotAs(find(flashNotAs==ax));
%     
%     axes('position',[0.08 0.9-yh2*(2*(ax))+yh1 0.1  yh3])%[0.08 0.9 0.8 0.02])
%     data = flash_corr(IDX,:);
%     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
%     plot(mean(data2(:,200:5900),1),'-k')
%     axis tight
%     box off
%     axis off
%     set(gca,'xtick',[])
%     
%     
%     axes('position',[0.08 0.9-yh2*(2*(ax)) 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     imagesc(1-data2(:,200:5900))
%     colormap gray
%     axis off
%     
%     axes('position',[0.2 0.9-yh2*(2*(ax))+yh1 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     data = comb_spike(IDX,:);
%     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
%     
%     plot(mean(data2,1),'-k')
%     
%     axis tight
%     axis off
%     box off
%     set(gca,'xtick',[])
%     
%     axes('position',[0.2 0.9-yh2*(2*(ax)) 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     imagesc(1-data2)
%     colormap gray
%     axis off
%     
%     
%     
%     axes('position',[0.03 0.9-yh2*(2*(ax)) 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
%     if ax == 6
%         text(0.05,0.9,'cluster 6 (ON)','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['put. ON parasol, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%     elseif ax == 7
%         text(0.05,0.9,'cluster 7 (OFF)','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['put. OFF parasol, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%     elseif ax ==8
%         text(0.05,0.9,'cluster 8 (ON-OFF)','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['putative ON midget, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%     end
%     text(0.05,0.7,['n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'],...
%         'fontname','Arial','fontsize',7)
%     
%     box off
%     axis off
%     
%     axes('position',[0.44 0.9-yh2*(2*(ax)) 0.06 yh3*2])%[0.08 0.9 0.8 0.02])
%     normP = [];
%     for n = 1:length(IDX)
%         tmp = find(togo == chrp(IDX(n)));
%         now = normPeaks(tmp,:,2);
%         normP = [normP;now(:)' ];
%     end
%     htmp = 1-normP; htmp = mean(htmp,1);
%     htmp=reshape(htmp,4,6);
%     htmp = flipud(fliplr(htmp));
%     
%     imagesc(htmp)
%     colormap gray
%     axis tight
%     axis off
%     set(gca,'plotboxaspectratio',[1 1/6*4 1])
% end
% 
% % --- cells without chrp response, but tested
% for ax = 9:9+CC-1
%     idx = IDnoChirp(find(CLUST(:,CC-1)==ax-8));
%     IDX = []; IDXf = [];
%     for ii =1:length(idx)
%         tmp = find(chrp_tested==idx(ii));
%         IDX = [IDX;tmp];
%         tmp = find(flash==idx(ii));
%         IDXf = [IDXf;tmp];
%     end
%     
%     
%     axes('position',[0.08 0.85-yh2*(2*(ax))+yh1 0.1  yh3])%[0.08 0.9 0.8 0.02])
%     data = flash_corr(IDXf,:);
%     SS = sum(data,2); tmp = find(SS~=0);
%     data = data(tmp,:);
%     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
%     plot(nanmean(data2(:,200:5900),1),'-k')
%     axis tight
%     box off
%     axis off
%     set(gca,'xtick',[])
%     
%     
%     axes('position',[0.08 0.85-yh2*(2*(ax)) 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     imagesc(1-data2(:,200:5900))
%     colormap gray
%     axis off
%     
%     axes('position',[0.2 0.85-yh2*(2*(ax))+yh1 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     data = comb_all(IDX,:);
%     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
%     plot(mean(data2,1),'-k')
%     
%     axis tight
%     axis off
%     box off
%     set(gca,'xtick',[])
%     
%     axes('position',[0.2 0.85-yh2*(2*(ax)) 0.1 yh3])%[0.08 0.9 0.8 0.02])
%     imagesc(1-data2)
%     colormap gray
%     axis off
%     
%     
%     
%     axes('position',[0.03 0.85-yh2*(2*(ax)) 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
%     if ax == 9
%         text(0.05,0.9,'cluster 9 (ON)','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['put. ON parasol, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%     elseif ax == 10
%         text(0.05,0.9,'cluster 10 (OFF)','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['put. OFF parasol, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%     elseif ax ==11
%         text(0.05,0.9,'cluster 11 (ON-OFF)','fontname','Arial','fontsize',7,'fontweight','bold')
%         %         title(['putative ON midget, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
%     elseif ax == 12
%         text(0.05,0.9,'cluster 12 (ON-OFF)','fontname','Arial','fontsize',7,'fontweight','bold')
%         
%     elseif ax == 13
%         text(0.05,0.9,'cluster 13 (ON-OFF)','fontname','Arial','fontsize',7,'fontweight','bold')
%         
%     elseif ax == 14
%         text(0.05,0.9,'cluster 14 (ON-OFF)','fontname','Arial','fontsize',7,'fontweight','bold')
%         
%     end
%     text(0.05,0.7,['n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'],...
%         'fontname','Arial','fontsize',7)
%     
%     box off
%     axis off
%     
%     axes('position',[0.44 0.85-yh2*(2*(ax)) 0.06 yh3*2])%[0.08 0.9 0.8 0.02])
%     normP = [];
%     for n = 1:length(IDX)
%         tmp = find(togo == idx(n));
%         now = normPeaks(tmp,:,2);
%         normP = [normP;now(:)' ];
%     end
%     htmp = 1-normP; htmp = mean(htmp,1);
%     htmp=reshape(htmp,4,6);
%     htmp = flipud(fliplr(htmp));
%     
%     imagesc(htmp)
%     colormap gray
%     axis tight
%     axis off
%     set(gca,'plotboxaspectratio',[1 1/6*4 1])
% end
%% main plot v2
% yh1 = 0.025; yh2 = 0.03; yh3 = 0.023;
% yst1 = 0.97; yst2 = 0.94; yst3 = 0.91;
yh1 = 0.022; yh2 = 0.024; yh3 = 0.018;
yst1 = 0.97; yst2 = 0.94; yst3 = 0.91;
x1 = 0.16; x2 = 0.5; x3 = 0.84; x4 = 0.9;
w1 = 0.28; w2 = 0.06;
col = [0 0 1 0.35;0 0.2 0.8 0.35;1 0.2 0 0.35;1 0 1 0.35;1 0.4 0 0.35];
ex_names = {'A_{1}','A_{2}','A_{3}','A_{4}','A_{5}',...
    'B_{1}','B_{2}','B_{3}','B_{4}','B_{5}','B_{6}','B_{7}','C_{1}','C_{2}','C_{3}'};


figure('position',[310.6000 924.2000 760.0000 1028],'color','w')
mxv = []; mxv2 =[];
for ax = 1:5
    if ax == 1
        IDX = ontr2;
        tm = 2;
    elseif ax == 2
        IDX = oftr2;
        tm = 4;
    elseif ax ==3
        IDX = onsu2;
        tm = 1;
    elseif ax == 4
        IDX = ofsu2;
        tm = 5;
    elseif ax ==5
        IDX = onof2;
        tm = 3;
    end
    axes('position',[x1 yst1-yh2*(2*(ax))+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    data = flash_corr(IDboth(IDX),:);
    idnow = find(m2==tm);
    idnow = chrp(IDnoflash2(idnow));
    idnow2 = [];
    for id = 1:length(idnow)
        tmp = find(R == idnow(id));
        idnow2 = [idnow2;tmp];
    end
    data_2nd = flash_all(idnow2,1:12100);
    data1b = [data;data_2nd];
    data2=data1b;
    plot(mean(data2(:,200:5900),1),'-','color',[0.6 0.6 0.6])
    hold on
    plot(mean(data(:,200:5900),1),'-k')
    mxv = [mxv;max(mean(data(:,200:5900),1))];
    axis([1 5700 0 40])
    box off
    axis off
    set(gca,'xtick',[])
    
    
    axes('position',[x1 yst1-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2(:,200:5900))
    colormap gray
    hold on
    plot([10 10],[0.5 length(IDX)+0.5],'-','color',col(tm,1:3),'linewidth',5)
    axis off
    
    axes('position',[x2 yst1-yh2*(2*(ax))+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    data = comb_spike(IDboth(IDX),:);
    idnow = find(m2==tm);
    ASSIGN(chrp(IDnoflash2(idnow)))=-tm;
    data_2nd = comb_spike(IDnoflash2(idnow),:);
    data = [data;data_2nd];
    %     data2 = data;
    data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    if ~isempty([IDX;idnow'])
        %     plot(templ(tm,:),'-','color',col(ax,:))
        hold on
        plot(mean(data2,1),'-k')
        mxv2 = [mxv2;max(mean(data2,1))];
    end
    %     axis([1 16502 0 30])
    axis tight
    axis off
    box off
    set(gca,'xtick',[])
    
    axes('position',[x2 yst1-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2)
    colormap gray
    axis off
    hold on
    plot([10 10],[0.5 length(IDX)+0.5],'-','color',col(tm,1:3),'linewidth',5)
    
    
    
    
    axes('position',[0.02 yst1-yh2*(2*(ax)) 0.035 yh3*2])%[0.08 0.9 0.8 0.02])
    if ax == 1
        text(0.05,0.9,'cluster 1','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'put. ON parasol','fontname','Arial','fontsize',7,'fontweight','bold')
        %         title(['put. ON parasol, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
    elseif ax == 2
        text(0.05,0.9,'cluster 2','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'put. OFF parasol','fontname','Arial','fontsize',7,'fontweight','bold')
        %         title(['put. OFF parasol, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
    elseif ax ==3
        text(0.05,0.9,'cluster 3','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'put. ON midget','fontname','Arial','fontsize',7,'fontweight','bold')
        %         title(['putative ON midget, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
        
    elseif ax ==4
        text(0.05,0.9,'cluster 4','fontname','Arial','fontsize',7,'fontweight','bold')
        
        text(0.05,0.65,'put. OFF midget','fontname','Arial','fontsize',7,'fontweight','bold')
        %         title(['putative OFF midget, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
        
    elseif ax ==5
        text(0.05,0.9,'cluster 5','fontname','Arial','fontsize',7,'fontweight','bold')
        
        text(0.05,0.65,'put. ON-OFF bistrat.','fontname','Arial','fontsize',7,'fontweight','bold')
        %         title(['putative ON-OFF bistratified, n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'])
        
    end
    text(0.05,0.4,['n = ',int2str(length(IDX)),...
        ', ',int2str(length(idnow))],'fontname','Arial','fontsize',7)
    %         text(0.05,0.3,['n_n_f = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'],...
    %         'fontname','Arial','fontsize',7)
    text(0.05,0.15,['total n = ',int2str(length(idnow)+length(IDX)),' / ',num2str(round(100/tot*(length(idnow)+length(IDX))*10)/10),'%'],...
        'fontname','Arial','fontsize',7)
    box off
    axis off
    
    axes('position',[x3 yst1-yh2*(2*(ax)) w2 yh3*2])%[0.08 0.9 0.8 0.02])
    normP = [];
    for n = 1:length(IDX)
        tmp = find(togo == chrp(IDboth(IDX(n))));
        now = normPeaks(tmp,:,2);
        normP = [normP;now(:)' ];
    end
    for n = 1:length(idnow)
        tmp = find(togo == chrp(IDnoflash(idnow(n))));
        now = normPeaks(tmp,:,2);
        normP = [normP;now(:)' ];
    end
    htmp = 1-normP; htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    
    if ax == 1
     axes('position',[0.015 yst1-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'A_1','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
      axes('position',[x2-0.02 yst1-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'A_2','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
      axes('position',[x3-0.02 yst1-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'A_3','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
    end
end

% --- cells with flash+chrp but not assigned to the top 5
for ax = 6:7
    IDX = IDflashNotAs(find(flashNotAs==ax));
    ASSIGN(chrp(IDX)) = ax;
    
    axes('position',[x1 yst2-yh2*(2*(ax))+yh1 w1  yh3])%[0.08 0.9 0.8 0.02])
    data = flash_corr(IDX,:);
%     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
data2 = data;    
plot(mean(data2(:,200:5900),1),'-k')
    axis([1 5700 0 40])
    box off
    axis off
    set(gca,'xtick',[])
    
    
    axes('position',[x1 yst2-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2(:,200:5900))
    colormap gray
    axis off
    
    axes('position',[x2 yst2-yh2*(2*(ax))+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    data = comb_spike(IDX,:);
    data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    
    plot(mean(data2,1),'-k')
    
    axis tight
    axis off
    box off
    set(gca,'xtick',[])
    
    axes('position',[x2 yst2-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2)
    colormap gray
    axis off
    
    
    
    axes('position',[0.02 yst2-yh2*(2*(ax)) 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
    if ax == 6
        text(0.05,0.9,'cluster 6','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'ON full-field','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 7
        text(0.05,0.9,'cluster 7','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'OFF full-field','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax ==8
        text(0.05,0.9,'cluster 8','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'ON-OFF full-field','fontname','Arial','fontsize',7,'fontweight','bold')
    end
    text(0.05,0.4,['n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'],...
        'fontname','Arial','fontsize',7)
    
    box off
    axis off
    
    axes('position',[x3 yst2-yh2*(2*(ax)) w2 yh3*2])%[0.08 0.9 0.8 0.02])
    normP = [];
    for n = 1:length(IDX)
        tmp = find(togo == chrp(IDX(n)));
        now = normPeaks(tmp,:,2);
        normP = [normP;now(:)' ];
    end
    htmp = 1-normP; htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    
      if ax == 6
     axes('position',[0.015 yst2-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'B_1','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
      axes('position',[x2-0.02 yst2-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'B_2','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
      axes('position',[x3-0.02 yst2-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'B_3','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
    end
end

% --- cells without flash, but chrp that was not assigned to the top 5
for ax = 8
    IDX = chrp(IDnoflashNoTop);
     ASSIGN(IDX) = ax;
    
    axes('position',[x1 yst2-yh2*(2*(ax))+yh1 w1  yh3])%[0.08 0.9 0.8 0.02])
    idnow2 = [];
    for id = 1:length(IDX)
        tmp = find(R == IDX(id));
        idnow2 = [idnow2;tmp];
    end
    data = flash_all(idnow2,1:12100);
    %     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    plot(mean(data(:,200:5900),1),'-','color',[0.6 0.6 0.6])
    axis([1 5700 0 40])
    box off
    axis off
    set(gca,'xtick',[])
    
    
    axes('position',[x1 yst2-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data(:,200:5900))
    colormap gray
    axis off
    
    axes('position',[x2 yst2-yh2*(2*(ax))+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    data = comb_spike(IDnoflashNoTop,:);
    data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    
    plot(mean(data2,1),'-k')
    
    axis tight
    axis off
    box off
    set(gca,'xtick',[])
    
    axes('position',[x2 yst2-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2)
    colormap gray
    axis off
    
    
    
    axes('position',[0.02 yst2-yh2*(2*(ax)) 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
    if ax == 8
        text(0.05,0.9,'cluster 8','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'only chirp','fontname','Arial','fontsize',7,'fontweight','bold')
    end
    text(0.05,0.4,['n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'],...
        'fontname','Arial','fontsize',7)
    
    box off
    axis off
    
    axes('position',[x3 yst2-yh2*(2*(ax)) w2 yh3*2])%[0.08 0.9 0.8 0.02])
    normP = [];
    for n = 1:length(IDX)
        tmp = find(togo == IDX(n));
        now = normPeaks(tmp,:,2);
        normP = [normP;now(:)' ];
    end
    htmp = 1-normP; htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
end

% --- cells without chrp response, but tested
tmp = find(NoC>=5);
cl = 1:CC;
cl = cl(tmp); cnt = 0;
for ax = 9:9+length(cl)-1
    cnt = cnt+1;
    idx = IDnoChirp(find(CLUST(:,CC-1)==cl(cnt)));
    IDX = []; IDXf = []; idx2 = [];
    for ii =1:length(idx)
        tmp = find(chrp_tested==idx(ii));
        IDX = [IDX;tmp];
        tmp = find(flash==idx(ii));
        if ~isempty(tmp)
            IDXf = [IDXf;tmp];
        else
            idx2 = [idx2;idx(ii)];
        end
    end
    ASSIGN(idx) = ax;
    
    
    axes('position',[x1 yst3-yh2*(2*(ax))+yh1 w1  yh3])%[0.08 0.9 0.8 0.02])
    data = flash_corr(IDXf,:);
    SS = sum(data,2); tmp = find(SS~=0);
    data = data(tmp,:);
    idnow2 = [];
    for id = 1:length(idx2)
        tmp = find(R == idx2(id));
        idnow2 = [idnow2;tmp];
    end
    data_2nd = flash_all(idnow2,1:12100);
    data1b = [data;data_2nd];
    %     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    data2=data1b;
    data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    plot(nanmean(data1b(:,200:5900),1),'-','color',[0.6 0.6 0.6])
    hold on
    %     plot(nanmean(data(:,200:5900),1),'-k')
    
    axis([1 5700 0 40])
    box off
    axis off
    set(gca,'xtick',[])
    
    
    axes('position',[x1 yst3-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data1b(:,200:5900))
    colormap gray
    hold on
    plot([10 10],[0.5 size(data,1)+0.5],'-','color','k','linewidth',5)
    axis off
    
    axes('position',[x2 yst3-yh2*(2*(ax))+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    data = comb_all(IDX,:);
    data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    plot(mean(data2,1),'-k')
    axis tight
    axis off
    box off
    set(gca,'xtick',[])
    
    axes('position',[x2 yst3-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2)
    colormap gray
    axis off
    
    
    
    axes('position',[0.02 yst3-yh2*(2*(ax)) 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
    if ax == 9
        text(0.05,0.9,'cluster 9','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'speed-tuned?','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 10
        text(0.05,0.9,'cluster 10','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'temporal tuning','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax ==11
        text(0.05,0.9,'cluster 11','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'fast-big','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 12
        text(0.05,0.9,'cluster 12','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'4Hz-big','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 13
        text(0.05,0.9,'cluster 13','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'slow','fontname','Arial','fontsize',7,'fontweight','bold')
        
    elseif ax == 14
        text(0.05,0.9,'cluster 14','fontname','Arial','fontsize',7,'fontweight','bold')
        text(0.05,0.65,'medium speed/size','fontname','Arial','fontsize',7,'fontweight','bold')
    end
    text(0.05,0.4,['n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'],...
        'fontname','Arial','fontsize',7)
    
    box off
    axis off
    
    axes('position',[x3 yst3-yh2*(2*(ax)) w2 yh3*2])%[0.08 0.9 0.8 0.02])
    normP = [];
    for n = 1:length(IDX)
        tmp = find(togo == idx(n));
        now = normPeaks(tmp,:,2);
        normP = [normP;now(:)' ];
    end
    htmp = 1-normP; htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    
      if ax == 9
     axes('position',[0.015 yst3-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'C_1','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
      axes('position',[x2-0.02 yst3-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'C_2','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
      axes('position',[x3-0.02 yst3-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
     text(0.05,0.9,'C_3','fontname','Arial','fontsize',10,'fontweight','bold')
     axis off
     box off
    end
end

% axes('position',[0.31 yst3-yh2*(2*(ax)) 0.015 yst1-(yst3-yh2*(2*(ax)))])
for c = 1:ax
    tmp = find(ASSIGN(COLL_ID)==c);
    if c == 1
        newcl = 3;
    elseif c == 2
        newcl = 1;
    elseif c == 3
        newcl = 5;
    elseif c == 4
        newcl = 2;
    elseif c == 5
        newcl = 4;
    else
        newcl = c;
    end
    if c <6
        axes('position',[x3-0.048 yst1-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
    elseif c<9
        axes('position',[x3-0.048 yst2-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x3-0.048 yst3-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
    end
    
    for t = 1:length(tmp)
        if t<4
        text(0.05,1.25-0.4*t,ex_names{tmp(t)},'fontname','Arial','fontsize',7)
        elseif t == 4
                  text(0.99,1.25-0.4,ex_names{tmp(t)},'fontname','Arial','fontsize',7)  
        else
            text(0.99,1.25-0.8,ex_names{tmp(t)},'fontname','Arial','fontsize',7)  
        end
    end
    axis off
    box off
end
%%
% list = {'sust on','trans on','on-off','trans off','sust off'};
% col = 'rbygm';
% sh = 4;
% for ii = 1:length(oftr)
%     pl = IDboth(oftr(ii));
%     close all
%     figure
%     subplot(3,1,1)
%     plot((comb_spike(pl,:)-min(comb_spike(pl,:)))./(max(comb_spike(pl,:))-min(comb_spike(pl,:))),'-c')
%     hold on
%     plot(ours(pl,:),'-k')
%     
%     subplot(3,1,2)
%     plot(ours(pl,:),'-k')
%     hold on
%     for t = sh
%         plot(templ(t,:),'-','color',col(t))
%     end
%     %     tmp = find(matches(pl,1:5)==matches(pl,6));
%     % [m1,m2] = sort(matches(s2(pl),1:5));
%     % plot(templ(tmp,:),'-r')
%     % plot(templ(m2(2),:),'-c')
%     % title([num2str(s1(pl)),' || ',list{tmp},' / ',list{m2(2)}])
%     title([num2str(matches(pl,sh)),' || ',list{sh}])
%     
%     subplot(3,1,3)
%     plot(flash_corr(pl,:))
%     w = waitforbuttonpress
% end




