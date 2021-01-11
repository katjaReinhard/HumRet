clear
close all
clc

%% pathes
path1 = 'C:\Users\Katja\OneDrive - imec\HumanRetina\PlosOne';
path2 = 'C:\Users\Katja\OneDrive - imec\HumRet';
pathcode = 'C:\Users\katja\OneDrive - imec\HumRet\_scripts';
pathSave = 'C:\Users\katja\OneDrive - imec\HumanRetina\PlosOne\resubmission';
addpath(genpath(pathcode))

recluster = 0;
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
%% stimulus info
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
%% adjut Cowan 2020 to the same temporal resolution as our data and cut relevant parts
vers = 2;
if vers == 1% old
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
else
    old = load(fullfile(path2,'Cowan2020'));
    load(fullfile(pathSave,'Cowan2020'))
    coords=cell(6,2);
    coords{1,1} = c1(:,1);  coords{1,2} = c1(:,2);
    coords{2,1} = c2(:,1);  coords{2,2} = c2(:,2);
    coords{3,1} = c3(:,1);  coords{3,2} = c3(:,2);
    coords{4,1} = c4(:,1);  coords{4,2} = c4(:,2);
    coords{5,1} = c5(:,1);  coords{5,2} = c5(:,2);
    coords{6,1} = old.coords{6,1};  coords{6,2} = old.coords{6,2};
    coords2 = [];
    mi = 100; ma = 0;
    for c = 1:6
        mi = min(mi,coords{c,1}(1));
        ma = max(ma,coords{c,1}(end));
    end
    mi = floor(mi); ma = ceil(ma);
    twosec = ma/13.5;
    ms = (ma-mi)/twosec*2*1000;
    xq = mi:(ma-mi)/ms:ma;
    for c = 1:5
        currx = coords{c,1};
        [s1,s2] = sort(currx);
        currx = s1;
        curry = coords{c,2}; curry = curry(s2);
        [u1 u2] = unique(currx); currx = u1;
        curry = curry(u2);
        vq1 = resample(curry,length(xq),length(currx));
        %         vq1 = interp1(currx,curry,xq);
        %         tmp = find(isnan(vq1(1:2000)));
        %         vq1(tmp) = vq1(tmp(end)+1);
        %         tmp = find(isnan(vq1(25000:end)));
        %         vq1(tmp) = vq1(tmp(1)-1);
        
        %         vq2 = medfilt1(vq1(9700:24370),600);
        %         vq1(9710:24360) = vq2(11:end-10);
        coords2 = [coords2;vq1(1:26000)'];
        
    end
    %
    co = [  coords2(:,8446:8446+250) coords2(:,9912:24360+ 250)];
    co2 = [];
    for c = 1:size(co,1)
        now = co(c,:);
        new = resample(now,16522,length(now));
        co2 = [co2;new(11:16512)];
    end
    coords3 = [ co2];
    
end
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
[chirpFFTnorm,chirpAmp,freqSmoothChrp,togoChrp]=get_ChirpFFTnorm(path2,info);

if vers == 1
    % find similarity to the 5 templates
    
    
    matches2 = matches(IDboth,:);
    
    ONTR = matches2(ontr,2);
    ONSU = matches2(onsu,1);
    OFTR = matches2(oftr,4);
    OFSU = matches2(ofsu,5);
    ONOF = matches2(onof,3);
    
    % only consider cells that have <0.25 average difference
    tmp = find(ONSU<chrpThr);
    onsu2 = onsu(tmp);
    tmp = find(ONTR<chrpThr);
    ontr2 = ontr(tmp);
    tmp = find(OFTR<chrpThr);
    oftr2 = oftr(tmp);
    tmp = find(OFSU<chrpThr);
    ofsu2 = ofsu(tmp);
    tmp = find(ONOF<chrpThr);
    onof2 = onof(tmp);
elseif vers == 2
    early = median(ours2(IDboth(ontr),500:2000),2);
    late = median(ours2(IDboth(ontr),5000:7000),2);
    di = early-late;
    tmp = find(di<0);
    ontr2 = ontr(tmp);
    early = median(ours2(IDboth(onsu),500:2000),2);
    late = median(ours2(IDboth(onsu),5000:7000),2);
    di = early-late;
    tmp = find(abs(di)<0.2);
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
else
    early = mean(chirpFFTnorm(:,1:22),2);
late = mean(chirpFFTnorm(:,49:60),2);
di = (early-late)./(early+late);
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
% bck = mean(ours2(:,1:400));
% early = mean(abs(ours2(:,2000:4000)),2);
% late = mean(abs(ours2(:,5000:7000)),2);
early = mean(chirpFFTnorm(:,1:20),2);
late = mean(chirpFFTnorm(:,50:60),2);
di = (late-early)./(early+late);
%         di = (early-late);

tmp = find(di<0);
Chrp_index = di;

ALL_CHRP = intersect(R,chrp_tested);
ASSIGN = zeros(length(ALL_CHRP),1);
ASSIGN(chrp(IDboth(onsu2))) = 1;
ASSIGN(chrp(IDboth(ontr2))) = 2;
ASSIGN(chrp(IDboth(onof2))) = 3;
ASSIGN(chrp(IDboth(oftr2))) = 4;
ASSIGN(chrp(IDboth(ofsu2))) = 5;

simil = [matches(IDboth(onsu2),6);matches(IDboth(ontr2),6);matches(IDboth(ofsu2),6);matches(IDboth(oftr2),6)];

    chrpThr = mean(simil)+std(simil);
IDnoflash =  find(S==0);
matchesO = matches(IDnoflash,:);
tmp = find(matchesO(:,6)<chrpThr);
IDnoflash2 = IDnoflash(tmp);
matchesO2 = matchesO(tmp,:);
IDnoflashNoTop = IDnoflash(find(matchesO(:,6)>=chrpThr));
IDnoflashNoTop = setdiff(IDnoflashNoTop,78);

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
[m1,m2] = min(matchesO2');

% figure
% subplot(6,3,2)
% data = comb_spike(IDboth(ontr2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(2,:),'-','color',[0 0 1 alph])
% hold on
% plot(mean(data2,1),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
%
% subplot(6,3,1)
% data = flash_corr(IDboth(ontr2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(mean(data2(:,200:5900),1),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
% title(['putative ON parasol n = ',int2str(length(ontr2)),' / ',num2str(round(100/tot*length(ontr2)*10)/10),'%'])
%
%
% subplot(6,3,3)
% idnow = find(m2==2);
% ASSIGN(chrp(IDnoflash(idnow)))=-2;
% data = comb_spike(IDnoflash(idnow),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(2,:),'-','color',[0 0 1 alph])
% hold on
% plot(mean(data2,1),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
% % title(['putative ON midget (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
% title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
%
%
% subplot(6,3,5)
% data = comb_spike(IDboth(oftr2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(4,:),'-','color',[0 1 0.3 alph])
% hold on
% plot(mean(data2),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
%
% subplot(6,3,4)
% data = flash_corr(IDboth(oftr2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(mean(data2(:,200:5900)),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
% title(['putative OFF parasol n = ',int2str(length(oftr2)),' / ',num2str(round(100/tot*length(oftr2)*10)/10),'%'])
%
%
% subplot(6,3,6)
% idnow = find(m2==4);
% ASSIGN(chrp(IDnoflash(idnow)))=-4;
% data = comb_spike(IDnoflash(idnow),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(4,:),'-','color',[0 1 0.3 alph])
% hold on
% plot(mean(data2),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
% % title(['putative ON parasol (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
% title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
%
%
% subplot(6,3,8)
% data = comb_spike(IDboth(onsu2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(1,:),'-','color',[1 0.2 0.2 alph])
% hold on
% plot(mean(data2,1),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
%
% subplot(6,3,7)
% data = flash_corr(IDboth(onsu2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(mean(data2(:,200:5900)),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
% title(['putative ON midget n = ',int2str(length(onsu2)),' / ',num2str(round(100/tot*length(onsu2)*10)/10),'%'])
%
%
% subplot(6,3,9)
% idnow = find(m2==1);
% ASSIGN(chrp(IDnoflash(idnow)))=-1;
% data = comb_spike(IDnoflash(idnow),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(1,:),'-','color',[1 0.2 0.2 alph])
% hold on
% plot(mean(data2),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
% % title(['putative ON parasol (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
% title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
%
%
% subplot(6,3,11)
% data = comb_spike(IDboth(ofsu2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(5,:),'-','color',[1 0 1 alph])
% hold on
% plot(mean(data2,1),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
%
% subplot(6,3,10)
% data = flash_corr(IDboth(ofsu2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(mean(data2(:,200:5900),1),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
% title(['putative OFF midget n = ',int2str(length(ofsu2)),' / ',num2str(round(100/tot*length(ofsu2)*10)/10),'%'])
%
%
% subplot(6,3,12)
% idnow = find(m2==5);
% ASSIGN(chrp(IDnoflash(idnow)))=-5;
% data = comb_spike(IDnoflash(idnow),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(5,:),'-','color',[1 0 1 alph])
% hold on
% plot(mean(data2),'-k')
% axis tight
% box off
% set(gca,'xtick',[])
% % title(['putative OFF parasol (no flash) n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
% title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
%
%
%
% subplot(6,3,14)
% data = comb_spike(IDboth(onof2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(4,:),'-','color',[1 0.4 0 alph])
% hold on
% plot(mean(data2),'-k')
% axis tight
% box off
% set(gca,'xtick',0:2000:16000,'xticklabel',0:2:16)
% xlabel('sec')
%
% subplot(6,3,13)
% data = flash_corr(IDboth(onof2),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(mean(data2(:,200:5900)),'-k')
% set(gca,'xtick',0:1000:5000,'xticklabel',0:1:5)
% xlabel('sec')
% axis tight
% box off
% title(['putative ON-OFF bistratified n = ',int2str(length(onof2)),' / ',num2str(round(100/tot*length(onof2)*10)/10),'%'])
%
%
% subplot(6,3,15)
% idnow = find(m2==4);
% ASSIGN(chrp(IDnoflash(idnow)))=-4;
% data = comb_spike(IDnoflash(idnow),:);
% data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
% plot(templ(4,:),'-','color',[1 0.4 0 alph])
% hold on
% plot(mean(data2,1),'-k')
% axis tight
% set(gca,'xtick',[])
% title(['n = ',int2str(length(idnow)),' / ',num2str(round(100/tot*length(idnow)*10)/10),'%'])
% box off
%
%
%
% subplot(6,3,16)
% plot([1 1984-200],[128 128],'-','color',[0.4 0.4 0.4],'linewidth',2)
% hold on
% plot([1984-200 1984-200],[128 0],'-','color',[0.4 0.4 0.4],'linewidth',2)
% plot([1984-200 3968-200],[0 0],'-','color',[0.4 0.4 0.4],'linewidth',2)
% plot([3968-200 3968-200],[0 128],'-','color',[0.4 0.4 0.4],'linewidth',2)
% plot([3968-200 5952-200],[128 128],'-','color',[0.4 0.4 0.4],'linewidth',2)
% axis([ 1 5952-200 0 255])
% box off
% set(gca,'xtick',[])
% ylabel('RGB')
%
% subplot(6,3,17)
% plot([prot2(3:481,4); prot2(485:970,4)],'-','color',[0.4 0.4 0.4])
% axis tight
% box off
% axis([ 1 965 0 255])
% set(gca,'xtick',[])
% ylabel('RGB')
%
% subplot(6,3,18)
% plot([prot2(3:481,4); prot2(485:970,4)],'-','color',[0.4 0.4 0.4])
% axis tight
% box off
% axis([ 1 965 0 255])
% set(gca,'xtick',[])
% ylabel('RGB')
%% how about cells that do not respond to chirp
load(fullfile(path2,'EXAMPLE_INFO'),'SUPER_COLL','CLUSTER_INFO')

IDnoChirp = setdiff(chrp_tested,chrp);
R2 = setdiff(R,158:200);
IDnoChirp = intersect(R2,IDnoChirp);

togo = togo([1:47 53:end]); %remove the 5 cells from retina 4

% --- 125 cells respond NOT to chirp but to DG ---
normP = []; ids = []; ids2 = [];
for n = 1:length(IDnoChirp)
    tmp = find(togo==IDnoChirp(n));
    if ~isempty(tmp)
        ids2 = [ids2;n];
    end
    ids = [ids;tmp];
    now = normPeaks(tmp,:,2);
    normP = [normP;now(:)' ];
end
keepDGids = ids2;

if recluster == 1
    CLUST = [];
    EVAL = [];
    for k = 2:30
        [idx,C] = kmeans(normP,k,'Replicates',1000);
        eva1 = evalclusters(normP,idx,'CalinskiHarabasz');
        eva2 = evalclusters(normP,idx,'DaviesBouldin');
        CLUST = [CLUST idx];
        evals = [eva1.CriterionValues; eva2.CriterionValues];
        EVAL = [EVAL evals];
    end
    save(fullfile(path2,'cluster_DG_noMPB'),'CLUST','EVAL','normP','ids')
else
    load(fullfile(path2,'cluster_DG_noMPB'))
end

CC = 8;
CL = CLUST(:,CC-1);
figure
HEATMAP = []; NoC = [];
for c = 1:CC
    tmp  =find(CL==c);
    disp(length(tmp))
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

load(fullfile(path2,'abs_ID_exampleCells'))
load(fullfile(path2,'EXAMPLE_INFO'))


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
HTMP = [];
yh1 = 0.022; yh2 = 0.024; yh3 = 0.018;
yst1 = 0.97; yst2 = 0.92; yst3 = 0.87;
x1 = 0.15; x2 = 0.3; x3 = 0.5;
w0 = 0.1; w1 = 0.15; w2 = 0.06;
xnow = [0.6 0.7 0.8 0.9 0.92];
w3 = 0.06; yh4 = 0.033;
% yh1 = 0.022; yh2 = 0.024; yh3 = 0.018;
% yst1 = 0.97; yst2 = 0.94; yst3 = 0.91;
% x1 = 0.16; x2 = 0.5; x3 = 0.84; x4 = 0.9;
% w1 = 0.28; w2 = 0.06;
% col = [0 0 1 0.35;0 0.2 0.8 0.35;1 0.4 0 0.35;1 0.2 0 0.35;1 0 1 0.35;1 0 1 0.35;];
% col = [0 0.2 0.8 0.35;0 0.8 0.2 0.35;1 0 0 0.35; 1 0 1 0.35; 1 0.2 0 0.35];
col = [1 0 0 0.35; 0 0.2 0.8 0.35; 0.9 0.5 0 0.35;0 0.8 0.2 0.35;1 0 1 0.35];
cols2 = [27,158,119;217,95,2;117,112,179;231,41,138;102,166,30];
% cols2 = [166,206,227;31,120,180;178,223,138;51,160,44;251,154,153];
cols2 = cols2./255;
ex_names = {'A_{1}','A_{2}','A_{3}','A_{4}','A_{5}',...
    'B_{1}','B_{2}','B_{3}','B_{4}','B_{5}','B_{6}','B_{7}','C_{1}','C_{2}','C_{3}'};


figure('position',[310.6000 924.2000 760.0000 1028],'color','w')
% figure('position',[488 70.6000 734.6000 691.4000],'color','w')

mxv = []; mxv2 =[];

axes('position',[x1 yst1 w0 0.01])%[0.08 0.9 0.8 0.02])
rectangle('position',[0 0 1 1],'facecolor',[0.6 0.6 0.6],'edgecolor','none')
hold on
rectangle('position',[1 0 1 1],'facecolor',[0 0 0 ],'edgecolor','none')
rectangle('position',[2 0 1 1],'facecolor',[0.6 0.6 0.6],'edgecolor','none')
plot([1 1.5],[-0.4 -0.4],'-k','linewidth',2)
text(1,-1.4,'1 s','fontname','arial','fontsize',6,'horizontalalignment','left')
axis off
axis([0 3 -inf inf])
box off

axes('position',[x2 yst1 w1 0.01])%[0.08 0.9 0.8 0.02])
plot(prot2(2:end,4),'-k')
hold on
plot([1 2000/16500*971],[-40 -40],'-k','linewidth',2)
text(0,-300,'2 s','fontname','arial','fontsize',6,'horizontalalignment','left')
axis([1 length(prot2)-1 -220 255])
axis off
box off
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
    
    axes('position',[x1 yst1-yh2*(2*(ax))+yh1 w0 yh3])%[0.08 0.9 0.8 0.02])
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
    
      if ax == 1
        plot([100 100],[20 40],'-k','linewidth',1)
        %          plot([150 1150],[20 20],'-k','linewidth',1)
        text(-400,-0,'20 sp/s','fontname','Arial','fontsize',6,'rotation',90)
        %         text(170,30,'0.5 s','fontname','Arial','fontsize',6)
    end
    
    axes('position',[x1 yst1-yh2*(2*(ax)) w0 yh3])%[0.08 0.9 0.8 0.02])
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
    axis([-inf inf 0 1])
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
    %     text(0.05,0.4,['n = ',int2str(length(IDX)),...
    %         ', ',int2str(length(idnow))],'fontname','Arial','fontsize',7)
    text(0.05,0.35,['n = ',int2str(length(idnow)+length(IDX)),' / ',num2str(round(100/tot*(length(idnow)+length(IDX))*10)/10),'%'],...
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
    HTMP = [HTMP;htmp(:)'];
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    hold on
    rectangle('position',[0.5 0.5  6 4],'edgecolor','k','facecolor','none')
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    
    if ax == 1
        axes('position',[0.015 yst1-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'A_1','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[x2-0.035 yst1-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'A_2','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[x3-0.02 yst1-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'A_3','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[xnow(1)-0.03 yst1-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'A_4','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
    end
    
end

% --- cells with flash+chrp but not assigned to the top 5
for ax = 6:7
    IDX = IDflashNotAs(find(flashNotAs==ax));
    ASSIGN(chrp(IDX)) = ax;
    
    axes('position',[x1 yst2-yh2*(2*(ax))+yh1 w0  yh3])%[0.08 0.9 0.8 0.02])
    data = flash_corr(IDX,:);
    %     data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    data2 = data;
    plot(mean(data2(:,200:5900),1),'-k')
    axis([1 5700 0 40])
    box off
    axis off
    set(gca,'xtick',[])
    
    
    axes('position',[x1 yst2-yh2*(2*(ax)) w0 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2(:,200:5900))
    colormap gray
    hold on
    plot([10 10],[0.5 length(IDX)+0.5],'-','color','k','linewidth',5)
    axis off
    
    axes('position',[x2 yst2-yh2*(2*(ax))+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    data = comb_spike(IDX,:);
    data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    plot(mean(data2,1),'-k')
    axis([-inf inf 0 1])
    axis off
    box off
    set(gca,'xtick',[])
    
    axes('position',[x2 yst2-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2)
    colormap gray
    hold on
    plot([10 10],[0.5 length(IDX)+0.5],'-','color','k','linewidth',5)
    axis off
    
    
    
    axes('position',[0.02 yst2-yh2*(2*(ax))-0.01 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
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
    HTMP = [HTMP;htmp(:)'];
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    hold on
    rectangle('position',[0.5 0.5  6 4],'edgecolor','k','facecolor','none')
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    
    if ax == 6
        axes('position',[0.015 yst2-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'B_1','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[x2-0.035 yst2-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'B_2','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[x3-0.02 yst2-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'B_3','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[xnow(1)-0.02 yst2-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'B_4','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
    end
end

% --- cells without flash, but chrp that was not assigned to the top 5
for ax = 8
    IDX = chrp(IDnoflashNoTop);
    ASSIGN(IDX) = ax;
    
    axes('position',[x1 yst2-yh2*(2*(ax))+yh1 w0  yh3])%[0.08 0.9 0.8 0.02])
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
    
    
    axes('position',[x1 yst2-yh2*(2*(ax)) w0 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data(:,200:5900))
    colormap gray
    axis off
    
    axes('position',[x2 yst2-yh2*(2*(ax))+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    data = comb_spike(IDnoflashNoTop,:);
    data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    
    plot(mean(data2,1),'-k')
    
    axis([-inf inf 0 1])
    axis off
    box off
    set(gca,'xtick',[])
    
    axes('position',[x2 yst2-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2)
    colormap gray
    axis off
    
    
    
    axes('position',[0.02 yst2-yh2*(2*(ax))-0.01 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
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
    HTMP = [HTMP;htmp(:)'];
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    hold on
    rectangle('position',[0.5 0.5  6 4],'edgecolor','k','facecolor','none')
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
end

% --- cells without chrp response, but tested
tmp = find(NoC>=5);
cl = 1:CC;
cl = cl(tmp); cnt = 0;
for ax = 9:9+length(cl)-1
    cnt = cnt+1;
    idx = IDnoChirp(keepDGids(find(CLUST(:,CC-1)==cl(cnt))));
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
    
    
    axes('position',[x1 yst3-yh2*(2*(ax))+yh1 w0  yh3])%[0.08 0.9 0.8 0.02])
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
    
    
    axes('position',[x1 yst3-yh2*(2*(ax)) w0 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data1b(:,200:5900))
    colormap gray
    hold on
    plot([10 10],[0.5 size(data,1)+0.5],'-','color','k','linewidth',5)
    axis off
    
    axes('position',[x2 yst3-yh2*(2*(ax))+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    data = comb_all(IDX,:);
    data2 = (data-repmat(min(data')',1,size(data,2)))./(repmat(max(data')',1,size(data,2))-repmat(min(data')',1,size(data,2)));
    plot(mean(data2,1),'-k')
    axis([-inf inf 0 1])
    axis off
    box off
    set(gca,'xtick',[])
    hold on
    %     if ax == 15
    %          plot([150 1150],[0 0],'-k','linewidth',1)
    %         text(170,0.2,'0.5 s','fontname','Arial','fontsize',6)
    %     end
    
    axes('position',[x2 yst3-yh2*(2*(ax)) w1 yh3])%[0.08 0.9 0.8 0.02])
    imagesc(1-data2)
    colormap gray
    axis off
    
    
    
    axes('position',[0.02 yst3-yh2*(2*(ax))-0.01 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
    if ax == 9
        text(0.05,0.9,'cluster 9','fontname','Arial','fontsize',7,'fontweight','bold')
        %         text(0.05,0.65,'speed-tuned?','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 10
        text(0.05,0.9,'cluster 10','fontname','Arial','fontsize',7,'fontweight','bold')
        %         text(0.05,0.65,'temporal tuning','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax ==11
        text(0.05,0.9,'cluster 11','fontname','Arial','fontsize',7,'fontweight','bold')
        %         text(0.05,0.65,'fast-big','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 12
        text(0.05,0.9,'cluster 12','fontname','Arial','fontsize',7,'fontweight','bold')
        %         text(0.05,0.65,'4Hz-big','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 13
        text(0.05,0.9,'cluster 13','fontname','Arial','fontsize',7,'fontweight','bold')
        %         text(0.05,0.65,'slow','fontname','Arial','fontsize',7,'fontweight','bold')
        
    elseif ax == 14
        text(0.05,0.9,'cluster 14','fontname','Arial','fontsize',7,'fontweight','bold')
        %         text(0.05,0.65,'medium speed/size','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 15
        text(0.05,0.9,'cluster 15','fontname','Arial','fontsize',7,'fontweight','bold')
    elseif ax == 16
        text(0.05,0.9,'cluster 16','fontname','Arial','fontsize',7,'fontweight','bold')
        
    end
    text(0.05,0.65,['n = ',int2str(length(IDX)),' / ',num2str(round(100/tot*length(IDX)*10)/10),'%'],...
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
    HTMP = [HTMP;htmp(:)'];
    imagesc(htmp)
    colormap gray
    axis tight
    hold on
    rectangle('position',[0.5 0.5  6 4],'edgecolor','k','facecolor','none')
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    if ax == 16
        set(gca,'xtick',1:6,'xticklabel',[4000 2000 1000 500 200 100],'fontname','Arial','fontsize',7)
        xtickangle(90)
        set(gca,'ytick',1:4,'yticklabel',[8 4 2 1],'fontname','Arial','fontsize',7)
        text(-2,2.5,'Hz','fontname','Arial','fontsize',7,'horizontalalignment','left')
        text(3,8,'\mum','fontname','Arial','fontsize',7,'horizontalalignment','left')
    else
        axis off
        
    end
    
    if ax == 9
        axes('position',[0.015 yst3-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'C_1','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[x2-0.035 yst3-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'C_2','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[x3-0.02 yst3-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'C_3','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
        axes('position',[xnow(1)-0.02 yst3-yh2*(2*(ax))+yh1 0.015 yh3*2])%[0.08 0.9 0.8 0.02])
        text(0.05,0.9,'C_4','fontname','Arial','fontsize',10,'fontweight','bold')
        axis off
        box off
    end
end

% --- add example cells from fig 2 ----
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
        axes('position',[x3-0.045 yst1-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
    elseif c<9
        axes('position',[x3-0.045 yst2-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x3-0.045 yst3-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
    end
    
    for t = 1:length(tmp)
        if t<4
            text(0.05,1.25-0.4*t,ex_names{tmp(t)},'fontname','Arial','fontsize',7,'color',cols2(t,:),'fontweight','bold')
        elseif t == 4
            text(0.99,1.25-0.4,ex_names{tmp(t)},'fontname','Arial','fontsize',7,'color',cols2(t,:),'fontweight','bold')
        else
            text(0.99,1.25-0.8,ex_names{tmp(t)},'fontname','Arial','fontsize',7,'color',cols2(t,:),'fontweight','bold')
        end
    end
    axis off
    box off
end
%% cluster remaining cells with DG but not tested for chrp
ASSIGN2 = ASSIGN;
for a = -5:-1
    tmp = find(ASSIGN2 == a);
    ASSIGN2(tmp) = -a;
end
CC = unique(ASSIGN2);
tmp = find(CC>0); CC = CC(tmp);
done = find(ASSIGN2~=0);
missing = setdiff(togo,done);
missing = setdiff(missing,chrp_tested);

CLUH = cell(1,1); ALL =[];
for c = 1:length(CC)
    idxc = find(ASSIGN2 == CC(c));
    ht = [];
    for n = 1:length(idxc)
        tmp = find(togo == idxc(n));
        if ~isempty(tmp)
            now = normPeaks(tmp,:,2);
            now = 1-now;
            now = flipud(fliplr(now));
            D = HTMP - repmat(now,16,1);
            S = sum(abs(D)');
            [m1,m2]=min(abs(S));
            ht = [ht;m1 m2];
        end
        
    end
    CLUH{c,1} = ht;
    ALL =[ALL;ht];
end
[s1,s2] = sort(ALL(:,1));
thrsh = s1(length(s1)*.95);

missingHTMP = []; SS = [];
for m = 1:length(missing)
    tmp = find(togo == missing(m));
    now = normPeaks(tmp,:,2);
    now = 1-now;
    now = flipud(fliplr(now));
    %    [cX,cY,oX,oY,dX,dY] = get_GaussHeatmaps(now);
    D = HTMP - repmat(now,16,1);
    S = sum(abs(D)');
    [m1,m2]=min(abs(S));
    SS = [SS;m1 m2];
    missingHTMP = [missingHTMP;now];
end

ASSIGN3 = ASSIGN2;
tmp = find(SS(:,1)<=thrsh);
ASSIGN3(missing(tmp)) = SS(tmp,2);

for ax = 1:16
    if ax == 1
        tm = 2;
    elseif ax == 2
        tm = 4;
    elseif ax ==3
        tm = 1;
    elseif ax == 4
        tm = 5;
    elseif ax ==5
        tm = 3;
    else
        tm = ax;
    end
    IDX = find(ASSIGN3==tm);
    if ax<6
        axes('position',[0.02 yst1-yh2*(2*(ax))-0.005 0.035 0.01])%[0.08 0.9 0.8 0.02])
    elseif ax<8
        axes('position',[0.02 yst2-yh2*(2*(ax))-0.008 0.035 0.01])%[0.08 0.9 0.8 0.02])
    elseif ax == 8
        axes('position',[0.02 yst2-yh2*(2*(ax))-0.008 0.035 0.01])%[0.08 0.9 0.8 0.02])
    elseif ax >8
        axes('position',[0.02 yst3-yh2*(2*(ax))+0.002 0.035 0.01])%[0.08 0.9 0.8 0.02])
    end
    % axes('position',[0.02 yst3-yh2*(2*(ax))-0.01 0.035 yh1*2])%[0.08 0.9 0.8 0.02])
    text(0.05,0.65,['n_all = ',int2str(length(IDX)),' / ',num2str(round(100/(tot+length(tmp))*length(IDX)*10)/10),'%'],...
        'fontname','Arial','fontsize',7,'interpreter','none')
    axis off
    box off
end
%% add distributions to plot




IDSdg = intersect(togo,find(ASSIGN2>0));
ids = [];
for i = 1:length(IDSdg)
    tmp = find(togo==IDSdg(i));
    ids = [ids;tmp];
end
normP = normPeaks(ids,:,2);
[cX,cY,oX,oY,dX,dY] = get_GaussHeatmaps(normPeaks(ids,:,:));

rangeF = length(chirpAmp)-chirpAmp(1);
steps = round(1:length(chirpAmp)/3:length(chirpAmp));
chk1=steps(1):steps(2)-1;
chk3=steps(3):length(chirpAmp);
cntrstIDX = [];
for ce = 1:size(chirpAmp,1)
    low = sum(abs(chirpAmp(ce,chk1)));
    high = sum(abs(chirpAmp(ce,chk3)));
    
    idx = (high-low)./(low+high);
    cntrstIDX = [cntrstIDX;idx];
end


spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
SP = unique(spat); FR =unique(freq);

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
    cc = CC(tm);
    tmp = find(ASSIGN2(IDSdg)==cc);
    dat = cX(tmp);
    %---spatial ---
    axes('position',[xnow(3) yst1-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat2 = [];
    for d = 1:length(dat)
        now = dat(d);
        in = floor(now);
        rest = now-in;
        if rest >0
            slp = (SP(in+1)-SP(in));
            um = slp*rest+SP(in);
            dat2 = [dat2;um];
        else
            dat2 = [dat2;SP(in)];
        end
    end
    dat2 = 4000-dat2;
    plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 0.5]},'xyOri','flipped')
    hold on
    plot([mean(dat2) mean(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    exnt = 0;
    for exc = 1:length(COLL_ID)
        if ~isempty(find(IDSdg(tmp)==COLL_ID(exc)))
            exnt = exnt+1;
            tp = find(IDSdg(tmp)==COLL_ID(exc));
            scatter(dat2(tp),1,15,cols2(exnt,:),'linewidth',1)
        end
    end
    
    axis([0 4000 0.3 1.7 ])
    set(gca,'xtick',[0 2000 4000],'xticklabel',[4 2 0])
    %     xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    if ax == 1
        title({'spatial', 'pref.'},'fontname','Arial','fontsize',7)
    end
    
    if ax == 5
        xlabel('mm','fontname','Arial','fontsize',7)
    end
    
    %---temporal ---
    axes('position',[xnow(4) yst1-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat = 5-cY(tmp);
    dat2 = [];
    for d = 1:length(dat)
        now = dat(d);
        in = floor(now);
        rest = now-in;
        if rest >0
            slp = (FR(in+1)-FR(in));
            um = slp*rest+FR(in);
            dat2 = [dat2;um];
        else
            dat2 = [dat2;FR(in)];
        end
    end
    plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(dat2) mean(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    
    axis([0 8 0.3 1.7 ])
    set(gca,'xtick',[0 4 8],'xticklabel',[0 4 8])
    %      xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    exnt = 0;
    for exc = 1:length(COLL_ID)
        if ~isempty(find(IDSdg(tmp)==COLL_ID(exc)))
            exnt = exnt+1;
            tp = find(IDSdg(tmp)==COLL_ID(exc));
            scatter(dat2(tp),1,15,cols2(exnt,:),'linewidth',1)
        end
    end
    if ax == 1
        title({'temporal', 'pref.'},'fontname','Arial','fontsize',7)
    end
    if ax == 5
        xlabel('Hz','fontname','Arial','fontsize',7)
    end
    
    
    %---broadness ---
    %     axes('position',[xnow(3) yst1-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    %     dat2 = (dX(tmp)+dY(tmp))./2;
    %     plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
    %         {'.'},'distributionColors',{[0 0 0 1]},'xyOri','flipped')
    %     hold on
    %     plot([median(dat2) median(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    %
    %     axis([0 4 0.3 1.7 ])
    %     set(gca,'xtick',[0 2 4],'xticklabel',[0 2 4])
    % %      xtickangle(90)
    %     set(gca,'ytick',[],'ycolor','w')
    %     if ax == 1
    %        title({'width of', 'S-T pref.'},'fontname','Arial','fontsize',7)
    %     end
    
    % ----chrp---
    axes('position',[xnow(1) yst1-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat = [];
    for t = 1:length(tmp)
        id = find(chrp==IDSdg(tmp(t)));
        dat = [dat;id];
    end
    plotSpread(Chrp_index(dat),'distributionIdx',repmat(1,length(dat),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(Chrp_index(dat)) mean(Chrp_index(dat))],[0.5 1.5],'-','color','k','linewidth',1)
    plot([0 0],[0.3 1.7],'--','color',[0.4 0.4 0.4])
    axis([-1 1 0.3 1.7 ])
    set(gca,'xtick',[-1 0 1],'xticklabel',[-1 0 1])
    %      xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    exnt = 0;
    for exc = 1:length(COLL_ID)
        if ~isempty(find(IDSdg(tmp)==COLL_ID(exc)))
            exnt = exnt+1;
            tp = find(IDSdg(tmp)==COLL_ID(exc));
            scatter(Chrp_index(dat(tp)),1,15,cols2(exnt,:),'linewidth',1)
        end
    end
    if ax == 1
        title({'chirp freq.', 'index'},'fontname','Arial','fontsize',7)
    end
    if ax == 5
        xlabel('a.u.','fontname','Arial','fontsize',7)
    end
    
    % --- chrp contrast ---
    axes('position',[xnow(2) yst1-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat = [];
    for t = 1:length(tmp)
        id = find(chrp==IDSdg(tmp(t)));
        dat = [dat;id];
    end
    plotSpread(cntrstIDX(dat),'distributionIdx',repmat(1,length(dat),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(cntrstIDX(dat)) mean(cntrstIDX(dat))],[0.5 1.5],'-','color','k','linewidth',1)
    plot([0 0],[0.3 1.7],'--','color',[0.4 0.4 0.4])
    axis([-1 1 0.3 1.7 ])
    set(gca,'xtick',[-1 0 1],'xticklabel',[-1 0 1])
    %      xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    exnt = 0;
    for exc = 1:length(COLL_ID)
        if ~isempty(find(IDSdg(tmp)==COLL_ID(exc)))
            exnt = exnt+1;
            tp = find(IDSdg(tmp)==COLL_ID(exc));
            scatter(cntrstIDX(dat(tp)),1,15,cols2(exnt,:),'linewidth',1)
        end
    end
    if ax == 1
        title({'contrast', 'pref. index'},'fontname','Arial','fontsize',7)
    end
    if ax == 5
        xlabel('a.u.','fontname','Arial','fontsize',7)
    end
    
end

for ax = 6:8
    cc = CC(ax);
    tmp = find(ASSIGN2(IDSdg)==cc);
    dat = cX(tmp);
    mx = median(cX(tmp)); my = median(cY(tmp));
    %---spatial ---
    axes('position',[xnow(3) yst2-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat2 = [];
    for d = 1:length(dat)
        now = dat(d);
        in = floor(now);
        rest = now-in;
        if rest >0
            slp = (SP(in+1)-SP(in));
            um = slp*rest+SP(in);
            dat2 = [dat2;um];
        else
            dat2 = [dat2;SP(in)];
        end
    end
    dat2 = 4000-dat2;
    
    plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(dat2) mean(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    exnt = 0;
    for exc = 1:length(COLL_ID)
        if ~isempty(find(IDSdg(tmp)==COLL_ID(exc)))
            exnt = exnt+1;
            tp = find(IDSdg(tmp)==COLL_ID(exc));
            scatter(dat2(tp),1,15,cols2(exnt,:),'linewidth',1)
        end
    end
    axis([0 4000 0.3 1.7 ])
    set(gca,'xtick',[0 2000 4000],'xticklabel',[4 2 0])
    %     xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    
    if ax == 8
        xlabel('mm','fontname','Arial','fontsize',7)
    end
    
    %---temporal ---
    axes('position',[xnow(4) yst2-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat = 5-cY(tmp);
    dat2 = [];
    for d = 1:length(dat)
        now = dat(d);
        in = floor(now);
        rest = now-in;
        if rest >0
            slp = (FR(in+1)-FR(in));
            um = slp*rest+FR(in);
            dat2 = [dat2;um];
        else
            dat2 = [dat2;FR(in)];
        end
    end
    plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(dat2) mean(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    exnt = 0;
    for exc = 1:length(COLL_ID)
        if ~isempty(find(IDSdg(tmp)==COLL_ID(exc)))
            exnt = exnt+1;
            tp = find(IDSdg(tmp)==COLL_ID(exc));
            scatter(dat2(tp),1,15,cols2(exnt,:),'linewidth',1)
        end
    end
    axis([0 8 0.3 1.7 ])
    set(gca,'xtick',[0 4 8],'xticklabel',[0 4 8])
    %  xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    if ax == 8
        xlabel('Hz','fontname','Arial','fontsize',7)
    end
    
    %---broadness ---
    %     axes('position',[xnow(3) yst2-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    %     dat2 = (dX(tmp)+dY(tmp))./2;
    %     plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
    %         {'.'},'distributionColors',{[0 0 0 1]},'xyOri','flipped')
    %     hold on
    %     plot([median(dat2) median(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    %
    %     axis([0 4 0.3 1.7 ])
    %     set(gca,'xtick',[0 2 4],'xticklabel',[0 2 4])
    % %    xtickangle(90)
    %     set(gca,'ytick',[],'ycolor','w')
    
    % ----chrp---
    axes('position',[xnow(1) yst2-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat = [];
    for t = 1:length(tmp)
        id = find(chrp==IDSdg(tmp(t)));
        dat = [dat;id];
    end
    plotSpread(Chrp_index(dat),'distributionIdx',repmat(1,length(dat),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(Chrp_index(dat)) mean(Chrp_index(dat))],[0.5 1.5],'-','color','k','linewidth',1)
    plot([0 0],[0.3 1.7],'--','color',[0.4 0.4 0.4])
    axis([-1 1 0.3 1.7 ])
    set(gca,'xtick',[-1 0 1],'xticklabel',[-1 0 1])
    %    xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    exnt = 0;
    for exc = 1:length(COLL_ID)
        if ~isempty(find(IDSdg(tmp)==COLL_ID(exc)))
            exnt = exnt+1;
            tp = find(IDSdg(tmp)==COLL_ID(exc));
            scatter(Chrp_index(dat(tp)),1,15,cols2(exnt,:),'linewidth',1)
        end
    end
    if ax == 8
        xlabel('a.u.','fontname','Arial','fontsize',7)
    end
    
    % --- chrp contrast ---
    axes('position',[xnow(2) yst2-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat = [];
    for t = 1:length(tmp)
        id = find(chrp==IDSdg(tmp(t)));
        dat = [dat;id];
    end
    plotSpread(cntrstIDX(dat),'distributionIdx',repmat(1,length(dat),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(cntrstIDX(dat)) mean(cntrstIDX(dat))],[0.5 1.5],'-','color','k','linewidth',1)
    plot([0 0],[0.3 1.7],'--','color',[0.4 0.4 0.4])
    axis([-1 1 0.3 1.7 ])
    set(gca,'xtick',[-1 0 1],'xticklabel',[-1 0 1])
    exnt = 0;
    for exc = 1:length(COLL_ID)
        if ~isempty(find(IDSdg(tmp)==COLL_ID(exc)))
            exnt = exnt+1;
            tp = find(IDSdg(tmp)==COLL_ID(exc));
            scatter(cntrstIDX(dat(tp)),1,15,cols2(exnt,:),'linewidth',1)
        end
    end
    %      xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    if ax == 8
        xlabel('a.u.','fontname','Arial','fontsize',7)
    end
end

for ax = 9:16
    cc = CC(ax);
    tmp = find(ASSIGN2(IDSdg)==cc);
    dat = cX(tmp);
    mx = median(cX(tmp)); my = median(cY(tmp));
    %---spatial ---
    axes('position',[xnow(3) yst3-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat2 = [];
    for d = 1:length(dat)
        now = dat(d);
        in = floor(now);
        rest = now-in;
        if rest >0
            slp = (SP(in+1)-SP(in));
            um = slp*rest+SP(in);
            dat2 = [dat2;um];
        else
            dat2 = [dat2;SP(in)];
        end
    end
    dat2 = 4000-dat2;
    
    plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(dat2) mean(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    
    axis([0 4000 0.3 1.7 ])
    set(gca,'xtick',[0 2000 4000],'xticklabel',[4 2 0])
    %     xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    if ax == 16
        xlabel('mm','fontname','Arial','fontsize',7)
    end
    
    %---temporal ---
    axes('position',[xnow(4) yst3-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    dat = 5-cY(tmp);
    dat2 = [];
    for d = 1:length(dat)
        now = dat(d);
        in = floor(now);
        rest = now-in;
        if rest >0
            slp = (FR(in+1)-FR(in));
            um = slp*rest+FR(in);
            dat2 = [dat2;um];
        else
            dat2 = [dat2;FR(in)];
        end
    end
    plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
        {'.'},'distributionColors',{[0.7 0.7 0.7 1]},'xyOri','flipped')
    hold on
    plot([mean(dat2) mean(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    
    axis([0 8 0.3 1.7 ])
    set(gca,'xtick',[0 4 8],'xticklabel',[0 4 8])
    %     xtickangle(90)
    set(gca,'ytick',[],'ycolor','w')
    if ax == 16
        xlabel('Hz','fontname','Arial','fontsize',7)
    end
    
    %---broadness ---
    %     axes('position',[xnow(3) yst3-yh2*(2*(ax)) w3 yh4])%[0.08 0.9 0.8 0.02])
    %     dat2 = (dX(tmp)+dY(tmp))./2;
    %     plotSpread(dat2,'distributionIdx',repmat(1,length(dat2),1),'distributionMarkers',...
    %         {'.'},'distributionColors',{[0 0 0 1]},'xyOri','flipped')
    %     hold on
    %     plot([median(dat2) median(dat2)],[0.5 1.5],'-','color','k','linewidth',1)
    %
    %     axis([0 4 0.3 1.7 ])
    %     set(gca,'xtick',[0 2 4],'xticklabel',[0 2 4])
    % %     xtickangle(90)
    %     set(gca,'ytick',[],'ycolor','w')
    %     if ax == 15
    %        xlabel('a.u.','fontname','Arial','fontsize',7)
    %     end
    
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


exportgraphics(gcf,fullfile(pathSave,'NEW_Fig5_clustering.png'),'Resolution',600)
exportgraphics(gcf,fullfile(pathSave,'NEW_Fig5_clustering.eps'),'Resolution',600)

