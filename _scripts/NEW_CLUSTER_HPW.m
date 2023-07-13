clear
close all
clc

%% pathes
path1 = 'D:\_data\all_species';
% path2 = 'C:\Users\Lenovo2020\OneDrive - imec\HumRet';
path2 = 'D:\_data\Human';
pathcode = 'C:\Users\Lenovo2020\OneDrive - imec\HumRet\_scripts';
pathSave = 'C:\Users\Lenovo2020\OneDrive - imec\HumanRetina\NEWPAPER';
addpath(genpath(pathcode))

recluster = 0;

%% stimulus info
prot=read_header_field_heka(fullfile('C:\Users\Lenovo2020\OneDrive - imec\HumRet\data\20130211c\HEKA'),'20130211_C1#0211_HumanPig_FFFlashBW%dur_2000_gap_2000%_FW0ND4.phys','Stimulus Protocol');
prot2=read_header_field_heka(fullfile('C:\Users\Lenovo2020\OneDrive - imec\HumRet\data\20130211c\HEKA'),'20130211_C1#0283_HumanPig_chirp%mean_128_amplitude_128_pause_500%_FW0ND4.phys','Stimulus Protocol');

%% HUMAN: get spike rates chrp and flash and transiency
load(fullfile(path1,'h_info'))
load(fullfile(path1,'h_list_date_responses'))

chrp_spike = [];
for g=1:length(chrp)
    
    currUnit=info{chrp(g),2};
    currDate=info{chrp(g),1};
    
    startF=info{chrp(g),4};
    stopF=info{chrp(g),5};
    
    %pathes%
    dpath=fullfile(path2,currDate);
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
    dpath=fullfile(path2,currDate);
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
    dpath=fullfile(path2,currDate);
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
    dpath=fullfile(path2,currDate);
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
%% HUMAN: collect data for responsive units
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

old = load(fullfile(path1,'Cowan2020_old'));
load(fullfile(path1,'Cowan2020'))
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

%% HUMAN normalize templates and our data for matching
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
%% HUMAN find cells with flash + chirp responses and collect polarity & transience

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
[chirpFFTnorm,chirpAmp,freqSmoothChrp,togoChrp]=get_ChirpFFTnorm_new(path2,info);


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

ofsu2 = ofsu;

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

%% HUMAN assign cells without matching chirp (6-7)
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

%% HUMAN how about cells that do not respond to chirp (cluster 9-16)
% load(fullfile(path2,'EXAMPLE_INFO'),'SUPER_COLL','CLUSTER_INFO')
load(fullfile(path1,'h_normPeaks'))
IDnoChirp = setdiff(chrp_tested,chrp);
R2 = setdiff(R,158:200);
IDnoChirp = intersect(R2,IDnoChirp);

keeptogo = togo;
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
    load(fullfile(path1,'cluster_DG_noMPB'))
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


%% main plot v2
[m1,m2] = min(matchesO2');
tot = 278;
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


figure('position',[795.4000 845 760 977.6000],'color','w')
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
% for c = 1:ax
%     tmp = find(ASSIGN(COLL_ID)==c);
%     if c == 1
%         newcl = 3;
%     elseif c == 2
%         newcl = 1;
%     elseif c == 3
%         newcl = 5;
%     elseif c == 4
%         newcl = 2;
%     elseif c == 5
%         newcl = 4;
%     else
%         newcl = c;
%     end
%     if c <6
%         axes('position',[x3-0.045 yst1-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
%     elseif c<9
%         axes('position',[x3-0.045 yst2-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
%     else
%         axes('position',[x3-0.045 yst3-yh2*(2*(newcl))+yh3/4 0.018 yh3*1.2])%[0.08 0.9 0.8 0.02])
%     end
%
%     for t = 1:length(tmp)
%         if t<4
%             text(0.05,1.25-0.4*t,ex_names{tmp(t)},'fontname','Arial','fontsize',7,'color',cols2(t,:),'fontweight','bold')
%         elseif t == 4
%             text(0.99,1.25-0.4,ex_names{tmp(t)},'fontname','Arial','fontsize',7,'color',cols2(t,:),'fontweight','bold')
%         else
%             text(0.99,1.25-0.8,ex_names{tmp(t)},'fontname','Arial','fontsize',7,'color',cols2(t,:),'fontweight','bold')
%         end
%     end
%     axis off
%     box off
% end
%% HUMAN cluster remaining cells with DG but not tested for chrp
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
H_ASSIGN = ASSIGN3;
H_CHRP = comb_spike;
H_CHRPsm = ours2;
H_CHRP_ID = chrp;
H_CHRP_ALL = comb_all;
H_CHRO_ALLID = chrp_tested;
H_FLASH_ID = R;
H_HTMP = HTMP;
H_FLASH = flash_all;
H_MATCH = matches;
H_MATCH2 = matchesO;
H_MATCH3 = matchesO2;
H_TOGO = togo;
H_DG = normPeaks([1:47 53:end],:,:);
H_CHRP_THR = chrpThr;
H_ChrpIndex = Chrp_index;
%% now add mice and pigs
% human stats
Hum = []; Hum2 = [];
for c = 1:size(H_CHRP,1)
    now = H_CHRP(c,:);
    [yupper,ylower] = envelope(now,200,'peak');
    now2 = movmedian(now,600);
    now2 = (now2-min(now2))./(max(now2)-min(now2));
    sub = now2(8468:9005);
    Hum = [Hum;now2];
    Hum2 = [Hum2;yupper];
end
IDS = H_ASSIGN(H_CHRP_ID);
templH = []; Hchrpvar = [];
for a = 1:8
    tmp = find(IDS==a);
    templH = [templH;mean(Hum(tmp,:))];
    mtchs = [];
    for b = 1:length(tmp)
        d1 = Hum(tmp(b),:);
        D = templH-d1;
        DS = sum(abs(D)')./size(d1,2);
        mtchs = [mtchs;DS(1) ];
    end
    Hchrpvar = [Hchrpvar;median(mtchs) mean(mtchs) min(mtchs) max(mtchs) quantile(mtchs,0.75)];
end

HDGvar = [];
for a = 1:16
    tmp = find(H_ASSIGN==a);
    tmp2 = intersect(tmp,H_TOGO);
    idsnow = [];
    for t = 1:length(tmp2)
        tmp = find(H_TOGO==tmp2(t));
        idsnow = [idsnow;tmp];
    end
    now = H_DG(idsnow,:,2);
    now = 1-now;
    now = flipud(fliplr(now));
    ht = [];
    for n = 1:size(now,1)
        D = HTMP(a,:) - (now(n,:));
        S = sum(abs(D)');
        ht = [ht;S];
    end
    HDGvar = [HDGvar;median(ht) mean(ht) min(ht) max(ht) quantile(ht,0.95)];
end

cCOLL = [];
cCOLL2 = [];
cCOLL3 = [];
ccomb_spike = [];
cflash_all = [];
ccomb_all = [];
cchrp_tested = [];
cflash_spike = [];
cflash_trans = [];
cCOLL2b = [];
cchirp_index = [];
cchrp = [];
cflash = [];
cR = [];
cDGmatch = [];
cmatches = [];
call = [];
ctogo = [];
cnormPeaks = [];
SPID_togo = [];
SPID_R = [];
SPID_chrp = [];
SPID_all = [];
for sp = {'w','p'}
    spnow = sp{1};
    if spnow == 'w'
        path2 = 'D:\_data\WT';
        spid = 1;
    else
        path2 = 'D:\_data\Pig';
        spid = 2;
    end
    %% get spike rates chrp and flash and transiency
    load(fullfile(path1,[spnow,'_info']))
    load(fullfile(path1,[spnow,'_list_date_responses']))
    
    chrp_spike = [];
    for g=1:length(chrp)
        
        currUnit=info{chrp(g),2};
        currDate=info{chrp(g),1};
        
        startF=info{chrp(g),4};
        stopF=info{chrp(g),5};
        
        %pathes%
        dpath=fullfile(path2,currDate);
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
        dpath=fullfile(path2,currDate);
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
        dpath=fullfile(path2,currDate);
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
        dpath=fullfile(path2,currDate);
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
    IInt = intersect(chrp,flash);
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
    
    %% normalize data
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
    
    
    
    %% only cluster cells with DG & chrp & flash responses
    load(fullfile(path1,[spnow,'_normPeaks.mat']))
    
    all = intersect(chrp,togo);
    
    % collect flash type
    S = sum(flash_corr');
    IDboth = find(S~=0);
    all = intersect(chrp(IDboth),all);
    flash_type = [];
    COLL = [];
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
        if ~isempty(find(chrp(IDboth(ii))==all))
            COLL = [COLL;flash_type(end,:)];
        end
    end
    
    % chirp matches
    matches = []; COLL2 = [];COLL2b = [];ChirpIndex = [];
    for d = 1:length(chrp)
        d1 = ours(d,:);
        D = templH-repmat(d1,8,1);
        DS = sum(abs(D)')./size(d1,2);
        matches = [matches;DS min(DS)];
        if ~isempty(find(chrp(d)==all))
            COLL2 = [COLL2;matches(end,:)];%best chirp match
        end
        
        early = median(ours2(d,500:2000),2);
        late = median(ours2(d,5000:7000),2);
        di = early-late;
di2 = (late-early)./(early+late);
   ChirpIndex = [ChirpIndex;di2];
        if ~isempty(find(chrp(d)==all))
            COLL2b = [COLL2b;di];%high/low frequency index
        end
        
    end
    
    % DG matches
    COLL3 = []; DGmatch = [];
    for a = 1:length(togo)
        now = normPeaks(a,:,2);
        now = 1-now;
        now = flipud(fliplr(now));
        D = H_HTMP - repmat(now,16,1);
        S = sum(abs(D)');
        
        [m1,m2]=min(abs(S));
        DGmatch = [DGmatch;S];
        if ~isempty(find(togo(a)==all))
            COLL3 = [COLL3;S];
        end
    end
    cCOLL = [cCOLL;COLL];
    cCOLL2 = [cCOLL2;COLL2];
    cCOLL2b = [cCOLL2b;COLL2b];
    cchirp_index = [cchirp_index;ChirpIndex];
    cCOLL3 = [cCOLL3;COLL3];
    ccomb_spike = [ccomb_spike;comb_spike];
    ccomb_all = [ccomb_all;comb_all];
    cchrp_tested = [cchrp_tested;chrp_tested];
    cflash_all = [cflash_all;flash_all];
    cflash_spike = [cflash_spike;flash_spike];
    cflash_trans = [cflash_trans;flash_trans];
    cflash = [cflash;flash];
    cchrp = [cchrp;chrp];
    cR = [cR;R];
    cDGmatch = [cDGmatch;DGmatch];
    cmatches = [cmatches;matches];
    call = [call;all];
    ctogo = [ctogo;togo];
    cnormPeaks = cat(1,cnormPeaks,normPeaks);
    SPID_togo = [SPID_togo;repmat(spid,length(togo),1)];
    SPID_R = [SPID_R;repmat(spid,length(R),1)];
    SPID_chrp = [SPID_chrp;repmat(spid,length(chrp),1)];
    SPID_all = [SPID_all;repmat(spid,length(all),1)];
end

%% now find cells that fit in one of the existing clusters 1-8
id1 = find(SPID_togo==1); mxt1 = ctogo(id1(end));
id2 = find(SPID_togo==2); mxt2 = ctogo(id2(end));
ASSIGN = zeros(length(SPID_togo),1);
tester = [];
for c = 1:length(cCOLL3)
    sp = SPID_all(c);
    if sp == 1
        tmp = find(ctogo(id1) == call(c));
        cid = tmp;
    else
        tmp = find(ctogo(id2) == call(c));
        cid = tmp+length(id1);
    end
    tester = [tester;ctogo(cid)];
    di = cCOLL2b(c);
    poss = 1:8;% only these clusters possible since the others have no flash response
    %check flash
    if cCOLL(c,1) == 2 && cCOLL(c,2)<400 %onoff
        if di<0
            poss = 3;
        else
            poss = [];
        end
    elseif cCOLL(c,1) == 1 && cCOLL(c,2) <400 %ON trans
        if di<0
            poss = [2 6 8];
        else
            poss = [6 8];
            % check chirp
            mn = cCOLL2(c,end);
            tmp = find(cCOLL2(c,1:end-1)==mn);
            H = Hchrpvar(poss,:);
            tmp = find(H(:,end)>=mn);
            poss = poss(tmp);
            
        end
    elseif cCOLL(c,1) == 1 && cCOLL(c,2) >=400%ON sus
        if abs(di)<0.2
            poss = [1 6 8];
        else
            poss = [6 8];
            % check chirp
            mn = cCOLL2(c,end);
            tmp = find(cCOLL2(c,1:end-1)==mn);
            H = Hchrpvar(poss,:);
            tmp = find(H(:,end)>=mn);
            poss = poss(tmp);
        end
    elseif cCOLL(c,1) == -1 && cCOLL(c,2) <400%OFF tr
        if di<0
            poss = [4 7];
        else
            poss = 7;
            % check chirp
            mn = cCOLL2(c,end);
            tmp = find(cCOLL2(c,1:end-1)==mn);
            H = Hchrpvar(poss,:);
            tmp = find(H(:,end)>=mn);
            poss = poss(tmp);
        end
    elseif cCOLL(c,1) == -1 && cCOLL(c,2) >=400%OFF sus
        poss = [5 7];
        % check chirp
        mn = cCOLL2(c,end);
        tmp = find(cCOLL2(c,1:end-1)==mn);
        H = Hchrpvar(poss,:);
        tmp = find(H(:,end)>=mn);
        poss = poss(tmp);
    else
        poss = [];
    end
    
    
    
    
    
    
    if ~isempty(poss)
        % check DG
        dg = cCOLL3(c,poss);
        HDGvar(poss,:);
        tmp = HDGvar(poss,end)-dg';
        tmp2 = find(tmp>=0);
        tmp = tmp(tmp2);
        poss = poss(tmp2);
        if ~isempty(poss)
            [m1,m2]=max(tmp);
            ASSIGN(cid) = poss(m2);
        end
    else
        ASSIGN(cid) = -1;
        
    end
end
%% are there cells with chirp but without flash? -> NO, so we don't worry about it
id11 = find(SPID_all==1);
id22 = find(SPID_all==2);

all1 = call(id11); flash1 = cR(SPID_R==1);
setdiff(all1,flash1)
all2 = call(id22); flash2 = cR(SPID_R==2);
setdiff(all2,flash2)
%% assign cells without flash/chirp to cluster 9-16
id111 = find(SPID_chrp==1);
id222 = find(SPID_chrp==2);
DD = diff(cflash); tp = find(DD<0);
id1111 = 1:tp;
id2222 = tp+1:length(cflash);

hdg = HDGvar(9:end,:);
for s = 1:2
    if s == 1
        %         tknow = setdiff(ctogo(id1),call(id11));
        tknow = setdiff(ctogo(id1),cchrp(id111));
        tknow = setdiff(tknow,cflash(id1111));
        tknow2 = [];tknow3 = [];
        for t = 1:length(tknow)
            tmp = find(ctogo(id1)==tknow(t));
            tknow2 = [tknow2;t];
            tknow3 = [tknow3;tmp];
        end
    else
        %         tknow = setdiff(ctogo(id2),call(id22));
        tknow = setdiff(ctogo(id2),cchrp(id222));
        tknow = setdiff(tknow,cflash(id2222));
        tknow2 = [];
        for t = 1:length(tknow)
            tmp = find(ctogo(id2)==tknow(t));
            tknow2 = [tknow2;t+length(id1)];
            tknow3 = [tknow3;tmp+length(id1)];
        end
    end
    
    DG = cDGmatch(tknow2,9:16);
    for d = 1:length(tknow2)
        [m1,m2] = min(DG(d,:));
        if (hdg(m2,end)-m1)>=0
            ASSIGN(tknow3(d)) = m2+8;
        end
    end
end

%% cluster remaining cells based on DG
ids = find(ASSIGN <1);
idsnow = [];
for i = 1:length(ids)
    sp = SPID_togo(ids(i));
    
    cid = ctogo(ids(i));
    
    if sp == 1
        tmp = find(call(id11)==cid);
    else
        tmp = find(call(id22)==cid);
        tmp = tmp + length(id11);
    end
    idsnow = [idsnow;tmp];
end

normP = [];
for s = 1:length(ids)
    now = cnormPeaks(ids(s),:,2);
    normP = [normP;now(:)' ];
end

clusteragain = 0;
if clusteragain == 1
    
    CLUST = [];
    EVAL = [];
    for k = 2:min(size(normP,1),30)
        [idx,C] = kmeans(normP,k,'Replicates',1000);
        eva1 = evalclusters(normP,idx,'CalinskiHarabasz');%should be high
        eva2 = evalclusters(normP,idx,'DaviesBouldin');%should be low
        eva3 = evalclusters(normP,idx,'Silhouette');%should be high
        CLUST = [CLUST idx];
        evals = [eva1.CriterionValues; eva2.CriterionValues;eva3.CriterionValues];
        EVAL = [EVAL evals];
    end
    figure
    for ev = 1:3
        subplot(1,3,ev)
        plot(2:length(EVAL)+1,EVAL(ev,:))
    end
else
    load(fullfile(path1,'DGclustering_PW'))
end

CC = 6;%4 or 6
% for c = 1:CC
%     figure
%     crr = find(CLUST(:,CC-1)==c);
%     rw = ceil(sqrt(length(crr)+1));
%     cl = ceil((length(crr)+1)/rw);
%     CO = zeros(1,4,6);
%     for cc = 1:length(crr)
%         subplot(rw,cl,cc)
%         htmp = 1-normP(crr(cc),:);
%         %      htmp = mean(htmp,1);
%         htmp=reshape(htmp,4,6);
%         htmp = flipud(fliplr(htmp));
%         if cc == 1
%             CO(1,:,:) = htmp;
%         else
%             htmp2 = zeros(1,4,6); htmp2(1,:,:)=htmp;
%             CO = cat(1,CO,htmp2);
%         end
%         imagesc(htmp)
%         colormap gray
%         axis tight
%         axis off
%     end
%     subplot(rw,cl,rw*cl)
%     imagesc(squeeze(mean(CO,1)))
%     colormap gray
%     axis tight
%     axis off
% end
% 
for c = 1:CC
    crr = find(CLUST(:,CC-1)==c);
    crr2 = ids(crr);
    if length(crr2)>2
        ASSIGN(crr2) = 16+c;
    end
end
%% test
KEEP=ASSIGN;
cls = 19;
crr = find(KEEP==cls);
rw = ceil(sqrt(length(crr)+1));
cl = ceil((length(crr)+1)/rw);
CO = zeros(1,4,6);
figure
for cc = 1:length(crr)
    subplot(rw,cl,cc)
    htmp = 1-cnormPeaks(crr(cc),:,2);
    %      htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    if cc == 1
        CO(1,:,:) = htmp;
    else
        htmp2 = zeros(1,4,6); htmp2(1,:,:)=htmp;
        CO = cat(1,CO,htmp2);
    end
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
end
figure
imagesc(squeeze(mean(CO,1)))
colormap gray
axis tight
axis off
%% let's check for consistency within those new clusters
flash_type = [];
for ii = 1:length(cflash_trans)
    curr = cflash_spike(ii,:);
    tr = cflash_trans(ii,:);
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
cflash_type2 = flash_type;

idn1 = ctogo(find(SPID_togo==1));
idn2 = ctogo(find(SPID_togo==2));
D = diff(cflash); swtch = find(D<0);
D = diff(ctogo); swtch2 = find(D<0);
fl1 = intersect(idn1,cflash(1:swtch));
fl2 = intersect(idn2,cflash(swtch+1:end));
cnt = 0;
clorder = []; origcl = [];
for c = 17:16+CC
    %     if c == 16+CC-1
    %         return
    %     end
    cnt = cnt+1;
    idn = find(ASSIGN==c);
    now1 = idn(find(idn<swtch2+1));%position in ASSIGN
    now2 = idn(find(idn>swtch2));
    now11 = intersect(ctogo(now1),cflash(1:swtch));%ID of cells
    now22 = intersect(ctogo(now2),cflash(swtch+1:end));
    left = setdiff(ctogo(now1),now11);%ID of cells
    left2 = setdiff(ctogo(now2),now22);
    
    tocheck = [];%IDs in cflash
    for n = 1:length(now11)
        tp = find(cflash(1:swtch)==now11(n));
        tocheck = [tocheck;tp];
    end
    for n = 1:length(now22)
        tp = find(cflash(swtch+1:end)==now22(n));
        tocheck = [tocheck;tp+swtch];
    end
    tocheck2 = [];%IDs in ASSIGN or ctogo
    for n = 1:length(now11)
        tp = find(ctogo(1:swtch2)==now11(n));
        tocheck2 = [tocheck2;tp];
    end
    for n = 1:length(now22)
        tp = find(ctogo(swtch2+1:end)==now22(n));
        tocheck2 = [tocheck2;tp+swtch2];
    end
    
    leftcheck = [];%IDs in ASSIGN or ctogo
    for n = 1:length(left)
        tp = find(ctogo(1:swtch2)==left(n));
        leftcheck = [leftcheck;tp];
    end
    for n = 1:length(left2)
        tp = find(ctogo(swtch2+1:end)==left2(n));
        leftcheck = [leftcheck;tp+swtch2];
    end
    
    
    comp = find(cflash_type2(tocheck,1)==1);
    tmp = find(cflash_type2(tocheck,2)<400);
    ontr = intersect(comp,tmp);
    comp = find(cflash_type2(tocheck,1)==1);
    tmp = find(cflash_type2(tocheck,2)>=400);
    onsu = intersect(comp,tmp);
    comp = find(cflash_type2(tocheck,1)==-1);
    tmp = find(cflash_type2(tocheck,2)<400);
    oftr = intersect(comp,tmp);
    comp = find(cflash_type2(tocheck,1)==-1);
    tmp = find(cflash_type2(tocheck,2)>=400);
    ofsu = intersect(comp,tmp);
    comp = find(cflash_type2(tocheck,1)==2);
    tmp = find(cflash_type2(tocheck,2)<400);
    ootr = intersect(comp,tmp);
    comp = find(cflash_type2(tocheck,1)==2);
    tmp = find(cflash_type2(tocheck,2)>=400);
    oosu = intersect(comp,tmp);
    
    allnow = [length(ontr);length(onsu);length(oftr);length(ofsu);length(ootr);length(oosu); length(leftcheck)];
    if allnow(end)<3
        ASSIGN(leftcheck) = 0;
    else
        clorder = [clorder;c];
        origcl = [origcl;c];
    end
    if length(find(allnow(1:end-1)>2))>0
        mv = find(allnow(1:end-1)>2);
        for m = 1:length(mv)
            if mv(m) == 1
                tk = ontr;
            elseif mv(m) == 2
                tk = onsu;
            elseif mv(m) == 3
                tk = oftr;
            elseif mv(m) == 4
                tk = ofsu;
            elseif mv(m) == 5
                tk = ootr;
            elseif mv(m) == 6
                tk = oosu;
                
            end
            idsmv = tocheck2(tk);
            mx = max(unique(ASSIGN));
            clorder = [clorder;mx+1];
            origcl = [origcl;c];
            ASSIGN(idsmv) = mx+1;
        end
    end
    mv = find(allnow(1:end-1)<3);
    for m = 1:length(mv)
        if mv(m) == 1
            tk = ontr;
        elseif mv(m) == 2
            tk = onsu;
        elseif mv(m) == 3
            tk = oftr;
        elseif mv(m) == 4
            tk = ofsu;
        elseif mv(m) == 5
            tk = ootr;
        elseif mv(m) == 6
            tk = oosu;
            
        end
        idsmv = tocheck2(tk);
        ASSIGN(idsmv) = 0;
    end
end
tmp = find(ASSIGN==-1);
ASSIGN(tmp) = 0;

figure
U = unique(ASSIGN);
for u = 1:length(U)
    tmp = find(ASSIGN==U(u));
    plot([0 length(tmp)],[U(u) U(u)],'-k')
    hold on
end

%% overview plot
clear bar
cols = 'rgb';

OVERVIEW = [];
ASSIGNtmp = ASSIGN;
tmp = find(ASSIGN==1); ASSIGNtmp(tmp) = 3;
tmp = find(ASSIGN==2); ASSIGNtmp(tmp) = 1;
tmp = find(ASSIGN==3); ASSIGNtmp(tmp) = 5;
tmp = find(ASSIGN==4); ASSIGNtmp(tmp) = 2;
tmp = find(ASSIGN==5); ASSIGNtmp(tmp) = 4;
ASSIGNF = ASSIGNtmp;
U = unique(ASSIGNtmp);

uu = find(U>16);
% cls = 17:length(uu)+16;
respcl = zeros(length(U)-1,1);
respcl(1:16) = 1:16;
for u = 1:length(uu)
    tmp = find(ASSIGNtmp==U(uu(u)));
    tmp2 = find(clorder == U(uu(u)));
    %     respcl = [respcl;origcl(tmp2)];
    ASSIGNF(tmp) = tmp2+16;
    respcl(tmp2+16) = origcl(tmp2);
end

HASSIGNtmp = H_ASSIGN;
tmp = find(H_ASSIGN==1); HASSIGNtmp(tmp) = 3;
tmp = find(H_ASSIGN==2); HASSIGNtmp(tmp) = 1;
tmp = find(H_ASSIGN==3); HASSIGNtmp(tmp) = 5;
tmp = find(H_ASSIGN==4); HASSIGNtmp(tmp) = 2;
tmp = find(H_ASSIGN==5); HASSIGNtmp(tmp) = 4;
H_ASSIGNF = HASSIGNtmp;


U = unique([ASSIGNF;H_ASSIGNF]);
half = 16;
totH = length( H_TOGO);
D = diff(ctogo);swch = find(D<0);
D = diff(cflash);swch2 = find(D<0);
D = diff(cchrp);swch3 = find(D<0);
D = diff(cR);swch4 = find(D<0);
D = diff(call);swch5 = find(D<0);
totW = length(find(ASSIGNF(1:swch)>0));
totP = length(find(ASSIGNF(swch+1:end)>0));


yh1 = 0.018; yh2 = 0.02; yh3 = 0.012;
% yh1 = 0.022; yh2 = 0.024; yh3 = 0.018;
yst1 = 0.97; yst2 = 0.92; yst3 = 0.87;
x1 = 0.05; x2 = 0.12; x3 = 0.25;x4 = 0.43; x5 = 0.55; x6 = 0.62; x7 = 0.75;x8 = 0.93;
w0 = 0.05; w1 = 0.1; w2 = 0.15;w3 = 0.06;
xnow = [0.6 0.7 0.8 0.9 0.92];
w3 = 0.06; yh4 = 0.033;

figure('position',[538.6000 918.6000 920.8000 926.4000])
ucnt = 0;
for u = 2:length(U)
    wp = find(ASSIGNF==U(u));
    H = find(H_ASSIGNF==U(u));
    
    if u<half+2
        sid = (u-2)*8+1;
        ucnt = ucnt+1;
    else
        sid = (u-half-2)*8+5;
        if u == half+2
            ucnt = 1;
        else
            ucnt = ucnt+1;
        end
    end
    
    
    
    tw = find(wp<=swch);%IDs for ASSIGN/ctogo
    W = wp(tw);
    tp = find(wp>swch);
    P = wp(tp);
    w_togo = ctogo(W); p_togo = ctogo(P);%cellIDs
    h_togo = H; h_togo = intersect(h_togo, H_TOGO);
    
    % percentage
    %     subplot(16,8,sid)
    if u<half+2
        axes('position',[x1 yst1-yh2*(2*ucnt) w0 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x5 yst1-yh2*(2*ucnt) w0 yh3*2])%[0.08 0.9 0.8 0.02])
    end
    perc = [100./totH*length(H) 100./totP*length(P) 100./totW*length(W)];
    abso = [length(H) length(P) length(W)];
    OVERVIEW = [OVERVIEW;abso];
    %     b1 = bar(perc');
    %     hold on
    %     b2 = bar(-abso);
    b=bar([perc' -abso'],'stacked')
    hold on
    text(-0.3,0,int2str(U(u)))
    axis([-0.5 3.5 -50 20])
    axis off
    % axis tight
    
    % flash
    %     subplot(16,8,sid+1)
    if u<half+2
        axes('position',[x2 yst1-yh2*(2*ucnt) w1 yh3])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x6 yst1-yh2*(2*ucnt) w1 yh3])%[0.08 0.9 0.8 0.02])
    end
    FSP = [];
    fid = [];
    for f = 1:length(h_togo)
        tmp = find(H_FLASH_ID==h_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLFLASH  = H_FLASH(fid,:);
    FSP = [FSP;length(fid)];
    
    fid = [];
    for f = 1:length(p_togo)
        tmp = find(cflash(swch+1:end)==p_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp+swch];
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_spike(fid,:)];
    FSP = [FSP;length(fid)];
    
    fid = [];
    for f = 1:length(w_togo)
        tmp = find(cflash(1:swch)==w_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_spike(fid,:)];
    FSP = [FSP;length(fid)];
    
    ALLFLASH2 = ALLFLASH-repmat(mean(ALLFLASH(:,1:200)')',1,size(ALLFLASH,2));
    ALLFLASH2 = ALLFLASH2./repmat(max(ALLFLASH2')',1,size(ALLFLASH2,2));
    
    imagesc(1-ALLFLASH(:,200:5900),[-100 0])
    %       cmap = cmocean('balance','pivot',0);
    %     colormap(gca,cmap)
    colormap(gca,'gray')
    axis off
    axis tight
    
    if u<half+2
        axes('position',[x2 yst1-yh2*(2*ucnt)+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x6 yst1-yh2*(2*ucnt)+yh1 w1 yh3])%[0.08 0.9 0.8 0.02])
    end
    % t1 = 1:FSP(1);
    % plot(mean(ALLFLASH(t1,200:5900),1),'-','color',cols(1))
    %  hold on
    %  t2 = FSP(1)+1:FSP(1)+FSP(2);
    %  plot(mean(ALLFLASH(t2,200:5900),1),'-','color',cols(2))
    %  t3 = FSP(1)+FSP(2)+1:sum(FSP);
    %  plot(mean(ALLFLASH(t3,200:5900),1),'-','color',cols(3))
    plot(mean(ALLFLASH(:,200:5900),1),'-k')
    axis off
    axis([0 5700 0 40])
    
    % chirp
    %     subplot(16,8,sid+2)
    if u<half+2
        axes('position',[x3 yst1-yh2*(2*ucnt) w2 yh3])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x7 yst1-yh2*(2*ucnt) w2 yh3])%[0.08 0.9 0.8 0.02])
    end
    CSP = [];
    fid = [];
    for f = 1:length(h_togo)
        tmp = find(H_CHRP_ID==h_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLCHRP  = H_CHRP(fid,:);
    CSP = [CSP;length(fid)];
    
    fid = [];
    for f = 1:length(p_togo)
        tmp = find(cchrp(swch3+1:end)==p_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp+swch3];
        end
    end
    ALLCHRP  = [ALLCHRP;ccomb_spike(fid,:)];
    CSP = [CSP;length(fid)];
    
    fid = [];
    for f = 1:length(w_togo)
        tmp = find(cchrp(1:swch3)==w_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLCHRP  = [ALLCHRP;ccomb_spike(fid,:)];
    CSP = [CSP;length(fid)];
    
    ALLCHRP2 = ALLCHRP-repmat(mean(ALLCHRP(:,1:200)')',1,size(ALLCHRP,2));
    ALLCHRP2 = ALLCHRP2./repmat(max(ALLCHRP2')',1,size(ALLCHRP2,2));
    ALLCHRP3 = ALLCHRP./repmat(max(ALLCHRP')',1,size(ALLCHRP,2));
    imagesc(1-ALLCHRP3)
    colormap gray
    axis off
    axis tight
    
    if u<half+2
        axes('position',[x3 yst1-yh2*(2*ucnt)+yh1 w2 yh3])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x7 yst1-yh2*(2*ucnt)+yh1 w2 yh3])%[0.08 0.9 0.8 0.02])
    end
    t1 = 1:CSP(1);
    plot(mean(ALLCHRP3(t1,:),1),'-','color',cols(1))
    hold on
    t2 = CSP(1)+1:CSP(1)+CSP(2);
    plot(mean(ALLCHRP3(t2,:),1),'-','color',cols(2))
    t3 = CSP(1)+CSP(2)+1:sum(CSP);
    plot(mean(ALLCHRP3(t3,:),1),'-','color',cols(3))
    axis off
    axis tight
    
    % DG
    %     subplot(16,8,sid+3)
    if u<half+2
        axes('position',[x4 yst1-yh2*(2*ucnt) w3 yh2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x8 yst1-yh2*(2*ucnt) w3 yh2])%[0.08 0.9 0.8 0.02])
    end
    DSP = [];
    fid = [];
    for f = 1:length(h_togo)
        tmp = find(H_TOGO==h_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLDG  = H_DG(fid,:,:);
    DSP = [DSP;length(fid)];
    
    fid = [];
    for f = 1:length(p_togo)
        tmp = find(ctogo(swch2+1:end)==p_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp+swch2];
        end
    end
    ALLDG  = cat(1,ALLDG,cnormPeaks(fid,:,:));
    DSP = [DSP;length(p_togo)];
    
    fid = [];
    for f = 1:length(w_togo)
        tmp = find(ctogo(1:swch)==w_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLDG  = cat(1,ALLDG,cnormPeaks(fid,:,:));
    DSP = [DSP;length(w_togo)];
    
    normP = [];
    for n = 1:size(ALLDG,1)
        now = ALLDG(n,:,2);
        normP = [normP;now(:)' ];
    end
    htmp = 1-normP; htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    hold on
    rectangle('position',[0.5 0.5  6 4],'edgecolor','k','facecolor','none')
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    
    
end

OVERVIEW = [OVERVIEW respcl];
%% overview plot2
noflash = [9:16 32];
yh1 = 0.016; yh2 = 0.018; yh3 = 0.012;
% yh1 = 0.022; yh2 = 0.024; yh3 = 0.018;
yst1 = 0.97; yst2 = 0.92; yst3 = 0.87;
x1 = 0.1; x2 = 0.135; x3 = 0.19;x4 = 0.26; xx1 = 0.33;xx2 = 0.34;xx3 = 0.45;xx4 = 0.39;
x5 = 0.6; x6 = 0.635; x7 = 0.69;x8 = 0.76; xx5 = 0.83; xx6 = 0.84; xx7 = 0.95; xx8 = 0.89;
w0 = 0.025; w1 = 0.04; w2 = 0.08;w3 = 0.06; ww1 = 0.04;
xnow = [0.6 0.7 0.8 0.9 0.92];
w3 = 0.06; yh4 = 0.033;

spat=[100 100 100 100 200 200 200 200 500 500 500 500 1000 1000 1000 1000 2000 2000 2000 2000 4000 4000 4000 4000];
freq=[1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8];
SP = unique(spat); FR =unique(freq);

clustertext = {'put. ON parasol','put. OFF parasol','put. ON midget','put. OFF midget',...,
    'put. ON-OFF bistrat.','H ON full-field','H OFF full-field','H full-field highTF',...,
    'H highTF/SF','H speed tuned','H highTF','H speed tuned 2','H highSF','H mediumTF', 'H 4Hz','H mediumTF/SF',...,
    'ONtr very lowTF','ONtr low TF','ONtr broad','ONsu very lowTF','ONsu low TF','ONsu broad','ON narrowTF/SF','ON high TF','ON very high TF',...,
    'OFFtr very lowTF','OFFtr low TF','OFFtr no contrast','OFFsu broad',...,
    'ON-OFF very lowTF','ON-OFF mixed','ON-OFF lowTF','ON-OFF narrowSF','ON-OFF highTF','(ON) 2000um'};

figure('position',[538.6000 918.6000 920.8000 926.4000])
ucnt = 0;
% Unow = [0 1:16 17 18 23 24 26 27 30 33 35 19 20 25 28 21 22 29 31 34 32];
Unow = [0 1:16 26 17 23 27 18 24 30 33 35 28 19 20 25 29 21 22 31 34 32];
PRC = [];
for u = 2:length(Unow)
    wp = find(ASSIGNF==Unow(u));
    H = find(H_ASSIGNF==Unow(u));
    
    if u<half+2
        sid = (u-2)*8+1;
        ucnt = ucnt+1;
    else
        sid = (u-half-2)*8+5;
        if u == half+2
            ucnt = 1;
        else
            ucnt = ucnt+1;
        end
    end
    
    
    
    tw = find(wp<=swch);%IDs for ASSIGN/ctogo
    W = wp(tw);
    tp = find(wp>swch);
    P = wp(tp);
    w_togo = ctogo(W); p_togo = ctogo(P);%cellIDs
    h_togo = H; h_togo = intersect(h_togo, H_TOGO);
    
    % percentage
    %     subplot(16,8,sid)
    if u<half+2
        axes('position',[x1 yst1-yh2*(2*ucnt) w0 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x5 yst1-yh2*(2*ucnt) w0 yh3*2])%[0.08 0.9 0.8 0.02])
    end
    perc = [100./totH*length(H) 100./totP*length(P) 100./totW*length(W)];
    PRC = [PRC;perc];
    abso = [length(H) length(P) length(W)];
    %     OVERVIEW = [OVERVIEW;abso];
    %     b1 = bar(perc');
    %     hold on
    %     b2 = bar(-abso);
    b=bar(perc)
    hold on
    text(-13,22,['cluster ',int2str(u-1)],'fontsize',7,'fontweight','bold')
    text(-13,10,clustertext{u-1},'fontsize',6)
    axis([-0.5 3.5 0 25])
    axis off
    % axis tight
    
    % flash
    %     subplot(16,8,sid+1)
    %    % latencies
    %       if u<half+2
    %         axes('position',[xx1 yst1-yh2*(2*ucnt) ww1 yh3*2])%[0.08 0.9 0.8 0.02])
    %     else
    %         axes('position',[xx5 yst1-yh2*(2*ucnt) ww1 yh3*2])%[0.08 0.9 0.8 0.02])
    %       end
    FSP = [];
    fid = [];
    for f = 1:length(h_togo)
        tmp = find(H_FLASH_ID==h_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLFLASH  = H_FLASH(fid,:);
    FSP = [FSP;length(fid)];
    TRANS = [];
    
    fid = [];tkn = [];flt = [];
    for f = 1:length(p_togo)
        tmp = find(cflash(swch+1:end)==p_togo(f));
        if ~isempty(tmp)
            tkn = [tkn;p_togo(f)];
            fid = [fid;tmp+swch];
            flt = [flt;cflash_type2(tmp+swch,2)];
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_spike(fid,:)];
    %     [N,edg] =histcounts(flt,0:100:2000);
    %     edg = edg+(edg(2)-edg(1))/2; edg = edg(1:end-1);
    %     plot(edg,N)
    hold on
    %     FSP = [FSP;length(fid)];
    fid2 = [];
    for f = 1:length(p_togo)
        if isempty(find(tkn==p_togo(f)))
            tmp = find(cR(swch4+1:end)==p_togo(f));
            if ~isempty(tmp)
                fid2 = [fid2;tmp+swch4];
            end
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_all(fid2,:)];
    FSP = [FSP;length(fid)+length(fid2)];
    
    fid = [];tkn = [];flt = [];
    for f = 1:length(w_togo)
        tmp = find(cflash(1:swch)==w_togo(f));
        if ~isempty(tmp)
            tkn = [tkn;w_togo(f)];
            fid = [fid;tmp];
            flt = [flt;cflash_type2(tmp,2)];
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_spike(fid,:)];
    %      [N,edg] =histcounts(flt,0:100:2000);
    %     edg = edg+(edg(2)-edg(1))/2; edg = edg(1:end-1);
    %     plot(edg,N)
    %     FSP = [FSP;length(fid)];
    fid2 = [];
    for f = 1:length(w_togo)
        if isempty(find(tkn==w_togo(f)))
            tmp = find(cR(1:swch4)==w_togo(f));
            if ~isempty(tmp)
                fid2 = [fid2;tmp];
            end
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_all(fid2,:)];
    FSP = [FSP;length(fid)+length(fid2)];
    
    ALLFLASH2 = ALLFLASH-repmat(mean(ALLFLASH(:,1:200)')',1,size(ALLFLASH,2));
    ALLFLASH2 = ALLFLASH2./repmat(max(ALLFLASH2')',1,size(ALLFLASH2,2));
    
    
    if u<half+2
        axes('position',[x2 yst1-yh2*(2*ucnt) w1 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x6 yst1-yh2*(2*ucnt) w1 yh3*2])%[0.08 0.9 0.8 0.02])
    end
    % t1 = 1:FSP(1);
    % plot(mean(ALLFLASH(t1,200:5900),1),'-','color',cols(1))
    %  hold on
    %  t2 = FSP(1)+1:FSP(1)+FSP(2);
    %  plot(mean(ALLFLASH(t2,200:5900),1),'-','color',cols(2))
    %  t3 = FSP(1)+FSP(2)+1:sum(FSP);
    %  plot(mean(ALLFLASH(t3,200:5900),1),'-','color',cols(3))
    if isempty(find(noflash==Unow(u)))
        plot(mean(ALLFLASH(:,200:5900),1),'-k')
    else
        plot(mean(ALLFLASH(:,200:5900),1),'-','color',[0.8 0.8 0.8])
    end
    axis off
    axis([0 5700 0 40])
    
    
    
    
    
    % chirp
    %     subplot(16,8,sid+2)
    
    CSP = [];
    fid = []; tkn = [];
    for f = 1:length(h_togo)
        tmp = find(H_CHRP_ID==h_togo(f));
        if ~isempty(tmp)
            tkn = [tkn;h_togo(f)];
            fid = [fid;tmp];
        end
    end
    ALLCHRP  = H_CHRP(fid,:);
    Hch1 = H_ChrpIndex(fid);
    %     CSP = [CSP;length(fid)];
    fid2 = [];
    for f = 1:length(h_togo)
        if isempty(find(tkn==h_togo(f)))
            tmp = find(H_CHRO_ALLID==h_togo(f));
            if ~isempty(tmp)
                fid2 = [fid2;tmp];
            end
        end
    end
    ALLCHRP  = H_CHRP_ALL(fid2,:);
    CSP = [CSP;length(fid)+length(fid2)];
    
    fid = [];
    for f = 1:length(p_togo)
        tmp = find(cchrp(swch3+1:end)==p_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp+swch3];
        end
    end
    ALLCHRP  = [ALLCHRP;ccomb_spike(fid,:)];
    CSP = [CSP;length(fid)];
        Pch1 = cchirp_index(fid);
        
    fid = []; 
    for f = 1:length(w_togo)
        tmp = find(cchrp(1:swch3)==w_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
        
    end
    ALLCHRP  = [ALLCHRP;ccomb_spike(fid,:)];
    CSP = [CSP;length(fid)];
     Wch1 = cchirp_index(fid);
    
    ALLCHRP2 = ALLCHRP-repmat(mean(ALLCHRP(:,1:200)')',1,size(ALLCHRP,2));
    ALLCHRP2 = ALLCHRP2./repmat(max(ALLCHRP2')',1,size(ALLCHRP2,2));
    ALLCHRP3 = ALLCHRP./repmat(max(ALLCHRP')',1,size(ALLCHRP,2));
    
    % chirp index
     if u<half+2
        axes('position',[xx1 yst1-yh2*(2*ucnt) ww1 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[xx5 yst1-yh2*(2*ucnt) ww1 yh3*2])%[0.08 0.9 0.8 0.02])
    end
    if ~isempty(Hch1)
        plotSpread(Hch1,'distributionIdx',repmat(1,length(Hch1),1),'distributionMarkers',...
            {'.'},'distributionColors',{[1 0 0 0]},'xyOri','flipped')
        hold on
         if length(Hch1)>2
        plot([mean(Hch1) mean(Hch1)],[0.2 1.8],'-','color','r','linewidth',1)
         end
    end
    if ~isempty(Pch1)
        plotSpread(Pch1,'distributionIdx',repmat(1,length(Pch1),1),'distributionMarkers',...
            {'.'},'distributionColors',{[0 1 0 0]},'xyOri','flipped')
        hold on
         if length(Pch1)>2
        plot([mean(Pch1) mean(Pch1)],[0.2 1.8],'-','color','g','linewidth',1)
         end
    end
     if ~isempty(Wch1)
        plotSpread(Wch1,'distributionIdx',repmat(1,length(Wch1),1),'distributionMarkers',...
            {'.'},'distributionColors',{[0 0 1 0]},'xyOri','flipped')
        hold on
         if length(Wch1)>2
        plot([mean(Wch1) mean(Wch1)],[0.2 1.8],'-','color','b','linewidth',1)
         end
    end
     axis([-1 1 0 2])
    plot([0 0],[0 2],'--','color',[0.3 0.3 0.3])
    set(gca,'xtick',[-1 0 1],'xticklabel',[])
    set(gca,'ytick',[],'ycolor','w')
    
    
    if u<half+2
        axes('position',[x3 yst1-yh2*(2*ucnt) w2 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x7 yst1-yh2*(2*ucnt) w2 yh3*2])%[0.08 0.9 0.8 0.02])
    end
    %     t1 = 1:CSP(1);
    %     plot(mean(ALLCHRP3(t1,:),1),'-','color',cols(1))
    %     hold on
    %     t2 = CSP(1)+1:CSP(1)+CSP(2);
    %     plot(mean(ALLCHRP3(t2,:),1),'-','color',cols(2))
    %     t3 = CSP(1)+CSP(2)+1:sum(CSP);
    %     plot(mean(ALLCHRP3(t3,:),1),'-','color',cols(3))
    plot(mean(ALLCHRP3,1),'-','color','k')
    axis off
    axis([0 16502 0 1])
    
    % DG
    %     subplot(16,8,sid+3)
    
    
    DSP = [];
    fid = [];
    for f = 1:length(h_togo)
        tmp = find(H_TOGO==h_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLDG  = H_DG(fid,:,:);
    DSP = [DSP;length(fid)];
    
    
    
    [cX,cY,oX,oY,dX,dY] = get_GaussHeatmaps(H_DG(fid,:,:));
    dat2 = [];
    for d = 1:length(cX)
        now = cX(d);
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
    Hsp = 4000-dat2;
    
    
    cY = 5-cY;
    dat2 = [];
    for d = 1:length(cY)
        now = cY(d);
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
    Htm = dat2;
    
    
    fid = [];
    for f = 1:length(p_togo)
        tmp = find(ctogo(swch2+1:end)==p_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp+swch2];
        end
    end
    ALLDG  = cat(1,ALLDG,cnormPeaks(fid,:,:));
    DSP = [DSP;length(p_togo)];
    [cX,cY,oX,oY,dX,dY] = get_GaussHeatmaps(cnormPeaks(fid,:,:));
    dat2 = [];
    for d = 1:length(cX)
        now = cX(d);
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
    Psp = 4000-dat2;
    cY = 5-cY;
    dat2 = [];
    for d = 1:length(cY)
        now = cY(d);
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
    Ptm = dat2;
    
    fid = [];
    for f = 1:length(w_togo)
        tmp = find(ctogo(1:swch)==w_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLDG  = cat(1,ALLDG,cnormPeaks(fid,:,:));
    DSP = [DSP;length(w_togo)];
    [cX,cY,oX,oY,dX,dY] = get_GaussHeatmaps(cnormPeaks(fid,:,:));
    dat2 = [];
    for d = 1:length(cX)
        now = cX(d);
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
    Wsp = 4000-dat2;
    cY = 5-cY;
    dat2 = [];
    for d = 1:length(cY)
        now = cY(d);
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
    Wtm = dat2;
    
    % spatial index
    if u<half+2
        axes('position',[xx3 yst1-yh2*(2*ucnt) ww1 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[xx7 yst1-yh2*(2*ucnt) ww1 yh3*2])%[0.08 0.9 0.8 0.02])
    end
    if ~isempty(Hsp)
        plotSpread(Hsp,'distributionIdx',repmat(1,length(Hsp),1),'distributionMarkers',...
            {'.'},'distributionColors',{[1 0 0 0]},'xyOri','flipped')
        hold on
         if length(Hsp)>2
        plot([mean(Hsp) mean(Hsp)],[0.2 1.8],'-','color','r','linewidth',1)
         end
    end
    if ~isempty(Psp)
        plotSpread(Psp,'distributionIdx',repmat(1,length(Psp),1),'distributionMarkers',...
            {'.'},'distributionColors',{[0 1 0 0]},'xyOri','flipped')
        hold on
         if length(Psp)>2
        plot([mean(Psp) mean(Psp)],[0.2 1.8],'-','color','g','linewidth',1)
         end
    end
    if ~isempty(Wsp)
        plotSpread(Wsp,'distributionIdx',repmat(1,length(Wsp),1),'distributionMarkers',...
            {'.'},'distributionColors',{[0 0 1 0]},'xyOri','flipped')
        hold on
        if length(Wsp)>2
        plot([mean(Wsp) mean(Wsp)],[0.2 1.8],'-','color','b','linewidth',1)
        end
    end
    axis([0 4000 0 2])
    plot([2000 2000],[0 2],'--','color',[0.3 0.3 0.3])
    set(gca,'xtick',[0 2000 4000],'xticklabel',[])
    set(gca,'ytick',[],'ycolor','w')
    
    
    if u<half+2
        axes('position',[xx4 yst1-yh2*(2*ucnt) ww1 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[xx8 yst1-yh2*(2*ucnt) ww1 yh3*2])%[0.08 0.9 0.8 0.02])
    end
    if ~isempty(Htm)
        plotSpread(Htm,'distributionIdx',repmat(1,length(Htm),1),'distributionMarkers',...
            {'.'},'distributionColors',{[1 0 0 0]},'xyOri','flipped')
        hold on
        if length(Htm)>2
        plot([mean(Htm) mean(Htm)],[0.2 1.8],'-','color','r','linewidth',1)
        end
    end
    if ~isempty(Ptm)
        plotSpread(Ptm,'distributionIdx',repmat(1,length(Ptm),1),'distributionMarkers',...
            {'.'},'distributionColors',{[0 1 0 0]},'xyOri','flipped')
        hold on
        if length(Ptm)>2
        plot([mean(Ptm) mean(Ptm)],[0.2 1.8],'-','color','g','linewidth',1)
        end
    end
    if ~isempty(Wtm)
        plotSpread(Wtm,'distributionIdx',repmat(1,length(Wtm),1),'distributionMarkers',...
            {'.'},'distributionColors',{[0 0 1 0]},'xyOri','flipped')
        hold on
        if length(Wtm)>2
        plot([mean(Wtm) mean(Wtm)],[0.2 1.8],'-','color','b','linewidth',1)
        end
    end
    plot([4 4],[0 2],'--','color',[0.3 0.3 0.3])
    set(gca,'xtick',[0 4 8],'xticklabel',[])
    set(gca,'ytick',[],'ycolor','w')
    axis([0 8 0 2])
    %      axis off
    
    
    normP = [];
    for n = 1:size(ALLDG,1)
        now = ALLDG(n,:,2);
        normP = [normP;now(:)' ];
    end
    htmp = 1-normP; htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    
    if u<half+2
        axes('position',[x4 yst1-yh2*(2*ucnt) w3 yh2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x8 yst1-yh2*(2*ucnt) w3 yh2])%[0.08 0.9 0.8 0.02])
    end
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    hold on
    rectangle('position',[0.5 0.5  6 4],'edgecolor','k','facecolor','none')
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    
    
end
exportgraphics(gcf,fullfile(pathSave,'NEW_Fig5_clustering.eps'),'Resolution',600)
%% overview plot3
noflash = [9:16 32];
yh1 = 0.018; yh2 = 0.02; yh3 = 0.012;
% yh1 = 0.022; yh2 = 0.024; yh3 = 0.018;
yst1 = 0.97; yst2 = 0.92; yst3 = 0.87;
x1 = 0.05; x2 = 0.12; x3 = 0.185;x4 = 0.34; 
x5 = 0.55; x6 = 0.62; x7 = 0.685;x8 = 0.84;
w0 = 0.05; w1 = 0.05; w2 = 0.15;w3 = 0.06;
xnow = [0.6 0.7 0.8 0.9 0.92];
w3 = 0.06; yh4 = 0.033;

figure('position',[538.6000 918.6000 920.8000 926.4000])
ucnt = 0;
% Unow = [0 1:16 17 18 23 24 26 27 30 33 35 19 20 25 28 21 22 29 31 34 32];
Unow = [0 1:16 26 17 23 27 18 24 30 33 35 28 19 20 25 29 21 22 31 34 32];
PRC = [];
for u = 2:length(Unow)
    wp = find(ASSIGNF==Unow(u));
    H = find(H_ASSIGNF==Unow(u));
    
    if u<half+2
        sid = (u-2)*8+1;
        ucnt = ucnt+1;
    else
        sid = (u-half-2)*8+5;
        if u == half+2
            ucnt = 1;
        else
            ucnt = ucnt+1;
        end
    end
    
    
    
    tw = find(wp<=swch);%IDs for ASSIGN/ctogo
    W = wp(tw);
    tp = find(wp>swch);
    P = wp(tp);
    w_togo = ctogo(W); p_togo = ctogo(P);%cellIDs
    h_togo = H; h_togo = intersect(h_togo, H_TOGO);
    
    
   
   
    
    
    % flash
   
    FSP = [];
    fid = [];
    for f = 1:length(h_togo)
        tmp = find(H_FLASH_ID==h_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLFLASH  = H_FLASH(fid,:);
    FSP = [FSP;length(fid)];
    
    
    fid = [];tkn = [];
    for f = 1:length(p_togo)
        tmp = find(cflash(swch+1:end)==p_togo(f));
        if ~isempty(tmp)
            tkn = [tkn;p_togo(f)];
            fid = [fid;tmp+swch];
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_spike(fid,:)];
%     FSP = [FSP;length(fid)];
    fid2 = [];
    for f = 1:length(p_togo)
        if isempty(find(tkn==p_togo(f)))
        tmp = find(cR(swch4+1:end)==p_togo(f));
        if ~isempty(tmp)
            fid2 = [fid2;tmp+swch4];
        end
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_all(fid2,:)];
    FSP = [FSP;length(fid)+length(fid2)];
    
    fid = [];tkn = [];
    for f = 1:length(w_togo)
        tmp = find(cflash(1:swch)==w_togo(f));
        if ~isempty(tmp)
              tkn = [tkn;w_togo(f)];
            fid = [fid;tmp];
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_spike(fid,:)];
%     FSP = [FSP;length(fid)];
 fid2 = [];
    for f = 1:length(w_togo)
        if isempty(find(tkn==w_togo(f)))
        tmp = find(cR(1:swch4)==w_togo(f));
        if ~isempty(tmp)
            fid2 = [fid2;tmp];
        end
        end
    end
    ALLFLASH  = [ALLFLASH;cflash_all(fid2,:)];
    FSP = [FSP;length(fid)+length(fid2)];
    
    ALLFLASH2 = ALLFLASH-repmat(mean(ALLFLASH(:,1:200)')',1,size(ALLFLASH,2));
    ALLFLASH2 = ALLFLASH2./repmat(max(ALLFLASH2')',1,size(ALLFLASH2,2));
        if u<half+2
        axes('position',[x2 yst1-yh2*(2*ucnt) w1 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x6 yst1-yh2*(2*ucnt) w1 yh3*2])%[0.08 0.9 0.8 0.02])
        end
    imagesc(1-ALLFLASH(:,200:5900),[-100 0])
    %       cmap = cmocean('balance','pivot',0);
    %     colormap(gca,cmap)
    colormap(gca,'gray')
    axis off
    axis tight
    
  
   
     if u<half+2
        axes('position',[x1 yst1-yh2*(2*ucnt) w0 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x5 yst1-yh2*(2*ucnt) w0 yh3*2])%[0.08 0.9 0.8 0.02])
    end
    t1 = 1:FSP(1);
    plot(mean(ALLFLASH(t1,200:5900),1),'-','color',cols(1))
     hold on
     t2 = FSP(1)+1:FSP(1)+FSP(2);
     plot(mean(ALLFLASH(t2,200:5900),1),'-','color',cols(2))
     t3 = FSP(1)+FSP(2)+1:sum(FSP);
     plot(mean(ALLFLASH(t3,200:5900),1),'-','color',cols(3))
%     if isempty(find(noflash==Unow(u)))
%     plot(mean(ALLFLASH(:,200:5900),1),'-k')
%     else
%         plot(mean(ALLFLASH(:,200:5900),1),'-','color',[0.8 0.8 0.8])
%     end
    axis off
    axis([0 5700 0 40])
    
    % chirp
    %     subplot(16,8,sid+2)
   
    CSP = [];
    fid = []; tkn = [];
    for f = 1:length(h_togo)
        tmp = find(H_CHRP_ID==h_togo(f));
        if ~isempty(tmp)
              tkn = [tkn;h_togo(f)];
            fid = [fid;tmp];
        end
    end
    ALLCHRP  = H_CHRP(fid,:);
%     CSP = [CSP;length(fid)];
       fid2 = []; 
    for f = 1:length(h_togo)
        if isempty(find(tkn==h_togo(f)))
        tmp = find(H_CHRO_ALLID==h_togo(f));
        if ~isempty(tmp)
            fid2 = [fid2;tmp];
        end
        end
    end
    ALLCHRP  = H_CHRP_ALL(fid2,:);
    CSP = [CSP;length(fid)+length(fid2)];
    
    fid = [];
    for f = 1:length(p_togo)
        tmp = find(cchrp(swch3+1:end)==p_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp+swch3];
        end
    end
    ALLCHRP  = [ALLCHRP;ccomb_spike(fid,:)];
    CSP = [CSP;length(fid)];
    
    fid = [];
    for f = 1:length(w_togo)
        tmp = find(cchrp(1:swch3)==w_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLCHRP  = [ALLCHRP;ccomb_spike(fid,:)];
    CSP = [CSP;length(fid)];
    
    ALLCHRP2 = ALLCHRP-repmat(mean(ALLCHRP(:,1:200)')',1,size(ALLCHRP,2));
    ALLCHRP2 = ALLCHRP2./repmat(max(ALLCHRP2')',1,size(ALLCHRP2,2));
    ALLCHRP3 = ALLCHRP./repmat(max(ALLCHRP')',1,size(ALLCHRP,2));
   
    if u<half+2
        axes('position',[x3 yst1-yh2*(2*ucnt) w2 yh3*2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x7 yst1-yh2*(2*ucnt) w2 yh3*2])%[0.08 0.9 0.8 0.02])
    end
%     t1 = 1:CSP(1);
%     plot(mean(ALLCHRP3(t1,:),1),'-','color',cols(1))
%     hold on
%     t2 = CSP(1)+1:CSP(1)+CSP(2);
%     plot(mean(ALLCHRP3(t2,:),1),'-','color',cols(2))
%     t3 = CSP(1)+CSP(2)+1:sum(CSP);
%     plot(mean(ALLCHRP3(t3,:),1),'-','color',cols(3))
 plot(mean(ALLCHRP3,1),'-','color','k')
    axis off
    axis([0 16502 0 1])
    
    % DG
    %     subplot(16,8,sid+3)
    if u<half+2
        axes('position',[x4 yst1-yh2*(2*ucnt) w3 yh2])%[0.08 0.9 0.8 0.02])
    else
        axes('position',[x8 yst1-yh2*(2*ucnt) w3 yh2])%[0.08 0.9 0.8 0.02])
    end
    DSP = [];
    fid = [];
    for f = 1:length(h_togo)
        tmp = find(H_TOGO==h_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLDG  = H_DG(fid,:,:);
    DSP = [DSP;length(fid)];
    
    fid = [];
    for f = 1:length(p_togo)
        tmp = find(ctogo(swch2+1:end)==p_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp+swch2];
        end
    end
    ALLDG  = cat(1,ALLDG,cnormPeaks(fid,:,:));
    DSP = [DSP;length(p_togo)];
    
    fid = [];
    for f = 1:length(w_togo)
        tmp = find(ctogo(1:swch)==w_togo(f));
        if ~isempty(tmp)
            fid = [fid;tmp];
        end
    end
    ALLDG  = cat(1,ALLDG,cnormPeaks(fid,:,:));
    DSP = [DSP;length(w_togo)];
    
    normP = [];
    for n = 1:size(ALLDG,1)
        now = ALLDG(n,:,2);
        normP = [normP;now(:)' ];
    end
    htmp = 1-normP; htmp = mean(htmp,1);
    htmp=reshape(htmp,4,6);
    htmp = flipud(fliplr(htmp));
    imagesc(htmp)
    colormap gray
    axis tight
    axis off
    hold on
    rectangle('position',[0.5 0.5  6 4],'edgecolor','k','facecolor','none')
    set(gca,'plotboxaspectratio',[1 1/6*4 1])
    
    
end

OVERVIEW = [OVERVIEW respcl];
