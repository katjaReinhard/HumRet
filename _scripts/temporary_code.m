
%% adjustable variables
clear
close all

Dtype='h';
% superpath='S:\data\Katja\allspecies';
superpath='S:\data\Katja\allspecies';
%% data

switch Dtype
    
    case 'h'
        % Human data
        dates={'20120210','20120224','20120320c','20120608b','20130124a','20130211a',...,
            '20130211c','20130211d','20130510a','20130510d','20130510e','20130718b','20130718c','20131210a','20131217f'};
        mainpath='S:\data\Katja\Human';
        starts=[255 245 8 91 1 97 192 97 1 1 1 173 1 1 88]; %index as in unit files
        stops=[374 365 277 280 191 190 449 285 246 164 246 429 172 86 260];
        NDs=[4 4 4 4 4 3 4 4 4 4 4 3 4 4 3];
    case 'p'
        % Pig data
        dates={'20121127','20121211a','20130613a','20130613b','20130613c','20140723d','20120213','20120319a','20120319b'};
        mainpath='S:\data\Katja\Pig';
        starts=[1 1 1 1 1 156 383 1 361];
        stops=[190 190 164 164 164 285 502 270 630];
        NDs=[4 4 4 4 4 4 3 3 3];
        
    case 'm'
        % Mouse data
        dates={'20120126','20120127','20120221','20120222'};
        mainpath='S:\data\Katja\OPA';
        starts=[363 363 363 363];
        stops=[543 543 543 543];
        NDs=[4 4 4 4];
end
%% 10) !!!!!!!!!!!!!!!!!!!!!!!!work on this for human data!!
% col18) fourier peak
load(fullfile(superpath,[Dtype,'_info']))
load(fullfile(superpath,[Dtype,'_DG']))
cd(superpath)


DGrates=cell(length(info),1);
togo=find(goodones>0);
for i=1:length(togo)
    i
    currID=togo(i);
    currd=info{currID,1};
    curru=info{currID,2};
    currstart=info{currID,4};
    currstop=info{currID,5};
    
    hpath=fullfile(mainpath,currd,'HEKA');
    hlist=dir(fullfile(hpath,'*phys'));
    upath=fullfile(mainpath,currd,'units');
    ulist=dir(fullfile(upath,'*mat'));
    
    found=0;m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(ulist(m).name,['unit_',curru]))
            found=1;
        end
    end
    load(fullfile(upath,ulist(m).name))
    
    hekalist=unit{1,2}(:,1);
    
    DG=[];
    for f=currstart:currstop
        if ~isempty(regexp(hekalist{f},'DriftGrat_Contr1'))
            DG=[DG;1];
        else
            DG=[DG;0];
        end
    end
    
    goodfiles=currstart:currstop;
    goodfiles=goodfiles(DG==1);
    tmp=regexp(hekalist{goodfiles(1)},'_');
    tmp=tmp(2);
    names=hekalist(goodfiles);
    clear names2
    for n=1:length(names)
        names2{n}=names{n}(tmp:end);
    end
    DGlist=unique(names2)';
    
    clear real
    cnt=0; freqs=[]; spats=[];
    for dg=1:length(DGlist)
        if isempty(regexp(DGlist{dg},'speed_0.25'))
            cnt=cnt+1;
            real{cnt}=DGlist{dg};
            tmp=regexp(DGlist{dg},'speed_');
            tmp1=tmp+6;
            tmp=regexp(DGlist{dg}(tmp1:end),'%');
            tmp2=tmp1+tmp-2;
            ff=DGlist{dg}(tmp1:tmp2);
            freqs=[freqs;str2num(ff)];
            tmp=regexp(DGlist{dg},'um%');
            tmp2=tmp-1;
            tmp=regexp(DGlist{dg}(1:tmp2),'_');
            tmp1=tmp(end)+1;
            sp=DGlist{dg}(tmp1:tmp2);
            spats=[spats;str2num(sp)];
        end
    end
    DGlist=real;
    
   
    rates=[];
    for g=1:length(DGlist)
        currfile=DGlist{g};
        DG=[];
        for f=goodfiles
            if ~isempty(regexp(hekalist{f},currfile))
                DG=[DG;1];
            else
                DG=[DG;0];
            end
        end
        goodnow=goodfiles(DG==1);
        if g==1
            first=goodnow(1);
            prot=read_header_field_heka(hpath,hekalist{first},'Stimulus Protocol');
            if prot(end,1)>19000
                gofor=0;
            else
                gofor=1;
            end
        end
        
        if gofor==1
            allrate=[];
            for f=goodnow
                spikes=unit{1,2}{f,2};
                rate=MEA_spikerates(spikes,40,14000);
                allrate=[allrate;rate];
            end
            meanrate=mean(allrate);
            rates=[rates;meanrate];
            
            stimFreq=freqs(g);
            
           
        end
    end
    
    
    if gofor==1
      
        for r=1:size(rates,1)
            DGrates{currID,r}=rates(r,:);
        end
      
    
    end
end
save(fullfile(superpath,[Dtype,'_DG']),'stimInfo','DGfouriers','DGrates')