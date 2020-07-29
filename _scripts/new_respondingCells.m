%% TABLE: which cell responds to what stimulus
% 1) flash, 2) LF, 3) DG, 4) chirp, 5) 6vel black, 6) 6vel white, 7) DS

clear
close all
superpath='/media/NERFFS01/Users/Katja/human_retina';

Dtype='h';

% idlist=[8 9 0 19 16 17 10/11];

d=1
load(fullfile(superpath,[Dtype(d),'_info']))
load(fullfile(superpath,[Dtype(d),'_normPeaks']))
response=zeros(size(info,1),7);
pols=zeros(size(info,1),2);
pols=pols-2;
dates=unique(info(:,1));
did=zeros(size(info,1),1);


for i=1:size(info,1)
    found=0; m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(dates{m},info{i,1}))
            found=1;
        end
    end
    did(i)=m;
    
    
    if ~isempty(info{i,8}) %flash
        if info{i,8}~=0 && ~isnan(info{i,8})
            response(i,1)=1;
            if info{i,8}==1;
                pols(i,1)=1;
            elseif info{i,8}==-1
                pols(i,1)=-1;
            elseif info{i,8}==2
                pols(i,1)=0;
            end
        end
    end
    if ~isempty(info{i,9}) %LF
        if info{i,9}~=0 && ~isnan(info{i,9})
            response(i,2)=1;
            if info{i,9}==1;
                pols(i,2)=1;
            elseif info{i,9}==-1
                pols(i,2)=-1;
            end
        end
    end
    if ~isempty(info{i,19})
        if info{i,19}~=0 && ~isnan(info{i,19})
            response(i,4)=1;
        end
    end
    if ~isempty(info{i,16})
        if info{i,16}~=0 && ~isnan(info{i,16})
            response(i,5)=1;
        end
    end
    if ~isempty(info{i,17})
        if info{i,17}~=0 && ~isnan(info{i,17})
            response(i,6)=1;
        end
    end
    
    %         if ~isempty(info{i,10}) ||  ~isempty(info{i,11})
    %             if info{i,10}~=0 || info{i,11}~=0
    %                 if ~isnan(info{i,10}) || ~isnan(info{i,11})
    %                     response(i,7)=1;
    %                 end
    %             end
    %         end
    if ~isempty(info{i,12})
        if info{i,12}~=0 && ~isnan(info{i,12})
            response(i,7)=1;
        end
    end
    
    if d==2
        if ~isempty(find(did(i)==[1 2 3 7]));
            response(i,4)=-1;
        end
        if ~isempty(find(did(i)==7));
            response(i,3)=-1;
        end
        if ~isempty(find(did(i)==9));
            response(i,[3 5 6])=-1;
        end
    end
    
    if d==3
        if ~isempty(find(did(i)==[1 2 3]));
            response(i,4)=-1;
        end
    end
end
response(togo,3)=1;
tmp=find(goodones==0);
tst=find(response==-1);
response(tmp,:)=-1;
response(tst)=-1;


resp=[];
for i=1:1067
    if ~isempty(find(response(i,:)==1))
        resp=[resp;i];
    end
end


save(fullfile('/media/NERFFS01/Users/Katja/human_retina','h_responding'),'response','resp','pols')




%% TABLE: which cell responds to what stimulus
% 1) flash, 2) LF, 3) DG, 4) chirp, 5) 6vel black, 6) 6vel white, 7) DS

clear
close all
superpath='R:\Users\Katja\human_retina';

Dtype='hpw';

% idlist=[8 9 0 19 16 17 10/11];

d=3
load(fullfile(superpath,[Dtype(d),'_info']))
load(fullfile(superpath,[Dtype(d),'_normPeaks']))
response=zeros(size(info,1),7);
pols=zeros(size(info,1),2);
pols=pols-2;
dates=unique(info(:,1));
did=zeros(size(info,1),1);


for i=1:size(info,1)
    found=0; m=0;
    while found==0
        m=m+1;
        if ~isempty(regexp(dates{m},info{i,1}))
            found=1;
        end
    end
    did(i)=m;
    
    
    if ~isempty(info{i,8}) %flash
        if info{i,8}~=0 && ~isnan(info{i,8})
            response(i,1)=1;
            if info{i,8}==1;
                pols(i,1)=1;
            elseif info{i,8}==-1
                pols(i,1)=-1;
            elseif info{i,8}==2
                pols(i,1)=0;
            end
        end
    end
    if ~isempty(info{i,9}) %LF
        if info{i,9}~=0 && ~isnan(info{i,9})
            response(i,2)=1;
            if info{i,9}==1;
                pols(i,2)=1;
            elseif info{i,9}==-1
                pols(i,2)=-1;
            end
        end
    end
    if ~isempty(info{i,19})
        if info{i,19}~=0 && ~isnan(info{i,19})
            response(i,4)=1;
        end
    end
    if ~isempty(info{i,16})
        if info{i,16}~=0 && ~isnan(info{i,16})
            response(i,5)=1;
        end
    end
    if ~isempty(info{i,17})
        if info{i,17}~=0 && ~isnan(info{i,17})
            response(i,6)=1;
        end
    end
    
    %         if ~isempty(info{i,10}) ||  ~isempty(info{i,11})
    %             if info{i,10}~=0 || info{i,11}~=0
    %                 if ~isnan(info{i,10}) || ~isnan(info{i,11})
    %                     response(i,7)=1;
    %                 end
    %             end
    %         end
    if ~isempty(info{i,12})
        if info{i,12}~=0 && ~isnan(info{i,12})
            response(i,7)=1;
        end
    end
    
    if d==2
        if ~isempty(find(did(i)==[1 2 3 7]));
            response(i,4)=-1;
        end
        if ~isempty(find(did(i)==7));
            response(i,3)=-1;
        end
        if ~isempty(find(did(i)==9));
            response(i,[3 5 6])=-1;
        end
    end
    
    if d==3
        if ~isempty(find(did(i)==[1 2 3]));
            response(i,4)=-1;
        end
    end
end
response(togo,3)=1;
tmp=find(goodones==0);
tst=find(response==-1);
response(tmp,:)=-1;
response(tst)=-1;


resp=[];
for i=1:231
    if ~isempty(find(response(i,:)==1))
        resp=[resp;i];
    end
end


save(fullfile('R:\Users\Katja\human_retina','w_responding'),'response','resp','pols')