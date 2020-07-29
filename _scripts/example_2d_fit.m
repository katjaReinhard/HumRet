
sf=[0.025 0.05 0.1 0.2 0.4];
tf=[0.5 1  2 4 8 16];

raw=[0 0 0 0 0;  % example matrix of raw response, column are SFs, rows are TFs;
    0 0 1 0 0;
    0 1 1 1 0;
    0 1 1 1 0;
    0 0 1 0 0;
    0 0 0 0 0];
raw=raw+rand(size(raw))*0.2;

[mat_fit,xhat]=fitsftf(raw,sf,tf); % xhat contains the parameters of the fit; mat_fit is the fit

%%
clf;
subplot(2,2,1);
imshow(raw);
title('raw');
subplot(2,2,2);
imshow(mat_fit);
title('fit')

