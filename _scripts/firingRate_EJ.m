stimLength=6000; %in ms

%collect data
allRates=[];
for f=1:length(file)
    spikes=unit{1,2}{file(f),2};
    spikes=round(spikes*10); %0.1 ms
    
    rate=zeros(1,stimLength*10);
    rate(spikes)=1;
    rate(spikes(diff(spikes)==0))=2;
    
    allRates=[allRates;rate];
end

%Gauss filter
sigma = 20; %sigma of 2 ms
flength = 5 * sigma;
x = linspace(-flength / 2, flength / 2, flength);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize

%Convolve
yfilt = filter (gaussFilter,1, mean(allRates));
figure;plot(0.1:0.1:length(yfilt)/10,yfilt*10000)