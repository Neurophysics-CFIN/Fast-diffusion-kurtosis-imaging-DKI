close all
clearvars
addpath('..')



% load data example (single slice 496 b-matrices 96x96x1x496)
load('data.mat','data','bmats')

% create appropriate mask
mask = data(:,:,:,1)>10;



% make fit
[DKIcylPars,sse] = fitDKIcyl(data,bmats,mask);



% calculate metrics
metrics = DKIcylMetrics(DKIcylPars);

% make some plots
titles = {'Dpar','Dper','Dbar','Wpar','Wper','Wbar','Kpar','Kper','Kbar'};
cRange = 0.9;
for n = 1:size(metrics,4)
    map = metrics(:,:,1,n);
    
    figure
    imagesc(map,getCscale(map(mask),cRange))
    title(titles{n})
end


