function Cscale = getCscale(data,target)
if nargin==1
    target = 0.98;
end

if target==1
    Cscale = [min(data(:)) max(data(:))];
    return
end

x0 = median(data(:));
Cscale(1) = fminsearch(@(x) errorFunc(x,data(:),target),x0);

targetHigh = 1-target;
Cscale(2) = fminsearch(@(x) errorFunc(x,data(:),targetHigh),x0);


function error = errorFunc(threshold,data,target)
error = sum(target - sum(data>threshold)/numel(data))^2;