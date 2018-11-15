function [ component_clean ] = cleanup_footprints(component,ctr_pixel,options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

medfilt_param = options.medfilt_param;
nrgthr = options.nrgthr;
close_elem = options.close_elem;

component_clean = zeros(size(component));
[d1,d2] = size(component);
d = d1*d2;

indf = 0;
valf = 0;

component = medfilt2(component,options.medfilt_param);

component = component(:);

[temp,ind] = sort(component(:).^2,'ascend');
temp = cumsum(temp);

ff = find(temp > (1-nrgthr)*temp(end),1,'first');

BW = zeros(d1,d2);
BW(ind(ff:d)) = 1;
BW = imclose(BW,close_elem);

[L,NUM] = bwlabel(BW);

if NUM > 0
    nrg = zeros(NUM,1);  
    for l = 1:NUM
        ff = (L == l);
        if ~isempty(intersect(ctr_pixel,find(ff)))
            indf = find(ff); break;
        else
            nrg(l) = sum(component(ff).^2);
        end
    end
    if indf == 0
        [~,indm] = max(nrg);
        indf = find(L==indm);
    end 
    valf = component(indf);
else
    valf = 0;
end

component_clean(indf) = valf;

end

