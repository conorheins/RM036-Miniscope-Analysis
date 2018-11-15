function img = circular_constraints(img,ctr_pixel,show_imgs)
% enforce a stronger prior to neuron's spatial component by assuming that
% the gradient at each pixel is pointing to the peak location
% It also removes isolated pixels.
if nargin<3
    show_imgs = false;
end
[tmp1, tmp2, ~] = find(img);
if isempty(tmp1)
    return;
end
rmin = min(tmp1);
rmax = max(tmp1);
cmin = min(tmp2);
cmax = max(tmp2);
[nr, nc] = size(img);   % image dimension
if (rmax-rmin<1) || (cmax-cmin<1)
    return;
end
% crop a small region
if rmin==1 && rmax==nr && cmin==1 && cmax==nc
    
    if show_imgs
        figure;
        subplot(121);
        imagesc(img);
        axis equal off tight;
    end
    
    [y0, x0] = ind2sub([nr, nc], ctr_pixel);
    [x, y] = meshgrid(1:nc, 1:nr);
    [fx, fy] = gradient(img);
    ind = ((fx.*(x0-x)+fy.*(y0-y)) < 0) & (img<img(ctr_pixel)/3);
    img(ind) = 0;
    
    % remove isolated pixels
    l = bwlabel(img, 4);
    ind = imdilate(l==l(ctr_pixel), strel('square', 3));
    img(~ind) = 0;
    
    if show_imgs
        subplot(122);
        imagesc(img);
        axis equal off tight;
        pause;
        close;
    end
else
    [y0, x0] = ind2sub([nr, nc], ctr_pixel);
    relative_r = y0 - rmin + 1;
    relative_c = x0 - cmin + 1;
    if relative_r < 1 || relative_c < 1 || relative_r > length(rmin:rmax) || relative_c > length(cmin:cmax) % in case original 'ctr_pixel' is now at a zero-pixel outside of rmin:rmax, cmin:cmax range
        [~,ctr_pixel] = max(reshape(img(rmin:rmax,cmin:cmax),[length(rmin:rmax)*length(cmin:cmax),1]));
    else
        ctr_pixel = sub2ind([length(rmin:rmax),length(cmin:cmax)],relative_r,relative_c);
    end
    tmp_img = circular_constraints(img(rmin:rmax, cmin:cmax), ctr_pixel,show_imgs);
    img(rmin:rmax, cmin:cmax) = tmp_img;
end

