function [ psf ] = create_filter( gSig,gSiz )
% Creates a gaussian filter with standard deviation gSig and total size
% gSiz, using 'fspecial' matlab function. Finishes by subtracting the mean
% and pushing low values to zero

psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;

end

