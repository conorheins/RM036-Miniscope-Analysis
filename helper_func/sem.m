function [ sem_output ] = sem( x,DIM )
%Standard error of the mean Function to calculate 1 SEM

if ~exist('DIM','var') || nargin < 2
    sem_output = std(x(:))./sqrt(length(x(:)));
else
    sem_output = std(x,0,DIM)./sqrt(size(x,DIM));
end



end

