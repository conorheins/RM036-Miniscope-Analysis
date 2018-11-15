function [Atrim,ind_small] = trimSpatial_noObj(A,Ysiz,thr,sz,minpixel)
% spatially trim components with image opening and thresholding

        % remove small nonzero pixels
        
        Atrim = zeros(size(A));
    
        if nargin<3;    thr = 0.02; end;
        if nargin<4;    sz = 2; end;
        if nargin<5;    minpixel = 5; end
            
        se = strel('disk', sz);
        ind_small = false(size(A,2),1);
        d1 = Ysiz(1); d2 = Ysiz(2);
            
        for m=1:size(A,2)
            ai = A(:,m);
            ai_open = imopen(reshape(ai,d1,d2), se);
            
            temp = ai_open>max(ai)*thr;
            l = bwlabel(temp, 8);   % remove disconnected components
            [~, ind_max] = max(ai_open(:));
            
            ai(l(:)~=l(ind_max)) = 0;
            ai(ai<max(ai)*0.3)=0;
           
            if sum(ai(:)>0) < minpixel %the ROI is too small
                ind_small(m) = true;
            end
            
            Atrim(:, m) = ai(:);
            
        end
        
        ind_small = find(ind_small);

end

