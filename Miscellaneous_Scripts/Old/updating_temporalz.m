
figure;
for i = 1:size(Atrim,2)
    
    subplot(121);
    [~,center] = max(Atrim(:,i));
    [x,y] = ind2sub([d1,d2],center);
    mask = reshape(Atrim(:,i),d1,d2);
    imagesc(mask((x-10):(x+10),(y-10):(y+10)));
    
    subplot(122);
    plot(C(i,:));
    
    % update spatial/temporal components of the merged neuron
    temp = Atrim(:, i)*C(i,:);
    
    ci = C(i,:);
    for miter=1:10
        
        ai = temp*ci'/(ci*ci');
        ci = ai'*temp/(ai'*ai);
        
        subplot(121);
        mask = reshape(ai,d1,d2);
        imagesc(mask((x-10):(x+10),(y-10):(y+10)));
        
        subplot(122);
        plot(ci);
        
        pause;
        
    end
    
end
    