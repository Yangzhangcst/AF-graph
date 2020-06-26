function [H] = histogram_smoothing(H, KDthreshold, smoothing)

    %Remove the first value of the histogram -- the black areas of the
    %image outside the subject
    H(1,:) = 0;
    idx = 1:length(H);
    idx = idx';
    %Remove noise from histogram based on kernel density probability
    [bandwidth,density,X,Y]=kde2d_m([idx H]);

    %Find closest x value on the density image to the values on the histogram
    C=bsxfun(@minus,X(1,:),idx);
    [~,index_X]=min(abs(C),[],2);

    %Find the closest y value on the density image to the values on the histogram
    C = bsxfun(@minus,Y(:,1)',H);
    [~,index_Y]=min(abs(C),[],2);

    %Find the density of the data on the histogram using the X and Y indexs and
    %the density matrix
    for i = 1:length(index_X)-1
       density_histogram(i) = density(index_Y(i), index_X(i));    
    end

    %Calculate the density cutoffs using the KDthreshold and the maximum density 
    density_cutoff = KDthreshold/100*(max(max(density_histogram(:,2))));


    %Remove the data points with a density less than the KDthreshold 
    %(assumed most likely to be noise)
    locations = find(density_histogram > density_cutoff);
    H1 = interp1(locations, H(locations), 1:length(H), 'pchip');
    H  = smoothts(H1, 'e', smoothing); 
    H = [idx H'];
    
end    