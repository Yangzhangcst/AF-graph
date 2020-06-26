function  histogram = smoothKDE(histogram,KDthreshold,smoothing)
% smooth features by the KDE method
%     cnormalize
    histogram1 = [];
    for j=1:size(histogram,2)
        histogram2 = histogram(:,j); 
        histogram2 = histogram_smoothing(histogram2,KDthreshold,smoothing);
        histogram1 = [histogram1,histogram2(:,2)];
    end
    histogram = histogram1';
    histogram(:,1) = mean(histogram,2);
    [m,~] = size(histogram);
    d = sqrt( sum(histogram.^2,1) );
    d(d<eps)=1;
    histogram = histogram ./ repmat( d, m,1 );