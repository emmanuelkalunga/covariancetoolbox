function th = threeshold_estimation_z(distance,method, z_score)
    if nargin < 2
        method = 'mean';
    end
    
    switch method
        case 'mean'
            %th = mean(distance)+2.5*std(distance);
            th = mean(distance)+z_score*std(distance);
        case 'median'
            %th = median(distance)+2.5*1.4826*mad(distance,1); 
            th = median(distance)+z_score*1.4826*mad(distance,1);
    end
