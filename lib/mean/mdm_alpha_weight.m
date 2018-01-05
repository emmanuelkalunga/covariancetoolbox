%similar to Barachant's. But takes alpha as an input parameters for alpha
%divervence distance and mean
%Author Emmanuel K. Kalunga

function [Ytest d C] = mdm_alpha_weight(COVtest,COVtrain,Ytrain,varargin)
    
    if isempty(varargin)
        method_mean = 'riemann';
        method_dist = 'riemann';
    else
        method_mean = varargin{1};
        method_dist = varargin{2};
        alpha = varargin{3};
    end
    
    if numel(varargin) == 4
        W = varargin{4};
    else
        W = 1;
    end
    
    labels = unique(Ytrain);
    Nclass = length(labels);
    C = cell(Nclass,1);
    
    % estimation of center
    for i=1:Nclass
        C{i} = mean_covariances_alpha_weight(COVtrain(:,:,Ytrain==labels(i)),method_mean,alpha,W(Ytrain==labels(i)));
    end

    % classification
    NTesttrial = size(COVtest,3);
    
    d = zeros(NTesttrial,Nclass);
    for j=1:NTesttrial
        for i=1:Nclass
            d(j,i) = distance_alpha(COVtest(:,:,j),C{i},method_dist,alpha);
        end
    end
    
    [~,ix] = min(d,[],2);
    Ytest = labels(ix);