% Desccription:
%     W_MEAN computes the weighte mean of scalars or arrays
%     W_MEAN(X) retuns the normal (non-weighted) mean
%     W_MEAN(X,DIM) retuns the normal (non-weighted) mean along dimension DIM
%     W_MEAN(X,DIM,W) returns the weighted mean with wight W
% Author:
%     Emmanuel K. Kalunga

function M = w_mean(X, varargin)

if isempty(varargin)
    M = mean(X); % Ordinary mean or average
    return;
elseif numel(varargin)==1
    dim = varargin{1};
    M = mean(X,dim); % Ordinary mean with specified dimension
    return;
elseif numel(varargin)==2
    dim = varargin{1};
    W = varargin{2};
    W = W(:); %making sure W is a column vector
    if dim == 1
        M = X'*W;
        M = M';
        return;
    elseif dim == 2
        M = X*W;
        return;
    else
        error('myApp:argChk','The mean can only be computed along dimension 1 or 2');
    end
end

end
    
