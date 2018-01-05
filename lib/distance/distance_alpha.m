% [] = XXXX()
%
% -----------------------------Definition---------------------------------%
%
% description
%
% usage : 
%
% -----------------------------Input--------------------------------------%
%
% XX : description + dimension
%
% -----------------------------Output-------------------------------------%
%
% XX : description + dimension
%
% -----------------------------References---------------------------------%
%
% [1] : XXX
%
%
%   Project : BCI-EEG
%
%   author : A. Barachant /Emmanuel K. Kalunga
%   date : 2011-XXXX
%   version : 1.0 
%   status : a terminer, termin�   
%   CEA/GRENOBLE-LETI/DTBS
%
%   See also distance_riemann, distance_kullback, distance_logeuclid, norm.

% [EOF: XXX.m]

function d = distance_alpha(C1,C2,method_dist, alpha)

if (nargin<3)||(isempty(method_dist))
    method_dist = 'euclid';
    alpha = 0;
end

switch method_dist
    case 'riemann'
        d = distance_riemann(C1,C2);
    case 'kullback'
        d = distance_kullback(C1,C2);
    case 'jeffreys'
	d = distance_jeffreys(C1,C2);
    case 'logeuclid'
        d = distance_logeuclid(C1,C2);
    case 'opttransp'
        d = distance_opttransp(C1,C2);
    case 'ld'
        d = distance_ld_alpha(C1,C2,alpha);
    case 'wasserstein'
        d = distance_wasserstein(C1,C2); %EK
    case 'sdivergence'
        d = distance_sdiv(C1,C2); %EK
    case 'harmonic'
        d = distance_harm(C1,C2); %EK
    otherwise
        d = sqrt(norm(C1-C2,'fro'));
end
