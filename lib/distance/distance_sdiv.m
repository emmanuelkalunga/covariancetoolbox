%SD = distance_sdiv(P,Q)
%Compute S-divergence between two distributions given by they
%covariance matrices
%Input: P and Q, two covariance matrice
%Output: SD, the S-divergence between P and Q
%Author: Emmanuel K. Kalunga

function SD = distance_sdiv(P,Q)
%m = min([min(min(Q)) min(min(P))]);
%P = P/m;
%Q = Q/m;
if trace(P*Q) < length(diag(P))
    A = ( (sqrt(trace(P*Q)))^(-1) ) * eye(size(P)); % Congruent trasformation matrix for normalization to Q>0 and P>0
    P = A*P*A';
    Q = A*Q*A';
end
SD = abs(log(det((P+Q)/2)) - 0.5*log(det(P*Q)));