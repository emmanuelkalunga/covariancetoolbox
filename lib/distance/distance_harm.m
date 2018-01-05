%dW = distance_wasserstein(P,Q)
%Compute the 2-Wasserstein disctance between two distributions given by they
%covariance matrices
%Input: P and Q, two covariance matrice
%Output: dW, the Wasserstein-2 distance between P and Q
%Author: Emmanuel K. Kalunga

function HD = distance_harm(P,Q)
%m = min([min(min(Q)) min(min(P))]);
%P = P/m;
%Q = Q/m;
% if trace(P*Q) < length(diag(P))
%     A = ( (sqrt(trace(P*Q)))^(-1) ) * eye(size(P)); % Congruent trasformation matrix for normalization to Q>0 and P>0
%     P = A*P*A';
%     Q = A*Q*A';
% end
HD = norm((P^(-1)-Q^(-1)),'fro');