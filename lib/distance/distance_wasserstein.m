%dW = distance_wasserstein(P,Q)
%Compute the 2-Wasserstein disctance between two distributions given by they
%covariance matrices
%Input: P and Q, two covariance matrice
%Output: dW, the Wasserstein-2 distance between P and Q
%Author: Emmanuel K. Kalunga

function dW = distance_wasserstein(P,Q)
%m = min([min(min(Q)) min(min(P))]);
%P = P/m;
%Q = Q/m;
dW_sq = trace(P + Q - 2*(P^(0.5)*Q*P^(0.5))^(0.5) );
dW = sqrt(dW_sq);
