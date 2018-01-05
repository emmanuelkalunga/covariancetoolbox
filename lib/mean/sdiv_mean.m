% Author: Emmanuel K. Kalunga
% A = sdiv_mean(B,epsilon,tol)
%
% Compute the S-Divergence Centroid.
% A : Centroid of K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le crit�re < tol


function [MS,Snp,niter, crit, Sn] = sdiv_mean(B)

K = size(B,3); % Number of matrices
epsilon = 1e-10; % stoping criteria threshold
N_itermax = 200;
niter = 0;
%Sn = eye(size(B(:,:,1)));% Initialize with S0 = I (identity matrix)
Sn = mean(B,3)^(-1);

% sumB = zeros(size(B(:,:,1)));
sumB = sum(B,3)/K;
while (niter<N_itermax)
    niter = niter+1;    
    summ = zeros(size(B(:,:,1)));
    for m = 1:K
        summ = summ + ((Sn^(-1) + B(:,:,m))/2)^(-1);
    end
    Snp =  (1/K)*summ;
    %crit = det(Snp)-det(Sn)
    crit = distance_sdiv(Snp^(-1),Sn^(-1));
    if crit < epsilon
        %MS = Sn^(-1);
        MS = Snp^(-1);
        break;
    end
    Sn = Snp;
end

if niter==N_itermax
    disp('Warning : The maximum number of iterations was reached');
    MS = Snp^(-1);
end
