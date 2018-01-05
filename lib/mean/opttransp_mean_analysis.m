% Author: Emmanuel K. Kalunga
function [Snp, sum_dis,tr_Sn,niter, crit, Sn] = opttransp_mean_analysis(B,args)

K = size(B,3); % Number of matrices
epsilon = 1e-10; % stoping criteria threshold
N_itermax = 200;
niter = 0;
Sn = eye(size(B(:,:,1)));% Initialize with S0 = I (identity matrix)
A = Sn;
sum_dis(niter+1) = 0;
for m = 1:size(B,3)
    sum_dis(niter+1) = sum_dis(niter+1) + (distance_opttransp(A,B(:,:,m)))^2;
end
while (niter<N_itermax)
    niter = niter+1;
    sum = zeros(size(B(:,:,1)));
    for m = 1:K
        sum = sum + (1/K)*(Sn^(1/2)*B(:,:,m)*Sn^(1/2))^(1/2);
    end
    Snp = Sn^(-1/2) * sum^2 * Sn^(-1/2);
    A = Snp;
    sum_dis(niter+1) = 0;
    for m = 1:size(B,3)
        sum_dis(niter+1) = sum_dis(niter+1) + (distance_opttransp(A,B(:,:,m)))^2;
    end
    tr_Sn(niter) = trace(Sn);
    crit = distance_wasserstein(Sn, Snp);
    if crit < epsilon
        %A = Snp;
        break;
    end
    Sn = Snp;
end

if niter==N_itermax
    disp('Warning : The maximum number of iterations was reached');
end


