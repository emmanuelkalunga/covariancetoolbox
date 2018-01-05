% from Means of Hermitian positive-definite matrices based on the log-determinant α-divergence function
% Zeineb Chebbi , Maher Moakher.


function [A curve] = logdet_mean_alpha(B,alpha)

K = size(B,3); % Nombre de matrices
A = mean(B,3);
% tol = 10^-3;
tol = 10^-10;
imp = tol + 1;
cnt = 0; %Emmanuel K. Kalunga
maxiter=100;
while (imp>tol) && (cnt < maxiter)
    cnt = cnt+1;
    fc = zeros(size(B,1));

    for i=1:K
        %fc = fc + inv(0.5*B(:,:,i) + 0.5*A);
        fc = fc + inv(((1-alpha)/2)*B(:,:,i) + ((1+alpha)/2)*A);
    end

    Anew = inv(fc/K);
    imp = distance_ld_alpha(Anew,A,alpha);
    A = Anew;
    curve(cnt) = imp;
end
% imp
% cnt

function B = mypseudoinv(A)
    [U S] = svd(A);
    k=1;
    K = size(A,1);
    B = zeros(K,K);
   while (k<=K) && (S(k,k)>10^-15)
        B = B+U(:,k)*(1/S(k,k))*U(:,k)';
        k=k+1;
    end
     

