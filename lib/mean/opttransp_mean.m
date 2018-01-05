% A = riemann_mean(B,epsilon,tol)
%
% Calcul du barycentre des matrice de covariances.
% A : baricentre des K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le crit�re < tol

%ORIGINAL CODE BY BARACHAT
function [A critere niter] = opttransp_mean(B,args)
Imat = size(B,3);
N_itermax = 200;
if (nargin<2)||(isempty(args))
    tol = 10^-3;
    %A = mean(B,3);
    A = eye(size(B,1));
else
    tol = args{1};
    A = args{2};
end

niter = 0;
fc = 0;
K = A^(0.5);

while (niter<N_itermax)
    niter = niter+1;
   
    Ktmp = 0;
    for i=1:Imat
       Ktmp = Ktmp + (K*B(:,:,i)*K)^(0.5); %Modifiee par Emmanuel 
    end
    Ktmp = (Ktmp)^(0.5);
    
    fcn = norm(Ktmp-K,'fro');
    K = Ktmp;
    % improvement
    conv = abs((fcn-fc)/fc);
    if conv<tol % break if the improvement is below the tolerance
       break; 
    end
    fc = fcn;
end

A = ((1/Imat)*K)^2;

if niter==N_itermax
    disp('Warning : Nombre d''itérations maximum atteint');
end

critere = fc;
% 
% function [A critere niter] = opttransp_mean(B,args)
% Imat = size(B,3);
% N_itermax = 50;
% if (nargin<2)||(isempty(args))
%     tol = 10^-4;
%     A = mean(B,3);
%     %A = eye(size(B,1));
% else
%     tol = args{1};
%     A = args{2};
% end
% 
% niter = 0;
% fcn = intmax('int64');
% K = A^(0.5);
% 
% while (niter<N_itermax)
%     niter = niter+1;
%    
%     Ktmp = 0;
%     for i=1:Imat
%        Ktmp = Ktmp + (1/Imat)*(K*B(:,:,i)*K)^(0.5);
%     end
%     Ktmp = (Ktmp)^(0.5);
%     
%     fcn = norm(Ktmp-K,'fro');
%     K = Ktmp;
%     if fcn<tol % break if the improvement is below the tolerance
%        break; 
%     end    
% end
% 
% %A = ((1/Imat)*K)^2;
% A = K^2;
% if niter==N_itermax
%     disp('Warning : Nombre d''itérations maximum atteint');
% end
% 
% critere = fcn;




% Author: Emmanuel K. Kalunga
% function [Snp,tr_Sn,niter, crit, Sn] = opttransp_mean(B,args)
% 
% K = size(B,3); % Number of matrices
% epsilon = 1e-10; % stoping criteria threshold
% N_itermax = 200;
% niter = 0;
% Sn = eye(size(B(:,:,1)));% Initialize with S0 = I (identity matrix)
% while (niter<N_itermax)
%     niter = niter+1;
%     sum = zeros(size(B(:,:,1)));
%     for m = 1:K
%         sum = sum + (1/K)*(Sn^(1/2)*B(:,:,m)*Sn^(1/2))^(1/2);
%     end
%     Snp = Sn^(-1/2) * sum^2 * Sn^(-1/2);
%     tr_Sn(niter) = trace(Sn);
%     crit = distance_wasserstein(Sn, Snp);
%     if crit < epsilon
%         A = Snp;
%         break;
%     end
%     Sn = Snp;
% end
% 
% if niter==N_itermax
%     disp('Warning : The maximum number of iterations was reached');
% end


