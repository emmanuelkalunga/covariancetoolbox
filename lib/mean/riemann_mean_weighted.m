% A = riemann_mean(B,epsilon,tol)
%
% Calcul du barycentre des matrice de covariances.
% A : baricentre des K matrices NxN
%
% B : Matrice NxNxK
% epsilon : Pas de la descente de gradient
% tol : arret de la descente si le crit�re < tol


function [A, critere, niter, bad_cond] = riemann_mean_weighted(B,W,args)

N_itermax = 200; %changed from 100 to 200 by E. Kalunga
if (nargin<3)||(isempty(args))
    tol = 10^-5;
    if size(B,3)>1  %Added by E. Kalunga
        A = mean(B,3);
    else            %Added by E. Kalunga
        A = B;      %Added by E. Kalunga
    end
else
    tol = args{1};
    A = args{2};
end

niter = 0;
fc = 0;
bad_cond = 0;
while (niter<N_itermax)
    niter = niter+1;
    % Tangent space mapping
    T = Tangent_space(B,A);
    % sum of the squared distance
    fcn = sum(sum(T.^2));
    % improvement
    conv = abs((fcn-fc)/fc);
    if conv<tol % break if the improvement is below the tolerance
       break; 
    end
    % arithmetic mean in tangent space
    %TA = mean(T,2);
    TA = w_mean(T,2,W); % function for weighted mean writen by E. Kalunga                  <===================== W is the weight obtained from Kullback divergence
    % back to the manifold
    A_1 = A;
    A = UnTangent_space(TA,A);
    %Stop when there is NaN and/or Inf in matrices in B %Added by E. Kalunga
    %**************************************************************
    if (sum(sum(isnan(A))) || sum(sum(isinf(A)))) 
        fprintf('Matrice not well conditioned when returned from tangent space\n');
        bad_cond = 1;
        A = A_1;
        break;
    end   
    
    fc = fcn;
end

if niter==N_itermax
    disp('Warning : Nombre d''iterations maximum atteint');
end

critere = fc;
