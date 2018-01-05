%Description:
%   Compute the Mahalanobis distance between the matrix P and the set
%   (group) of matrices Q
%
%Parameters:
%   P: a square matrix of dimension D
%   Q: a set of N square matrices of dimemsion D belonging to one class
%
%Return:
%   Return a real number (R^+), the Mahalanobis distance between P and Q

function MD = distance_mahalanobis(P,Q,COVMAT,labels)
N = size(Q,3);
D = size(Q,1);
S_Q = zeros(D*(D+1)/2,N);
S_P = zeros(D*(D+1)/2);

index = reshape(triu(ones(D)),D*D,1)==1;
for i=1:N
    tmp = reshape(sqrt(2)*triu(Q(:,:,i),1)+diag(diag(Q(:,:,i))),D*D,1);
    S_Q(:,i) = tmp(index);
end
tmp = reshape(sqrt(2)*triu(P,1)+diag(diag(P)),D*D,1);
S_P = tmp(index);

S1 = S_P;
S2 = mean(S_Q,2);
% MQ = mean(Q,3);
%MQ = riemann_mean(Q); %mean(Q,3);
%tmp = reshape(sqrt(2)*triu(MQ,1)+diag(diag(MQ)),D*D,1);
%S2 = tmp(index);
%C = shcovft(S_Q');
%C = (1/N)*(bsxfun(@minus,S_Q,S2))*(bsxfun(@minus,S_Q,S2))';
C = get_A(COVMAT, labels);
MD = sqrt((S1-S2)'*inv(C)*(S1-S2));

