function A = get_A(Q,labels)
cl = unique(labels);
N_cl = length(cl);

N = size(Q,3);
D = size(Q,1);
index = reshape(triu(ones(D)),D*D,1)==1;
S = zeros(D*(D+1)/2,N);
SS = zeros(D*(D+1)/2);
SD = zeros(D*(D+1)/2);
for i=1:N
    tmp = reshape(sqrt(2)*triu(Q(:,:,i),1)+diag(diag(Q(:,:,i))),D*D,1);
    S(:,i) = tmp(index);
end

i = 1;
for cl_i = labels
	j = 1;
	for cl_j = labels
		if cl_i == cl_j
			SS = SS+(S(:,i)-S(:,j))*(S(:,i)-S(:,j))';
		else
			SD = SD+(S(:,i)-S(:,j))*(S(:,i)-S(:,j))';
		end
		j = j+1;
	end
	i = i+1;
end

A = SS^(-0.5) * (SS^(0.5) * SD * SS^(0.5))^(0.5) * SS^(-0.5);

