% KLS = Kullback_Leibler_sym(P,Q)
% Kullback_leibler distance.

function KL = distance_kullback(P,Q)

%KL = sqrt( 0.5* trace( P/Q + P\Q - 2*eye(size(P)))); %Barachant Implementation of symmetric KL
%KL = 0.5*(trace(P\Q) - 2*eye(size(P)) ); %-- Formula from Kang2009 and wikipedia (Kullback–Leibler divergence for multivariate normal distributions)
%KL = 0.5*( trace(Q\P) - log(det(P)/det(Q)) - size(P,1) );
KL = 0.5*( trace(Q\P) - log(det(Q\P)) - size(P,1) ); %This is equal to the previous expression
%KL = trace(Q\P) - 0.5*(log(det(P)/det(Q))) - size(P,1) ;
%KL = sqrt(0.5*( trace(Q\P) - log(det(P)/det(Q)) - size(P,1) ) + 0.5*( trace(P\Q) - log(det(Q)/det(P)) - size(P,1) ))

