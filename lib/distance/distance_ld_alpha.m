function a = distance_ld_alpha(A,B,alpha)

%a = sqrt(log(det((A+B)/2))-0.5*log(det(A*B)));

% Author: Emmanuel K. Kalunga
%alpha = 0.5;
% alpha = -0.9;
a = (4/(1-alpha^2))*log( (det(((1-alpha)/2)*A + ((1+alpha)/2)*B)) / (det(A)^((1-alpha)/2)*det(B)^((1+alpha)/2)) );
