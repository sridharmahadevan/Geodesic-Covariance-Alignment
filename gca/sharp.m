function [G,mA,mB]=sharp(A,B,t)

% [G,ma,mB]=SHARP(A,B,t) computes the point G=g(t) of the geodesic g joining 
%  A and B such that g(0)=A and g(1)=B, i.e., A(A^{-1}B)^t=B(B^{-1}A)^{1-t}
%
% A,B: positive definite matrices
% G: the point on the geodesic
% t: real parameter
% ma,mb: condition numbers of thematrices A and B, respectively
%
% References
% [1] B. Iannazzo, The geometric mean of two matrices from a computational
% viewpoint, arXiv preprint arXiv:1201.0101, 2011.

mA=cond(A);mB=cond(B);
if (mA>mB) % swap A and B if B is better conditioned
  C=A;A=B;B=C;t=1-t;
end

RA=chol(A);RB=chol(B);
Z=RB/RA;
[U V]=eig(Z'*Z);
T=diag(diag(V).^(t/2))*U'*RA;
G=T'*T;
