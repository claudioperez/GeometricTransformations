function Xi = ddTanSO3(theta, a, b)
%
%  [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric 
%      Transformations of Finite Elements" 
%      Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------


% Xi_{D^2T}
theta = reshape(theta, [3 1]);
a     = reshape(a    , [3 1]);
b     = reshape(b    , [3 1]);

[a0, a1, a2, a3, b1, b2, b3, c1, c2, c3] = GibSO3(theta);

I   = eye(3);
Xi = a3*(a*b' + b*a') + b1*(a'*b)*I ...
   + b2*(cross(a,b)*theta' + theta*cross(a,b)' +(cross(theta, a)'*b)*I) ...
   + b3*((theta'*a)*(b*theta' + theta*b') + (theta'*b)*(a*theta' + theta*a') + (theta'*a)*(theta'*b)*I) ...
   + (c1*(a'*b) + c2*(cross(theta,a)'*b) + c3*(theta'*a)*(theta'*b))*theta*theta';

