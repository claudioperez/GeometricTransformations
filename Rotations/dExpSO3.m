function T = dExpSO3(th, dth)
%
%  [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric 
%      Transformations of Finite Elements" 
%      Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------


% Force input to be a column vector
th  = reshape(th,  [3 1]);
% Form first Gib coefficients
[~,a1, a2, a3] = GibSO3(th);
% Form skew-symmetric matrix from 3-vector th
Th = Spin(th);
%
if nargin == 1
% T = eye(3) + a2*Th + a3*Th*Th;
  T = a1*eye(3) + a2*Th + a3*th*th';
  return
else
% T = dth + a2*Th*dth + a3*Th*Th*dth;
  T = a1*dth + a2*cross(th, dth) + a3*th*(th'*dth);
  return
end

