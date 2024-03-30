function R = ExpSO3(th) %, other)
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

% Force input to be a column vector
th = reshape(th, [3 1]);
% Form the first Gib coefficients
[a0, a1, a2] = GibSO3(th);
% Form 3x3 skew-symmetric matrix Th from axial vector th
Th = Spin(th);
%
R = eye(3) + a1*Th + a2*Th*Th;

