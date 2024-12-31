function Xi = dTanSO3(th, a, repr)
%  repr     'L' or 'R' indicating left or right representation, 
%           respectively, for the tangent space of SO(3)
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

if nargin < 3
  repr = 'L';
end

th = reshape(th, [3 1]);
a  = reshape(a , [3 1]);


[~, a1, a2, a3, b1, b2, b3] = GibSO3(th);

I   = eye(3);
switch repr
  case 'R'
    Xi = - a2*Spin(a) + a3*(th'*a)*I + a3*th*a' ...
         + b1*a*th' + b2*cross(th,a)*th' + b3*th'*a*(th*th');
  case 'L'
    Xi =   a2*Spin(a) + a3*(th'*a)*I + a3*th*a' ...
         + b1*a*th' - b2*cross(th,a)*th' + b3*th'*a*(th*th');
end

