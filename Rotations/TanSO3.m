function t = TanSO3(rot)
%
%   Compute right differential of the exponential.
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

%   Compute norm of rot
    angle2 = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3);

%   Check for angle near 0
    if angle2 < 1e-08
      fac1  = 1.0     - angle2*(1.0/6.0   ...
                      - angle2*(1.0/120.0 ...
                      - angle2/5040.0));
      fac2  = 0.5     - angle2*(1.0/24.0  ...
                      - angle2*(1.0/720.0 ...
                      - angle2/40320.0));
      fac3  = 1.0/6.0 - angle2/(1.0/120.0 ...
                      - angle2/(1.0/5040.0 ...
                      - angle2/362880.0));
    else
      angle = sqrt (angle2);
      fac1  = sin(angle)        /angle;   % a1
      fac2  = (1.0 - cos(angle))/angle2;  % a2
      fac3  = (1.0 - fac1      )/angle2;  % a3
    end

%   Assemble tangent
    t = zeros(3,3);
    t(1,1)  =         fac1 + fac3*rot(1)*rot(1);
    t(1,2)  = -rot(3)*fac2 + fac3*rot(1)*rot(2);
    t(1,3)  =  rot(2)*fac2 + fac3*rot(1)*rot(3);
    t(2,1)  =  rot(3)*fac2 + fac3*rot(2)*rot(1);
    t(2,2)  =         fac1 + fac3*rot(2)*rot(2);
    t(2,3)  = -rot(1)*fac2 + fac3*rot(2)*rot(3);
    t(3,1)  = -rot(2)*fac2 + fac3*rot(3)*rot(1);
    t(3,2)  =  rot(1)*fac2 + fac3*rot(3)*rot(2);
    t(3,3)  =         fac1 + fac3*rot(3)*rot(3);

end

