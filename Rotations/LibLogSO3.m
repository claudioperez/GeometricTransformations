function vect = LogSO3(R,option,small_option)
%LOGSO3 Inverse of the exponential map on SO(3).
%  vect = LogSO3(R)
%  vect = LogSO3(R,option)
%  Returns the axial parameters associated with the rotation `R`. The result
%  should satisfy the following equality for any 3-vector, `v`:
%
%        LogSO3(expm(spin(v))) == v
%
%  where `expm` is the Matlab built-in matrix exponential, and `spin` is a function
%  which produces the skew-symmetric 3x3 matrix associated with vector `v`.
%
%  Parameters
%    R       (3x3)   Rotation matrix.
%    option  string  Algorithm option; default is 'Quat'. See reference below.
%                    All other options are for internal use.
%
%  Remarks
%
%  - Does not check if input is really a rotation. If you arent sure, use
%    Matlab's `logm`. This is slower, but more general. Note however that`logm`
%    may not be as accurate in corner cases and sometimes returns complex-valued
%    results.
%  - The angle corresponding to the returned vector is always in the interval [0,pi].
%
%
%  References
%  1. Nurlanov Z (2021) Exploring SO(3) logarithmic map: degeneracies and 
%     derivatives.
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

if nargin < 2
  option = 'Quat';
end
% option = 'Quat';

switch option
  case 'Quat'
    vect = Quat2Axis(Rmat2Quat(R));

  case 'Vect'
    trR = trace(R);
    % Treat corner cases
    if  abs(trR - 3) <= 1e-18
      % "R is near the identity, result is near 0"
      angle = 0;
      vect  = [0; 0; 0];
    elseif abs(trR - 3)  <= 1e-6
      if nargin > 2
        vect = LogSO3(R, small_option);
      else
        vect = LogSO3(R, '4');
      end
    elseif abs(trR + 1) <= eps(1.0)
      % "Rotation by odd multiple of Pi."
      [m,idx] = max(diag(R));
      axis = (R(:,idx) + (1:3==idx)') / sqrt(2*(1+m));
      angle = pi;
      vect = axis*pi;
    else
      angle = real(acos(complex(0.5*(trR-1))));
      vect = angle/sin(angle)*askew(R);
    end

  case 'Old'
    trR = trace(R);
    % Treat corner cases
    if  abs(trR - 3) <= eps
      % "R is near the identity, result is near 0"
      angle   = 0;
      vect    = [0; 0; 0];
    elseif abs(trR + 1) <= eps
      % "Rotation by odd multiple of Pi."
      [m,idx] = max(diag(R));
      axis    = (R(:,idx) + (1:3==idx)') / sqrt(2*(1+m));
      angle   = pi;
      vect    = axis*pi;
    else
      angle   = acos(0.5*(trR - 1));
      vect    = angle/sin(angle)*askew(R);
    end

  case '0'
    % Rough approximation
    vect = askew(R);

  case '4'
    % Power series approximation
    trR   = trace(R);
    cs    = max(min(0.5 * (trR - 1), 1), -1);
    sn    = 0.5 * sqrt(max(0, (3 - trR) * (1 + trR)));
    angle = atan2(sn, cs);
    vect  = (1 + angle.^2/6 + 7/360 * (angle.^4))*askew(R);

  case 'B2'
    %
    vect = 2*askew(R);
    ang  = norm(vect);
    if ang == 0
      vect = [0; 0; 0];
    else
      vect = asin(ang/2)/ang*vect;
    end

  case 'C1'

    vect = asin(askew(R));

  case 'pade'
    % Matlab's Padé approximant to the log on GL(3)
    vect = Axial(logm(R));

end % switch

end % function LogSO3

%%  function askew ---------------------------------------------------------------------
function S = askew(A)
%  Axial vector (3x1) of the skew-symmetric part of a 3x3 matrix.
%  Mathematically, this can be expressed as
%      w = 0.5 e[A - A']
%        = 0.5*axial(A - A');
%  where `e` is the alternating tensor.

  S = 0.5*[ A(3,2) - A(2,3)
            A(1,3) - A(3,1)
            A(2,1) - A(1,2) ];

end

