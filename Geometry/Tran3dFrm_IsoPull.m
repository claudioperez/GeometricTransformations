function [Tr, AuxMat] = Tran3dFrm_IsoPull(T, e1, option, T0, incr)
%  Computation of a mean triad that is aligned with vector e1.
%  
%  Parameters
%    T  (3,3,n)  3D array containing orthogonal matrices to be averaged.
%    e1 (3,1)    Unit vector to align the result with.
%
%  =========================================================================================
%  function by Claudio Perez                                                         02-2023
%  -----------------------------------------------------------------------------------------
%
e1 = e1(:);       % ensure column vector
Ti = T(:,:,  1);  % Node i rotation relative to previous corot. frame
Tj = T(:,:,end);  % Node j rotation "
switch option
  case 'None'
    AuxMat = [];
    Tr = T0;

  case {'2D'}
    Tr     = [e1  [-e1(2); e1(1)]];
    AuxMat = [];

  case {'SFIN'}
    c      = 0.5;
    nn     = size(T,3);
    i      = floor(0.5*(nn + 1));
    j      = floor(0.5*(nn + 2));
    Ti     = T(:,:,i);
    Tj     = T(:,:,j);
    Phi    = LogSO3(Ti'*Tj);
    Tr     = Ti*ExpSO3(Phi*c);
    AuxMat = [Phi, -Phi];

  case {'R1', 'R2'}
    % Original method of Rankin and Nour-Omid
    e2_tr  = Ti(:,2);
    e3     = cross(e1,e2_tr);
    e3     = e3/norm(e3);
    e2     = cross(e3,e1);
    Tr     = [e1 e2 e3];
    AuxMat = Ti;

  case {'E1', 'E2'}
    AuxMat = ExpSO3(0.5*(incr(4:6,1) + incr(4:6,2)))*T0;
    e2_tr  = AuxMat(:,2);
    e3     = cross(e1,e2_tr);
    e3     = e3/norm(e3);
    e2     = cross(e3,e1);
    Tr     = [e1 e2 e3];
    AuxMat = [e2_tr, Ti(:,2), Tj(:,2)];

  case {'B1', 'B2'}
    % Pacoste method of determining y-axis for 2-node element,
    % implementation by Veronique LeCorvec
    e2_tr  = 0.5*(Ti(:,2) + Tj(:,2));
    e3     = cross(e1,e2_tr);
    e3     = e3/norm(e3);
    e2     = cross(e3,e1);
    Tr     = [e1 e2 e3];
    AuxMat = [e2_tr, Ti(:,2), Tj(:,2)];

  case {'C1', 'C2', 'C3', 'C4'}
    gammaw = Quat2Gibb(Rmat2Quat(Tj*Ti'));

    AuxMat = Gibb2Rmat(gammaw/2)*Ti;

    % Rotate the mean rotation matrix AuxMat 
    % on to e1 to obtain e2 and e3
    r1 = AuxMat(:,1);
    r2 = AuxMat(:,2);
    r3 = AuxMat(:,3);

    if strcmp(option, 'C1') ||  strcmp(option, 'C3')
      %  use the 'mid-point' procedure (Approximately orthogonal)
      e2 = r2 - (r2'*e1)*(e1 + r1)/2;
      e3 = r3 - (r3'*e1)*(e1 + r1)/2;
    else
      e2 = r2 - (r2'*e1)*(e1 + r1)/(1 + e1'*r1);
      e3 = r3 - (r3'*e1)*(e1 + r1)/(1 + e1'*r1);
    end

    Tr = [e1 e2 e3];


  case {'mean'}
    Tr = rotmat(meanrot([rot2quat(Ti); rot2quat(Ti)]));
    AuxMat = [];

  case {'slerp', 'lerp', 'nlerp'}
    % Spherical linear interpolation (SLERP) using quaternions.
    AuxMat = quat2rotm(quatinterp(rotm2quat(Ti), rotm2quat(Tj), 0.5, option));

    % Rotate r1 into e1
    r1 = AuxMat(1:3,1);
    r2 = AuxMat(1:3,2);
    r3 = AuxMat(1:3,3);
    e2 = r2 - (r2' * e1)*(e1 + r1)/(1+e1'*r1);
    e3 = r3 - (r3' * e1)*(e1 + r1)/(1+e1'*r1);
    Tr = [e1 e2 e3];

  case 'KM2'
    % Karcher mean via modified Weiszfeld algorithm
    % - First interpolate Ti and Tj as quaternions
    % - Rotate result into e1
    % - Perform corrective iterations about e1
    tol = 1e-10;
    maxiter = 15;

%   AuxMat = quat2rotm(quatinterp(rotm2quat(Ti), rotm2quat(Tj), 0.5, 'lerp'));
%   AuxMat = Ti;
%   r1 = AuxMat(1:3,1);
%   r2 = AuxMat(1:3,2);
%   r3 = AuxMat(1:3,3);
%   e2 = r2 - (r2'*e1)*(e1 + r1)/(1+e1'*r1);
%   e3 = r3 - (r3'*e1)*(e1 + r1)/(1+e1'*r1);
%   Tr = [e1 e2 e3];

    e2_tr = Ti(:,2);
    e3 = cross(e1,e2_tr);
    e3 = e3/norm(e3);
    e2 = cross(e3,e1);
    Tr  = [e1 e2 e3];
    AuxMat = Ti;

    % Corrective iterations
    r = e1*0.5*(e1'*(LogSO3(Tr'*Ti) + LogSO3(Tr'*Tj)));

    for i = 0:maxiter
      if norm(r) < tol, break; end
      Tr = ExpSO3(r)*Tr;
      r = e1*(0.5*e1'*(LogSO3(Tr'*Ti) + LogSO3(Tr'*Tj)));
    end

  case 'KM'
      % Karcher mean via Weiszfeld algorithm
      tol = 1e-13;

      AuxMat = eye(3);

      r = 0.5*(logm(AuxMat'*Ti) + logm(AuxMat'*Tj));

      while norm(Axial(r)) > tol
        AuxMat = expm(r)*AuxMat;
        r = 0.5*(logm(AuxMat'*Ti) + logm(AuxMat'*Tj));
      end

      r1 = AuxMat(:,1);
      r2 = AuxMat(:,2);
      r3 = AuxMat(:,3);
      e2 = r2 - (r2'*e1)*(e1 + r1)/(1 + e1'*r1);
      e3 = r3 - (r3'*e1)*(e1 + r1)/(1 + e1'*r1);
      Tr = [e1 e2 e3];

end


end

%%  function Rw ------------------------------------------------------------------------
function R = Gibb2Rmat(w)
% Rotation matrix in terms of the tangent-scaled pseudo-vector

I  = eye(3);
St = Spin(w);
R  = I + (St + St*St/2)/(1 + w'*w/4);

end

%%  function Quat2Gibb -----------------------------------------------------------------
function w = Quat2Gibb(q)
% get the tangent-scaled pseudo-vector w from the quaternion q

w = 2*q(1:3)/q(4);

end

