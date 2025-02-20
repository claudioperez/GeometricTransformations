function [p, k] = GeomTran3dFrm_Push(CoroData, LocalState, ElemState, CoroState, action, ElemType)
% ==========================================================================================
%
%  [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric Transformations of Finite
%      Elements" Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
% ==========================================================================================
% function by Claudio M. Perez                                                          2023
% ------------------------------------------------------------------------------------------
%
  nen  = length(LocalState.u)/6;
  ndf  = 6*nen;

  p = LocalState.q;
  if strcmp(action, 'stif')
    k = LocalState.k;
  else
    k = zeros(ndf);
  end

if strcmp(action, 'forc')
  %
  % do g loop
  %
% (a) Update Rotation Parameters
%---------------------------------------------------------------------------
  ath = G3_Push(LocalState, ElemType, p);
  p = ath'*p;
% (b) Isometric Transformation
%---------------------------------------------------------------------------
  CoroState.u = ExtrReshu(LocalState, 6, nen);
  [ar, ap] = G2_Push(CoroState, CoroData, p);
  if ~CoroData.Petrov
    p = ap'*p;
    p = ar'*p;
  else
    p = ar'*p;
  end
% (c) Interpolation Parameters
%---------------------------------------------------------------------------
  aa = G1_Push(ElemState, CoroData, p);
  p = aa'*p; %%

else
  %
  % do g loop
  %
%% (a) Update Rotation Parameters
%---------------------------------------------------------------------------
  [ath, kg] = G3_Push(LocalState, ElemType, p);
  p = ath'*p;
  k = ath'*k*ath + kg;

%% (b) Isometric Transformation
%---------------------------------------------------------------------------
  CoroState.u = ExtrReshu(LocalState, 6, nen);
  [ar, ap, kr, kp] = G2_Push(CoroState, CoroData, p);
  if ~CoroData.Petrov
    p = ap'*p;
    k = ap'*k*ap + kp;

    p = ar'*p;
    k = ar'*k*ar + kr;

  elseif false
    k = k*ap;

    p = ar'*p;
    k = ar'*k*ar + kr;
  else
%   k = k*ap;

    p = ar'*p;
    k = ar'*k*ar + kr;
  end

%% (c) Interpolation Parameters
%---------------------------------------------------------------------------
  [aa, kg] = G1_Push(ElemState, CoroData, p);
  p = aa'*p; %%
  k = aa'*k*aa + kg;
  % end g loop
end

end % function GeomTran3dFrm_Push


%  (c) Element interpolation
%---------------------------------------------------------------------------
function [ag, kg] = G3_Push(LocalState, ElemType, p)

  ndf = length(p);
  nen  = ndf/6;
  ag  = eye(ndf);
  kg  = zeros(ndf);

  if strcmp(ElemType.Solve, 'Iter') || ElemType.Petrov
    return
  end

% if ~strcmp(ElemType.Solve, 'Iter') && ~ElemType.Petrov

  % Extract element rotational parameter type
  switch ElemType.Solve
    case 'Init'
      [uloc, ~, ~] = ExtrReshu(LocalState, 6, nen);
    case 'Incr'
      [~, uloc, ~] = ExtrReshu(LocalState, 6, nen);
  end

  % Loop over rotational blocks
  j = 1;
  for i=4:6:ndf
    ag(i:i+2, i:i+2) = dLogSO3(uloc(i:i+2));
%   ag(i:i+2, i:i+2) = dLogC90(uloc(i:i+2), LocalState.Triad(j).Rmat);
    j = j + 1;
  end

  if nargout > 1
    % Loop over rotational blocks
    for i=4:6:ndf
      kg(i:i+2, i:i+2) = ddExpInvSO3(uloc(i:i+2), p(i:i+2));
    end
  end
% end
end


%% (b) Isometric projector
function [ar, ap, kr, kp] = G2_Push(CoroState, CoroData, p)
%
% Note: ag = ap*ar
%
  ndf = length(p);
  nen = ndf/6;

  Tr  = CoroState.CorotTriad;
  %CoroState.u = ExtrReshu(LocalState, 6, nen);
  %  spin-fitter and spin lever
  [W, Lth] = Tran3dFrm_IsoPush(CoroData, CoroState);

  % Projector
  ap = eye(ndf) - Lth*W;

  % Reference coordinate transormation
  ari = repmat({Tr'},nen*2,1);
  ar  = blkdiag(ari{:});
  
  if nargout > 2
  % Stiffness kg
    if ~CoroData.Petrov
      pbar = ap'*p;
      Pn = zeros(nen*6,3);
      for i=1:nen
        Pn(6*(i-1)+1:6*i-3,1:3) = Spin(pbar(6*(i-1)+1:6*i-3));
      end

      Pnm = Pn;
      for i=1:nen
        Pnm(6*(i-1)+4:6*i,1:3) = Spin(pbar(6*(i-1)+4:6*i));
      end
      kp = - W'*Pn'*ap;
      kr = - ar'*Pnm*W*ar;

    else
      kp = zeros(ndf);
      kr = zeros(ndf);
    end
  end
end

%% (a) (re)parameterization of rotation update
function [aa, kg] = G1_Push(ElemState, CoroData, p)
  ndf = length(p);
  nen  = ndf/6;

  aa = eye(ndf);
  kg  = zeros(ndf);

% return
  if strcmp(CoroData.RotName, 'Iter')
    return
  end

  switch CoroData.RotName
    case 'Init'
      [urot, ~, ~] = ExtrReshu(ElemState, 6, nen);
    case 'Incr'
      [~, urot, ~] = ExtrReshu(ElemState, 6, nen);
  end

  if CoroData.Quaternion
    qI  = Qvec2Quat(urot( 4: 6)');
    qJ  = Qvec2Quat(urot(10:12)');
    aa(4:6,4:6)    = (ExpSU2(qI) + eye(3))/qI(4);
    aa(10:12,10:12)= (ExpSU2(qJ) + eye(3))/qJ(4);
  else
    % Loop over rotational blocks
    for i=4:6:ndf
      aa(i:i+2, i:i+2) = dExpSO3(urot(i:i+2));
    end
  end

  if nargout > 1
    % Stiffness
    if CoroData.Quaternion
      kg( 4: 6, 4: 6) = kg( 4: 6, 4: 6) + ddExpSU2(qI, p( 4: 6));
      kg(10:12,10:12) = kg(10:12,10:12) + ddExpSU2(qJ, p(10:12));

    else
      % Loop over rotational blocks
      for i=4:6:ndf
        kg(i:i+2, i:i+2) = kg(i:i+2, i:i+2) + dTanSO3(urot(i:i+2), p(i:i+2), 'L');
      end
    end
  end
end

