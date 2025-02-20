function GState = GeomTran_Push(CoroData, LocalState, ElemState, CoroState, action, ElemType)
%  =========================================================================================
%
nn  =  length(LocalState.u)/6;
ndf  = 6*nn;

p = LocalState.q;
if strcmp(action, 'stif')
  k = LocalState.k;
else
  k = zeros(ndf);
end

[ath, kg] = G3_Push(LocalState, ElemType, p);
p = ath'*p;
k = ath'*k*ai + kg;

CoroState.u = ExtrReshu(LocalState, 6, nn);
[ar, ap, kr, kp] = G2_Push(CoroState, CoroData, p);
p  = ar'*p;
k  = ap'*k*ap + kp;
k  = ar'*k*ar + kr;

[aa, kg] = G1_Push(ElemState, CoroData, p);
p = ap'*p; %%
k = aa'*k*aa + kg;

GState.p  = p;
%GState.ag = ag;

if strcmp(action, 'stif')
  GState.ke = k;
end
end

function [ath, kg] = G3_Push(LocalState, ElemType, p)
ndf = length(p);
nn  = ndf/6;
%% (c) Element parameterization
ath              = eye(ndf);
% Extract element rotational parameter type
switch ElemType.Solve
  case 'Init'
    [uloc, ~, ~] = ExtrReshu(LocalState, 6, nn);
  case 'Incr'
    [~, uloc, ~] = ExtrReshu(LocalState, 6, nn);
  case 'Iter'
    [~, ~, uloc] = ExtrReshu(LocalState, 6, nn);
end

% Loop over rotational blocks
if ~strcmp(ElemType.Solve, 'Iter') && ~ElemType.Petrov
  for i=4:6:ndf
    ath(i:i+2, i:i+2) = dLogSO3(uloc(i:i+2));
  end
end

%% (c) Deformation variation (kgth)
  kg  = zeros(ndf);
  if ~strcmp(ElemType.Solve, 'Iter')
    % Loop over rotational blocks
    for i=4:6:ndf
      kg(i:i+2, i:i+2) = ddExpInvSO3(uloc(i:i+2), p(i:i+2));
    end
  end
end


function [ar, ap, kr, kp] = G2_Push(CoroState, CoroData, pbar)
ndf = length(pbar);
nn = ndf/6;
%% (b) Corotational projector
Tr  = CoroState.CorotTriad;
%CoroState.u = ExtrReshu(LocalState, 6, nn);
%pbar = CoroState.p;
%  spin-fitter and spin lever
[W, Lth] = Tran3dFrm_IsoPush(CoroData, CoroState);
% Projector
ap = eye(ndf) - Lth*W;

% Reference coordinate transormation
ari = repmat({Tr'},nn*2,1);
switch CoroData.Translate
  case {0}
    [ari{1:2:end}] = deal(eye(3));
    ar  = blkdiag(ari{:});
  otherwise
    ar  = blkdiag(ari{:});
end
%% (b) Corotational projector
  if ~CoroData.Petrov
    pbar = ap'*p;
    Pn = zeros(nn*6,3);
    for i=1:nn
      Pn(6*(i-1)+1:6*i-3,1:3) = Spin(pbar(6*(i-1)+1:6*i-3));
    end

    Pnm = Pn;
    for i=1:nn
      Pnm(6*(i-1)+4:6*i,1:3) = Spin(pbar(6*(i-1)+4:6*i));
    end

%   kp = - ar'*W'*Pn'*ap*ar;
    kp = - W'*Pn'*ap;
    kr = - ar'*Pnm*W*ar;

  else
    kp = zeros(ndf);
    kr = zeros(ndf);
  end
  % Reference coordinate transformation 
% p  = ar'*pbar;
% kg = ar'*kg*ar; % OG

end


function [aa, kg] = G1_Push(ElemState, CoroData, p)
ndf = length(p);
nn  = ndf/6;
%% (a) (re)parameterization of rotation update

aa = eye(ndf);
if ~strcmp(CoroData.RotName, 'Iter')
  switch CoroData.RotName
    case 'Init'
      [urot, ~, ~] = ExtrReshu(ElemState, 6, nn);
    case 'Incr'
      [~, urot, ~] = ExtrReshu(ElemState, 6, nn);
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

% ag  = ag*aa;
end

kg  = zeros(ndf);
%% (a) (re)parameterization of rotation update
  if ~strcmp(CoroData.RotName, 'Iter')
%   p = aa'*p;
%   kg  = aa'*kg*aa;

    if CoroData.Quaternion
      kg( 4: 6, 4: 6) = kg( 4: 6, 4: 6) + ddExpSU2(qI, p( 4: 6));
      kg(10:12,10:12) = kg(10:12,10:12) + ddExpSU2(qJ, p(10:12));

    else
      % Loop over rotational blocks
      for i=4:6:ndf
        kg(i:i+2, i:i+2) = kg(i:i+2, i:i+2) + ddExpSO3(urot(i:i+2), p(i:i+2));
      end
    end
  end
end

