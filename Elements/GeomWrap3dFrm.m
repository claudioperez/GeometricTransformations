function ElemResp = GeomWrap3dFrm (action,el_no,xyz,ElemData,ElemState)
%  ==========================================================================================
%
%  3D Beam Wrapper
%
%  [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric 
%      Transformations of Finite Elements" 
%      Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
%  ==========================================================================================
ElemType.ndm   =       3;
ElemType.ndf   =       6;
ElemType.nsr   =       3;
ElemType.nqv   =       6;
ElemType.nen   =       2;
ElemType.Shear =   false;
ElemType.Solve =  'Init';
ElemType.Petrov =  false;

switch action
  case 'type'
    ElemResp = ElemType;
    return

  case 'size'
    ElemResp = ones(size(xyz,2),6); % num_nodes x ndf
    return

  case 'chec'

    if  isfield(ElemData,'Geom') && ~strcmp(ElemData.Geom, 'basic')
      warning('off','backtrace');
      warning('E:W',['>> Element ',num2str(el_no,'%i'),': Field Geom="',ElemData.Geom,'" overridden with "basic"']);
      warning('on','backtrace');
      ElemData.Geom = 'basic';
      ElemData.CoroData.Translate = 5;

    elseif isfield(ElemData,'Geom') && strcmp(ElemData.Geom, 'basic') ...
       && ~isfield(ElemData.CoroData, 'Translate')
      ElemData.CoroData.Translate = 5;
    end
    if ~isfield(ElemData.CoroData, 'Translate')
      ElemData.CoroData.Translate = 4;
    end

    % Set basic element type
    if ~isfield(ElemData, 'ElemName')
      error([
          '>> Element ', num2str(el_no,'%i'), ':  ', mfilename, ' missing field "ElemName".'
      ]);
    end
    BasicElement = ElemData.ElemName;
    ElemData.ElemName = convertStringsToChars(BasicElement);

    % initial reference frame
    T0 = RotationInit(xyz, ElemData);
    switch ElemData.CoroData.Translate
      case {1, 3, 5, 6}
      ElemData.XYZ = T0(:,:,1)'*(xyz - xyz(:,1));
      case {2, 4}
      ElemData.XYZ = T0(:,:,1)'*xyz;
      otherwise
      ElemData.XYZ = xyz;
    end

    % call 'chec' on  wrapped element
    XYZ   = ElemData.XYZ;
    ElemData = feval(BasicElement, 'chec', el_no, XYZ, ElemData);
    ElemData.ElemType = feval(BasicElement, 'type', el_no, XYZ, ElemData, ElemType);
    if isempty(ElemData.ElemType)
      ElemData.ElemType = ElemType;
    end
    
    % Fill in remaining frame element defaults
    ElemResp = Check3dFrm(el_no, xyz, ElemData, ElemData.ElemType);
    return

  case 'defo'
    if size(xyz,2) == 2
      ElemResp = str2func('DeformShape3dFrm');
    else
      ElemResp = feval(ElemData.ElemName, 'defo', el_no, xyz, ElemData, ElemType);
    end
    return

  case 'post'
    %% post-processing information - coordinates of deformed shape
    ElemPost.ve = zeros(6,1);
    ElemResp = ElemPost;
    return
end


BasicElement = ElemData.ElemName;

switch action

  case 'init'
%% initialization and specification of history variables

    % TODO: Clean, and move to GeomTran3dFrm_Init
    ElemState.Pres.Rotations = RotationInit(xyz, ElemData);
    % set up transformation matrix from global to local coordinates
    ElemState.Pres.Geom.CorotTriad = ElemState.Pres.Rotations(:,:,1);

    % TODO: Clean this up
    XYZ   = ElemData.XYZ; % ElemState.Pres.Geom.CorotTriad'*(xyz - xyz(:,1));
    ElemData.yornt = [0; 1; 0];
    ElemState.Pres.Elem = feval(BasicElement, 'init', el_no, XYZ, ElemData).Pres;

    nen = size(xyz,2);
    ElemState.Pres.Local.v    = zeros(ElemType.nqv,1);
    ElemState.Pres.Local.w    = zeros(1,6*nen);
    ElemState.Pres.Local.u    = zeros(1,6*nen);
    ElemState.Pres.Local.Du   = zeros(1,6*nen);
    ElemState.Pres.Local.DDu  = zeros(1,6*nen);

    for i=1:nen
      ElemState.Pres.Local.Triad(i).Rmat = eye(3);
    end

    ElemResp      = ElemState;
    ElemResp.Past = ElemState.Pres;
    return

  case {'stif','forc'}
%% State determination

    % Update and transform state 
    [Pres, Local] = GeomTran3dFrm_Pull(ElemData.CoroData, xyz, ElemState, ElemData.ElemType);
    Pres.Local = Local;
    Local.Pres = ElemState.Pres.Elem;
    Local.Past = ElemState.Past.Elem;

    % Invoke the wrapped element
    Local = feval(BasicElement, action, el_no, Local.xyz, ElemData, Local);
    Local.q     = Local.p;
    if strcmp(action, 'stif')
      Local.k     = Local.ke;
    end

    % Push response back to current configuration in global coordinates
    [p, ke] = GeomTran3dFrm_Push(ElemData.CoroData, Local, ElemState, Pres.Geom, action, ElemData.ElemType);

    % Copy over response quantities
    ElemState.p        = p;
    ElemState.ConvFlag = Local.ConvFlag;
    if strcmp(action,'stif')
      ElemState.ke = ke;
    end

    % Commit update to history
    if ~strcmp(ElemData.ElemType.Solve, 'Init')
      Pres.Elem = Local.Pres;
    end

    ElemState.Pres = Pres;
    if strcmp(ElemData.CoroData.RotName, 'Init')
      ElemState.Pres.Rotations = ElemState.Past.Rotations;
    end
    ElemResp = ElemState;

  otherwise
    % Unknown action

end % switch action

end % function


%% ---- function BasicShearLE3dFrm --------------------------------------------------------------
function BElemResp = BasicShearLE3dFrm(L,ElemData,BElemState)
% Basic response for linear elastic frame with shear.
% =======================================================================
% function by Claudio M. Perez                                       2023
% -----------------------------------------------------------------------

  v  = BElemState.v;

  % Extract element properties
  A  = ElemData.A;
  Iy = ElemData.Iy;
  Iz = ElemData.Iz;
  J  = ElemData.J;
  E  = ElemData.E;
  G  = ElemData.G;
  Az = ElemData.A; % TODO: This should be scaled down by correction factor
  Ay = ElemData.A;

  EA  = E*A;
  EIy = E*Iy;
  EIz = E*Iz;
  GJ  = G*J;
  az  = 12*EIz/(G*Az*L^2);
  ay  = 12*EIy/(G*Ay*L^2);

  kzii = (4+az)/(1+az)*EIz/L;
  kyii = (4+ay)/(1+ay)*EIy/L;
  kzij = (2-az)/(1+az)*EIz/L;
  kyij = (2-ay)/(1+ay)*EIy/L;

  % Set up stiffness matrix in basic system
  %                i      j              i       j
  k   = [ EA/L     0      0      0       0       0 ;
            0    kzii   kzij     0       0       0 ;
            0    kzij   kzii     0       0       0 ;
            0      0      0    GJ/L      0       0 ;
            0      0      0      0     kyii    kyij;
            0      0      0      0     kyij    kyii];

  %% compatibility matrix in the presence of axial and/or moment releases
  % introduce release indices MR: 0 indicates no hinge, 1 indicates hinge
  MR = zeros(6,1);
  if isfield(ElemData,'Release'), MR(ElemData.Release==1) = 1; end

  ah1 = [ 1-MR(1)       0                    0;
           0            1-MR(2)        -0.5*(1-MR(3))*MR(2);
           0      -0.5*(1-MR(2))*MR(3)       1-MR(3)        ];
  ah2 = [ 1-MR(4)       0                    0;
           0            1-MR(5)        -0.5*(1-MR(6))*MR(5);
           0      -0.5*(1-MR(5))*MR(6)       1-MR(6)        ];
  ah  = [ah1 zeros(3); zeros(3) ah2];

  % transform stiffness matrix for the presence of releases
  BElemResp.k   = ah'*k*ah;

  % basic force-deformation relation
  BElemResp.q = k*v;
  BElemResp.ConvFlag = true;
end


function Local = BasicEulerLE3dFrm(BasicElement, action, el_no, xyz, ElemData, Local);
% =======================================================================
% function by Claudio M. Perez                                       2023
% -----------------------------------------------------------------------

  % Extract element properties
  E   = ElemData.E;       % Young modulus
  G   = ElemData.G;       % Shear modulus
  GJ  = G*ElemData.J;     % Torsional
  EA  = E*ElemData.A;     % Axial
  EIy = E*ElemData.Iy;    % Flexural
  EIz = E*ElemData.Iz;    % Flexural
  L   = ElmLenOr(xyz);

  av = [-1    0    0    0    0    0    1    0    0    0    0    0
         0    0    0    0    0    1    0    0    0    0    0    0
         0    0    0    0    0    0    0    0    0    0    0    1
         0    0    0   -1    0    0    0    0    0    1    0    0
         0    0    0    0    1    0    0    0    0    0    0    0
         0    0    0    0    0    0    0    0    0    0    1    0];

  % set up the local stiffness matrix
  Local.k = [ EA/L       0          0       0       0         0    ;
                0    4*EIz/L    2*EIz/L     0       0         0    ;
                0    2*EIz/L    4*EIz/L     0       0         0    ;
                0        0          0      GJ/L     0         0    ;
                0        0          0       0    4*EIy/L   2*EIy/L ;
                0        0          0       0    2*EIy/L   4*EIy/L ];

  % Compute the local forces
  Local.q = av'*Local.k*Local.v;
  Local.k = av'*Local.k*av;
  Local.ConvFlag = true;
end

