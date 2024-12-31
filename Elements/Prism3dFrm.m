function ElemResp = Prism3dFrm(action,el_no,xyz,ElemData,ElemState)
%  =========================================================================================
% 
%    Internal functions:
%
%       [k,v0] = EulerOrder00(L,ElemData,v);
%       [k,v0] = ShearOrder00(L,ElemData,v)
%       [k,v0] = AxialPrism3D(L,ElemData,v);
%
%       k = TangentRelease(k,ElemData);
%
%
%  =========================================================================================
%%
ndf = 6;              % no of element DOFs per node (2-node, 3d frame element)
ndm = size(xyz,1);    % element dimension
ElemType.ndm   =       3;
ElemType.ndf   =       6;
ElemType.nsr   =       3;
ElemType.nqv   =       4;
ElemType.nen   =       2;
ElemType.Shear =   false;
ElemType.Release =  true;
ElemType.Type  = 'Frame';
ElemType.Solve =   'Init';

%% report size of element arrays, or check element data;
switch action
%% Report size of element arrays
  case 'size'
    ElemResp = ones(2,ndf);
    return

%% Return element type information 
  case 'type'
    ElemResp = ElemType;
    if isfield(ElemData, "Petrov")
      ElemResp.Petrov = ElemData.Petrov;
    else
      ElemResp.Petrov = false;
    end
    return

%% check element data; assign default values, if necessary
  case 'chec'
    if ~isfield(ElemData,'w'),       ElemData.w       = zeros(ndm,1); end
    if ~isfield(ElemData,'Form'),    ElemData.Form    = 0; end
    if ~isfield(ElemData,'e0'),      ElemData.e0      = zeros(3,1);   end
    if ~isfield(ElemData,'Release'), ElemData.Release = zeros(6,1);   end

    if isfield(ElemData, 'Shear')
      ElemType.Shear = ElemData.Shear;
    end

    ElemData = Check3dFrm(el_no, xyz, ElemData, ElemType);
    ElemResp = ElemData;
    return

  case 'defo'
    ElemResp = str2func('DeformShape3dFrm');    
    return

end

% put joint offset and element orientation information to GeomData
if isfield(ElemData, 'GeomData')
  GeomData = ElemData.GeomData;
end
GeomData.JntOff = ElemData.JntOff;
GeomData.yornt  = ElemData.yornt;
% Extract element loading value w
w = ElemData.w;

% ==========================================================================================

%% Element actions
ElemResp = [];  % if not otherwise specified, ElemResp is empty
switch action
  case 'init'
    %% initialization and specification of history variables
    ElemState.Pres = []; % history array is empty for linear element
    ElemResp       = ElemState;

% ==========================================================================================
  case {'basic'}
    %% basic force-deformation
    L     = ElmLenOr(xyz+GeomData.JntOff);
    v     = ElemState;
    lamda = 1.0; % TODO NaN;
    switch ElemData.Shear
      case 0
        switch ElemData.Form
          case 2
            [k,v0] = AxialPrism3D(L,ElemData,v,lamda);
          otherwise
            [k,v0] = EulerOrder00(L,ElemData,v);
        end
      case 1
        switch ElemData.Form
          case 2
            [k,v0] = AxialPrism3D(L,ElemData,v,lamda);
          case 0   % Standard Timoshenko
            [k,v0] = ShearOrder00(L,ElemData,v);
        end
    end

    % Condense out hinges
    k = TangentRelease(k,ElemData);
    % basic force-deformation relation
    q = k*(v - v0);
    ElemResp.q  = q;
    ElemResp.k  = k;
    return

% ==========================================================================================
  case {'stif','forc'}
    %% state determination
    % undeformed element length
    L   = ElmLenOr(xyz+GeomData.JntOff);
    % extract displacements from ElemState and reshape to array
    nen = size(xyz,2);
    u   = ExtrReshu(ElemState,ndf,nen);
    % transform end displacements from global reference to basic system
    [ag,bg,ab,v] = GeomTran_3dFrm(ElemData.Geom,xyz,GeomData,u);
    %% basic force-deformation
    switch ElemData.Shear
      case 0
        switch ElemData.Form
          case 2
            [k,v0] = AxialPrism3D(L,ElemData,v,ElemState.lamda);

          otherwise
            [k,v0] = EulerOrder00(L,ElemData,v);
        end
      case 1
        switch ElemData.Form
          case 2
            [k,v0] = AxialPrism3D(L,ElemData,v,ElemState.lamda);
          case 0   % Standard Timoshenko
            [k,v0] = ShearOrder00(L,ElemData,v);
        end
    end

    % Condense out hinges
    k = TangentRelease(k,ElemData);

    % basic force-deformation relation
    q = k*(v - v0);

    %% Transform stiffness and forces of basic system to global coordinates
    % Determine equilibrium forces of basic system under element loads
    pbw = [-w(1)*L; -w(2)*L/2; -w(3)*L/2; 0; 0; 0; 
                 0; -w(2)*L/2; -w(3)*L/2; 0; 0; 0];
    % Transform basic forces to global coordinates and add end forces due to w
    p = bg*q + ab'*pbw;
    if strcmp(action,'stif')
      % Determine consistent geometric stiffness matrix
      kg = kg_3dFrm(ElemData.Geom, xyz, GeomData, u, q);
      % Transform stiffness matrix to global coordinates and add geometric stiffness
      ke = ag' * k * ag + kg;
      ElemState.ke = ke;
    end
    ElemState.p = p;
    ElemState.ConvFlag = true;       % element does not involve iterations
    ElemResp = ElemState;

% ==========================================================================================
  case 'mass'
    %% lumped mass vector and consistent mass matrix
    A   = ElemData.A;
    rho = ElemData.rho;
    % Determine element length and orientation (direction cosines)
    L   = ElmLenOr(xyz+GeomData.JntOff);
    % Lumped mass matrix
    tm  = 0.5*rho*A*L.*ones(2*ndm,1);
    ml  = zeros(2*ndf,1);
    ml ([1:3 7:9]) = tm;
    Ip = ElemData.Iy + ElemData.Iz;
    mc = Mass4Displ3dFrm(rho,L,A,Ip);
    ElemMass.ml = ml;
    ElemMass.mc = mc;
    ElemResp    = ElemMass;

% ==========================================================================================
  case 'post'
    %% post-processing information - coordinates of deformed shape
    % undeformed element length
    L   = ElmLenOr(xyz+GeomData.JntOff);
    % extract displacements from ElemState and reshape to array
    nen = size(xyz,2);
    u   = reshape(ElemState.u,ndf,nen);
    % TODO [u,Du,DDu] = ExtrReshu(ElemState,ndf,nen);
    % transform end displacements from global reference to basic system
    [~,~,~,v] = GeomTran_3dFrm(ElemData.Geom,xyz,GeomData,u);
    % basic force-deformation
    switch ElemData.Form
      case 2
        [k,v0] = AxialPrism3D(L,ElemData,v,ElemState.lamda);

      otherwise
        [k,v0] = EulerOrder00(L,ElemData,v);
    end
    % Condense out hinges
    k = TangentRelease(k,ElemData);

    % basic force-deformation relation
    q = k*(v - v0);
    % determine deformations ve in the presence of releases, if any
    A  = ElemData.A;
    Iy = ElemData.Iy;
    Iz = ElemData.Iz;
    J  = ElemData.J;
    E  = ElemData.E;
    G  = ElemData.G;
    % flexibility matrix
    f = blkdiag(L/(E*A),...
                L/(6*E*Iz).*[2 -1; -1 2],...
                L/(G*J),...
                L/(6*E*Iy).*[2 -1; -1 2]);

    % element deformations ve
    ve = f*q + v0;
    % post-processing information
    ElemPost.v  = v;
    ElemPost.q  = q;
    ElemPost.ve = ve;
    ElemResp = ElemPost;
% ==========================================================================================
  otherwise
    %% other actions not supported
    warning('off','backtrace');
    warning('E:W',['>> Element ',num2str(el_no,'%i'),...
            ': Action "',action,'" not supported for ' mfilename ' element']);
    warning('on','backtrace');
end
end

%% ---- function  EulerOrder00 --------------------------------------------------------------
function [k,v0] = EulerOrder00(L,ElemData,v)
% state determination of basic 3d frame element
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  A   = ElemData.A;
  Iy  = ElemData.Iy;
  Iz  = ElemData.Iz;
  Iyz = ElemData.Iyz;
  J   = ElemData.J;
  E   = ElemData.E;
  G   = ElemData.G;

  EA  = E * A;
  EIy = E * Iy;
  EIz = E * Iz;
  GJ  = G * J;

  % set up stiffness matrix in basic system

  k   = [ EA/L       0          0       0       0         0    ;
            0    4*EIz/L    2*EIz/L     0       0         0    ;
            0    2*EIz/L    4*EIz/L     0       0         0    ;
            0        0          0     GJ/L      0         0    ;
            0        0          0       0    4*EIy/L   2*EIy/L ;
            0        0          0       0    2*EIy/L   4*EIy/L ];


  %% initial element deformations due to element loading and non-mechanical effects
  w  = ElemData.w;
  e0 = ElemData.e0;
  v0      = zeros(6,1);
  v0(1)   = w(1)*L*L/(2*EA) + e0(1)*L;
  v0(2:3) = w(2)*L^3/(24*EIz).*[ 1;-1] + e0(2)*L/2.*[-1; 1];
  v0(5:6) = w(3)*L^3/(24*EIy).*[-1; 1] + e0(3)*L/2.*[ 1;-1];
end

%% ---- function BasicShearLE3dFrm --------------------------------------------------------------
function [k,v0] = ShearOrder00(L,ElemData,v)
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%           J  = polar moment of inertia
%           E  = modulus of elasticity
%           G  = shear modulus
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  A  = ElemData.A;   % A  = cross sectional area
  Iy = ElemData.Iy;  % Iy = moment of inertia about y-axis
  Iz = ElemData.Iz;  % Iz = moment of inertia about z-axis
  J  = ElemData.J;
  E  = ElemData.E;
  G  = ElemData.G;
  Az = ElemData.A; % TODO
  Ay = ElemData.A;

  EA  = E * A;
  EIy = E * Iy;
  EIz = E * Iz;
  GJ  = G * J;
  az = 12*EIz/(G*Az*L^2);
  ay = 12*EIy/(G*Ay*L^2);

  % Set up stiffness matrix in basic system
  kzii = (4+az)/(1+az)*EIz/L;
  kyii = (4+ay)/(1+ay)*EIy/L;
  kzij = (2-az)/(1+az)*EIz/L;
  kyij = (2-ay)/(1+ay)*EIy/L;

  %                i      j              i       j
  k   = [ EA/L     0      0      0       0       0 ;
            0    kzii   kzij     0       0       0 ;
            0    kzij   kzii     0       0       0 ;
            0      0      0    GJ/L      0       0 ;
            0      0      0      0     kyii    kyij;
            0      0      0      0     kyij    kyii];


  v0 = zeros(6,1);
end

%% ---- function AxialPrism3D --------------------------------------------------------------
function [k,v0] = AxialPrism3D(L,ElemData,v,lamda)
% state determination of basic 3d frame element
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  % Extract element parameters
  EA = ElemData.E*ElemData.A;
  Iy = ElemData.Iy;
  Iz = ElemData.Iz;
  J  = ElemData.J;
  E  = ElemData.E;
  Iyz = ElemData.Iyz;

  EIy = E * Iy;
  EIz = E * Iz;
  GJ  = ElemData.G * ElemData.J;

  % Axial force
  if ~isfield(ElemData, "N")
    qa = EA/L.*v(1);
  else
    qa = ElemData.N;
  end

  % Torque
  if ~isfield(ElemData, "T")
    T =  GJ.*v(4)/L;
  else
    T = ElemData.T*lamda;
  end

  %
  I = eye(2);
  O = zeros(2,2);

  Dx = [0 -1;
        1  0];

  % Form flexural section stiffness
  km = E*[Iy  Iyz;
          Iyz Iz ];

  % Modify for shear flexibility
  if  ElemData.Shear
    GAy = ElemData.G*ElemData.Ay;
    GAz = ElemData.G*ElemData.Az;
    Phi = [1/GAy  0;
            0    1/GAz];
  else
    Phi = O;
  end

  A  = Dx*(I + qa*Phi);

  C  = Dx*km*A;

  % Set up generalized "Psi" matrix

  P = [1  0; 
       0  1]*qa;
  
  Ci = inv(C);

  Y = [O  I     O      O 
       O  O     I      O
       O  O     O      I
       O  O  -Ci*P  -Ci*Dx*T];

  % Compute matrix exponential and extract submatrices
  eY  = expm(Y*L);
  E12 = eY(1:2,3:4);
  E13 = eY(1:2,5:6);
  E14 = eY(1:2,7:8);

  E22 = eY(3:4,3:4);
  E23 = eY(3:4,5:6);
  E24 = eY(3:4,7:8);

  E32 = eY(5:6,3:4);
  E33 = eY(5:6,5:6);
  E34 = eY(5:6,7:8);

  E42 = eY(7:8,3:4);
  E43 = eY(7:8,5:6);
  E44 = eY(7:8,7:8);

  %  Note that the following blocks should be zero:
  %     eY(7:8,3:4) == 0
  %     eY(7:8,5:6) == 0

  % Condense out c2 using boundary condition; Note
  % c2 = B3 c3 + B4 c4
  B3  = -E12\E13;
  B4  = -E12\E14;

  %
  % Set up stiffness matrix in basic system


  %        c3            c4
  Fa  = [  B3            B4               % u'(0)
           O             I                % u'''(0)
          E22*B3+E23   E22*B4+E24         % u'(L)
          E42*B3+E43   E42*B4+E44];       % u'''(L)

  %        u'     u'''(0)       u'(L)   u'''(L)
  Fb  = [  Dx  -Dx*Phi*Dx*km*A   O      O                % theta(0)
            O      O            Dx    -Dx*Phi*Dx*km*A];  % theta(L)

  Fc = Fb*Fa;

  Kc  =  [-Dx*C                O 
          -Dx*C*(E32*B3+E33) -Dx*C*(E32*B4+E34)];

  kb  = Kc*inv(Fc);


  %     |      |       theta_z     | theta_x |     theta_y       |
  %     |      |     i          j  |         |    i         j    |
  k   = [ EA/L       0          0       0         0         0    ;
            0   -kb(2,2)   -kb(2,4)     0     -kb(2,1)  -kb(2,3) ;  % i theta_z
            0    kb(4,2)    kb(4,4)     0      kb(4,1)   kb(4,3) ;  % j
            0        0          0     GJ/L        0         0    ;  %   theta_x
            0   -kb(1,2)   -kb(1,4)     0     -kb(1,1)  -kb(1,3) ;  % i theta_y
            0    kb(3,2)    kb(3,4)     0      kb(3,1)   kb(3,3) ]; % j


  % TODO
  %% initial element deformations due to element loading and non-mechanical effects
  %    w  = uniform element load ( w(1) = longitudinal, w(2),w(3) = transverse in y and z, resp.)
  %    e0 = initial deformations ( e(1) = axial strain, e(2),e(3) = curvature about y and z, resp.)
  w  = ElemData.w;
  e0 = ElemData.e0;
  v0      = zeros(6,1);
  v0(1)   = w(1)*L*L/(2*EA) + e0(1)*L;
  v0(2:3) = w(2)*L^3/(24*EIz).*[ 1;-1] + e0(2)*L/2.*[-1; 1];
  v0(5:6) = w(3)*L^3/(24*EIy).*[-1; 1] + e0(3)*L/2.*[ 1;-1];

  v0(4) =            0.5*(T/GJ)*L; %*L*L;

end


function k = TangentRelease(k, ElemData)
  %%    compatibility matrix in the presence of axial and/or moment releases
  % introduce release indices MR: 0 indicates no hinge, 1 indicates hinge
  MR = zeros(6,1);
  if isfield(ElemData,'Release'), MR(ElemData.Release==1) = 1; end

  ahi = [ 1-MR(1)       0                    0;
           0            1-MR(2)        -0.5*(1-MR(3))*MR(2);
           0      -0.5*(1-MR(2))*MR(3)       1-MR(3)        ];

  ahj = [ 1-MR(4)       0                    0;
           0            1-MR(5)        -0.5*(1-MR(6))*MR(5);
           0      -0.5*(1-MR(5))*MR(6)       1-MR(6)        ];

  ah  = [  ahi    zeros(3);
         zeros(3)   ahj   ];

  % transform stiffness matrix for the presence of releases
  k   = ah'*k*ah;
end

