function ElemResp = DisplShear3dFrm_wCS(action,el_no,xyz,ElemData,ElemState)
%  Geometrically exact 3D frame element.
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  When the character variable ACTION has one of the following values,
%  the function performs the listed operations and returns the results in ELEMRESP:
%  ACTION = 'size': report size of element arrays
%           'chec': check element property data for omissions and assign default values
%           'defo': report function handle for deformed shape

%           'init': initialize element history variables
%           'forc': report element resisting forces
%           'stif': report element stiffness matrix and resisting forces
%           'post': report post-processing information
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  The data structure ELEMRESP stands for the following data object(s) for each ACTION:
%  ELEMRESP = ARSZ        for action = 'size' 
%  ELEMRESP = ELEMDATA    for action = 'chec'
%  ELEMRESP = FunHandle   for action = 'defo'
%
%  ELEMRESP = ELEMSTATE   for action = 'init'
%  ELEMRESP = ELEMSTATE   for action = 'stif'
%  ELEMRESP = ELEMSTATE   for action = 'forc'
%  ELEMRESP = ELEMPOST    for action = 'post'
%  ELEMRESP is empty      for unsupported keywords
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  ELEMDATA is a data structure with element property information in fields
%         Geom       = character variable for geometric transformation of node variables
%                      (linear, PDelta or corotational) (default=linear)
%         w          = uniform element load ( w(1) = longitudinal, w(2) = transverse )
%         JntOff     = rigid joint offsets in global X and Y at element ends;
%                      column 1 for node i, column 2 for node j
%         nIP        = number of integration points
%         IntTyp     = function name for element integration
%         SecName    = function name for section s-e response
%         SecData{i} = section property data at integration point i (see function with SecName)
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  ELEMSTATE is a data structure with the current element state; it has the fields
%         u     = vector of total element displacements in global reference
%         Du    = vector of element displacement increments from last convergence
%         DDu   = vector of element displacement increments from last iteration
%         ke    = element stiffness matrix in global reference; 
%                 updated under ACTION = 'stif'
%         p     = element resisting force vector in global reference; 
%                 updated under ACTION = 'stif' or 'forc'
%         Past  = element history variables at last converged state
%         Pres  = current element history variables
%         lamda = row vector of current load factor(s)
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Limitations
%       - Static analysis only, mass is not implemented
%       - Does not yet properly handle initially curved geometry. Use of
%         interpolation functions assumes nodes are uniformly spaced.
%       - No joint offsets
%
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%  [1] Simo JC (1985) A finite strain beam formulation. The three-dimensional
%      dynamic problem. Part I.
%      Computer Methods in Applied Mechanics and Engineering, 49(1):55–70.
%      https://doi.org/10.1016/0045-7825(85)90050-7
%
%  [2] Simo JC, Vu-Quoc L (1986) A three-dimensional finite-strain rod model
%      Part II: Computational aspects.
%      Computer Methods in Applied Mechanics and Engineering, 58(1):79–116.
%      https://doi.org/10/b8wd4z
%
%  [3] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric Transformations of Finite
%      Elements" Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------
%
ElemType.ndm   =           3; % Element formulated in 3D
ElemType.ndf   =           6; % Number DOF per node
ElemType.nsr   =           6; % Number of section resultants
ElemType.isr   =         1:6; % Section resultants are [ N Vy Vz T  My  Mz ]
ElemType.nqv   =          12; % 
ElemType.iqv   =        1:12; % q : [ *change iqv* ]
ElemType.nen   =         1:4;
ElemType.Shear =        true; % 
ElemType.Field =     'Displ'; % Displacement formulation

ndf = 6;
% 

switch action
%% Report size of element arrays
  case 'size'
    ElemResp = ones(size(xyz,2), ElemType.ndf);
    return

%% Report element type information 
  case 'type'
    % Delare element type information
    ElemResp = ElemType;
    ElemResp.Solve  = ElemData.Update;
    ElemResp.Petrov = ElemData.Petrov;
    return

%% Check element data; assign default values, if necessary
  case 'chec'
    % Validate user input and 
    ElemResp = Check(el_no, xyz, ElemData, ElemType);
    return

  case 'init'
    LState = Init(xyz, ElemData);
    ElemResp = LState;

    if nargin > 4 && isfield(ElemState, 'u')
      ElemResp.u = ElemState.u;
      if ~isfield(ElemState, 'Du')
        ElemResp.Du  = ElemState.u;
        ElemResp.DDu = ElemState.u;
      else
        ElemResp.Du  = ElemState.Du;
        ElemResp.DDu = ElemState.DDu;
      end
    end
    return

%% State determination; report element resisting forces and/or stiffness
  case {'forc', 'stif'}
    ElemResp = BasicUpdate(action, el_no, xyz, ElemData, ElemState);
    return

%% Report function handle for deformed shape
  case 'defo'
    ElemResp = @DeformedShape;
    return

%% Report post-processing information
  case 'post'
    ElemPost.ve = zeros(6,1);
    ElemResp = ElemPost;
    return

  otherwise
    warning(['Unknown action in ' mfilename ': "' action '"']);

end % switch
end % function

%%
function ElemData = Check(el_no, xyz, ElemData, ElemType)
% Check element configuration and fill data that is constant
% over element's lifetime.
% =======================================================================
% function by Claudio M. Perez                                       2023
% -----------------------------------------------------------------------

  ElemData.nel = size(xyz,2);

  % Uniformly reduced integration by default
  if ~isfield(ElemData, 'nIP'),    ElemData.nIP    = ElemData.nel-1; end
  if ~isfield(ElemData, 'IntTyp'), ElemData.IntTyp = 'Gauss';        end
  if ~isfield(ElemData, 'Update'), ElemData.Update = 'Iter';         end
  if ~isfield(ElemData, 'Petrov'), ElemData.Petrov =  false;         end

  % TODO: Forcing reduced integration
% ElemData.nIP = ElemData.nel - 1;

  %
  % Pre-compute interpolation functions and Gauss points/weights.
  %
  [xIP,  wIP] = feval(ElemData.IntTyp, ElemData.nIP);

  ElemData.xIP = xIP;
  ElemData.wIP = wIP;

  ElemData = Check3dFrm(el_no, xyz, ElemData, ElemType);
  % Create data structure to store shape/Gauss point data
  Shape = cell(ElemData.nIP, ElemData.nel);

  % Indices to re-order output of "Lagrange" function
  idx = [1  3:ElemData.nel  2];

  for l = 1:ElemData.nIP
    val = Lagrange(ElemData.nel-1, 0, ElemData.xIP(l));
    grd = Lagrange(ElemData.nel-1, 1, ElemData.xIP(l));
    
    j = 0;
    for i = ElemData.nodix
      j = j + 1;
      Shape{l,i}.Value  = val(idx(j));
      Shape{l,i}.Grad   = grd(idx(j));
    end
  end

  ElemData.Shape = Shape;

end % function Check

%% ----------------------------------------------------------------------
function [ElemState, ElemData] = Init(xyz, ElemData)
%  Initialize ElemState
% =======================================================================
% function by Claudio M. Perez                                       2023
% -----------------------------------------------------------------------
%
  [~, R0] = DefGeom_3dFrm(xyz(:,[1 end]),ElemData);
  State.Curvature = zeros(3,ElemData.nIP);
% State.Rotations = repmat(R0 , 1, 1, ElemData.nIP);
  for l=1:ElemData.nIP
    State.Rotations{l} = R0; % quaternion(R0, 'rotmat', 'point');
  end

  ElemState.Past = State;
  ElemState.Pres = State;

end % function Init

%% ----------------------------------------------------------------------
function Resp = BasicUpdate(action, ~, XYZ, ElemData, ElemState)
%  Main element state determination procedure
% =======================================================================
% function by Claudio M. Perez                                       2023
% -----------------------------------------------------------------------
%
  ndm = 3;
  ndf = 6;
  Petrov = ElemData.Petrov;

  % extract element properties
  nIP     = ElemData.nIP;
  wIP     = ElemData.wIP;
  nel     = ElemData.nel;
  Shape   = ElemData.Shape;
  SecHndl = ElemData.SecHndl;

  % Permute section resultants
  %      N   Mz  Vz  T  My  Vy
  % Ps = [ 1   0   0   0   0   0 ; % N
  %        0   0   0   0   0   1 ; % Vy
  %        0   0   1   0   0   0 ; % Vz
  %        0   0   0   1   0   0 ; % T
  %        0   0   0   0   1   0 ; % My
  %        0   1   0   0   0   0]; % Mz
  Ps = eye(6);

  % Allocate working arrays
  p      = zeros(ndf*nel, 1);
  kt     = zeros(6*nel, 6*nel);
  Ba     = zeros(6,6,ElemData.nel);
  Bb     = zeros(6,6,ElemData.nel);
  theta  = zeros(3, nIP);
  dtheta = zeros(3, nIP);
  dr     = zeros(3, nIP);


  %
  % Interpolate configuration variables ([2] eqn. 5.1)
  %
  [u, Du, DDu] = ExtrReshu(ElemState, ndf, nel);

  % Compute original element length and Jacobian factor
  L   = ElmLenOr(XYZ(:,[1 end]));
  jac = L/2;

  % Compute deformed node positions
  xyz    = XYZ + u(1:ndm,:);
  switch ElemData.Update
    case 'Iter'
      Past  = ElemState.Pres;

      % Compute dr, theta, dtheta at Gauss points
      for l = 1:nIP
        for i = ElemData.nodix
          dr(:,l)     = dr(:,l)     + Shape{l,i}.Grad./jac.*xyz(:,i);
          theta(:,l)  = theta(:,l)  + Shape{l,i}.Value.*DDu(4:6,i);
          dtheta(:,l) = dtheta(:,l) + Shape{l,i}.Grad./jac.*DDu(4:6,i);
        end
      end

    otherwise
      Past  = ElemState.Past;
      if strcmp(ElemData.Update, 'Init')
        Du = u;
      end

      % Compute dr, theta, dtheta at Gauss points
      for l = 1:nIP
        for i = ElemData.nodix
          dr(:,l)     = dr(:,l)     + Shape{l,i}.Grad./jac.*xyz(:,i);
          theta(:,l)  = theta(:,l)  + Shape{l,i}.Value.*Du(4:6,i);
          dtheta(:,l) = dtheta(:,l) + Shape{l,i}.Grad./jac.*Du(4:6,i);
        end
      end
  end

  %
  % Gauss quadrature loop
  %
  for l = 1:nIP

    % Update the configuration
%   Q = quaternion(theta(:,l)',"rotvec");
    DR = ExpSO3(theta(:,l));
    Pres.Rotations{l} = DR*Past.Rotations{l}; % Q*Past.Rotations{l}; %
    R = Pres.Rotations{l}; % .rotmat('point');

    omega = DR*Past.Curvature(:,l);

    % Spatial curvature, kappa
    Pres.Curvature(:,l)   = omega + dExpSO3(theta(:,l))*dtheta(:,l); % Q.rotmat("point")*Past.Curvature(:,l); % 

    kappa = Pres.Curvature(:,l);

    % Material strain measures
    Gamma = R'*dr(:,l) - [1; 0; 0];
    Kappa = R'*kappa;

    % Section resultants
    SecState.e = Ps*[Gamma; Kappa];
    SecResp = SecHndl('stif',l,ndm,ElemData.SecData{l},SecState);
    Ks = SecResp.ks;
    S  = SecResp.s;

    % Push-forward material resultants
    RR = blkdiag(R,R)*Ps;                    % [2], eqn 4.5c
    s = RR*S;

    %
    % Form response and tangent
    %
    Ba(:,:,:) = 0.0;
    Bb(:,:,:) = 0.0;
    for i = ElemData.nodix
      if Petrov || strcmp(ElemData.Update, 'Iter')
        Ba(:,:,i) = B_nat(Shape{l,i}, dr(:,l), jac);
        if ~strcmp(ElemData.Update, 'Iter')
          Bb(:,:,i) = B_log(Shape{l,i}, dr(:,l), jac,   theta(:,l),dtheta(:,l));
        else
          Bb(:,:,i) = Ba(:,:,i);
        end
      else
        Ba(:,:,i) = B_log(Shape{l,i}, dr(:,l), jac,   theta(:,l),dtheta(:,l));
        Bb(:,:,i) = Ba(:,:,i);
      end
%     (Ba(:,:,i)*RR)'
      p(6*(i-1)+1:6*(i)) = p(6*(i-1)+1:6*(i)) + Ba(:,:,i)*s .* wIP(l);
    end


    if strcmp(action,'stif')
      % Material tangent, [2] eqn 5.6
      % --------------------------------------------------------------

      % Spatial moduli; see [2] in text after eqn 5.6
      ks = RR*Ks*RR';

      for i = ElemData.nodix
        for j = ElemData.nodix
          I = 6*(i-1) + 1;
          J = 6*(j-1) + 1;
          kt(I:I+5, J:J+5) = kt(I:I+5, J:J+5) ...
                           + Ba(:,:,i)*ks*Bb(:,:,j)'.*wIP(l);
        end
      end

      % Geometric tangent
      % --------------------------------------------------------------
      sn = Spin(s(1:3));   % n^ : skew matrix from axial resultants
      sm = Spin(s(4:end)); % m^ : skew matrix from moment resultants
      kg = zeros(6, 6);    %      Geometric stiffness

      % Pre-compute some frequently used quantities
      th = theta(:,l);
      T  = dExpSO3(theta(:,l));
      Tl = T';
      smT = sm*T;
      snT = sn*T;

      for i = ElemData.nodix
        % indices for node i
        ia = 6*(i-1)+1:6*i;

        for j = ElemData.nodix
          kg(:,:)         = 0.0;

          if Petrov || strcmp(ElemData.Update, 'Iter')
            kg(1:3,4:end)   = -sn .* Shape{l,i}.Grad./jac .* Shape{l,j}.Value;
            kg(4:end,1:3)   =  sn .* Shape{l,j}.Grad./jac .* Shape{l,i}.Value;
            kg(4:end,4:end) = ( ...
                -sm .* Shape{l,i}.Grad./jac .* Shape{l,j}.Value + ...
                Spin(dr(:,l)) * sn .* Shape{l,i}.Value .* Shape{l,j}.Value ...
            );

            if Petrov
              kg = kg*[ eye(3)  zeros(3)
                       zeros(3)    T    ];
            end

          else

            kg(1:3,4:end)   = -snT .* Shape{l,i}.Grad./jac * Shape{l,j}.Value;
            kg(4:end,1:3)   = Tl*sn .* Shape{l,j}.Grad./jac * Shape{l,i}.Value;
            kg(4:end,4:end) = ( ...
              (Tl*Spin(dr(:,l))*sn*T)*Shape{l,j}.Value.*Shape{l,i}.Value ...
              - Tl*smT.*Shape{l,i}.Grad./jac.*Shape{l,j}.Value ...
              ...
              + Shape{l,i}.Value.*(dTanSO3(th, sn*dr(:,l), 'L').*Shape{l,j}.Value) ...
              + Shape{l,i}.Grad./jac.*(dTanSO3(th, s(4:6), 'L').*Shape{l,j}.Value) ...
              + Shape{l,i}.Value.*(...
                    dTanSO3(th, S(4:6), 'R')'.*Shape{l,j}.Grad./jac ...
                 + ddTanSO3(th, dtheta(:,l), S(4:6)).*Shape{l,j}.Value ...
              ) ...
            );

          end

          kt(ia, 6*(j-1)+1:6*j) = kt(ia, 6*(j-1)+1:6*j) + wIP(l) .* kg;
        end
      end

    end % action == 'stif'
  end
  %
  % END OF GAUSS LOOP
  %

  Resp = ElemState;
  
  if strcmp(action,'stif')
    Resp.ke = kt.*jac;
  end
  
  Resp.p = p.*jac;
  if ~strcmp(ElemData.Update, 'Init')
    Resp.Pres = Pres;
  end
  Resp.ConvFlag = true;

end % function BasicUpdate

%%  function B_nat ------------------------------------------------------
function B = B_nat(Shape, dx, jac)
% Natural discrete strain matrix, Simo and Vu-Quoc (1986)
% =======================================================================
% function by Claudio M. Perez                                       2023
% -----------------------------------------------------------------------
  B = eye(6) .* Shape.Grad./jac;
  B(4:end,1:3) = -Shape.Value * Spin(dx);
end

%%  function B_log ------------------------------------------------------
function B = B_log(Shape, dx, jac,   th,dth)
% Logarithmic discrete strain matrix, Ibrahimbegovic (1995)
% =======================================================================
% function by Claudio M. Perez                                       2023
% -----------------------------------------------------------------------
  B = eye(6) .* Shape.Grad./jac;
  B(4:end,1:3) = -Shape.Value * Spin(dx);
  Tl = dExpSO3(th)';
  B(4:end,1:3)   = Tl*B(4:end,1:3);
  B(4:end,4:end) = Shape.Value*dTanSO3(th, dth, 'L')'*ExpSO3(th)' ... dot{Tl}
                 + Shape.Grad/jac*Tl; %
% B(4:end,4:end) = Shape.Value*(dTanSO3(th,dth,'R')' - Tt*Spin(kappa)) + Shape.Grad/jac*Tt; % Helix - 45
% B(4:end,4:end) = Shape.Value*(dTanSO3(th,dth,'L')  - Tt*Spin(kappa)) + Shape.Grad/jac*Tt; % Helix - 42
% B(4:end,4:end) = Shape.Value*(ddExpSO3(th,dth)'    - Tt*Spin(kappa)) + Shape.Grad/jac*Tt;
end

%% - function B_lft ------------------------------------------------------------------
function B = B_lft(Shape, dx, jac,R, Kappa, th,dth)
  if nargin < 6
    A = [
      Shape.Grad./jac*R'   Shape.Value.*Spin(R'*dx);
           zeros(3)        Shape.Grad./jac.*eye(3) + Shape.Value.*Spin(Kappa)
    ];
  else
    Tl = dExpSO3(th)';
    A = [
      Shape.Grad./jac*R'   Shape.Value.*Spin(R'*dx)*Tl;
           zeros(3)        (Shape.Grad./jac.*Tl + Shape.Value.*(Spin(Kappa)*Tl + dTanSO3(th, dth)'))
       ... zeros(3)        (Shape.Grad./jac.*Tl + Shape.Value.*ExpSO3(th)*dTanSO3(th, dth, 'L'))'
    ];
  end
  B = A';
end


%% --- function DeformedShape -------------------------------------------
function XYZd = DeformedShape(XYZ,~,u,~,MagF)
% Compute points along the deformed element shape for rendering.
% =======================================================================
% function by Claudio M. Perez                                       2023
% -----------------------------------------------------------------------
  xyz = XYZ + u(1:3,:);

  nn = size(XYZ,2);         % number of element nodes
  idx = [1  3:nn  2];       % node index map for Lagrange function
  nsub = 100;               % number of subdivisions
  XYZd = zeros(3, nsub);    % deformed shape
  s = linspace(-1,1,nsub);  % Coordinates in element's parent domain
  N = Lagrange(nn-1, 0, s); % Shape functions

  % Compute the isoparametric mapping
  for i=1:3
    XYZd(i,:) = N(idx,:)'*xyz(i,:)'*MagF;
  end
end


%%
function Resp = LeftUpdate(action, XYZ, ElemData, ElemState)
  ndm = 3;
  ndf = 6;
  Petrov = false;

  % extract element properties
  nIP     = ElemData.nIP;
  xIP     = ElemData.xIP;
  wIP     = ElemData.wIP;
  nel     = ElemData.nel;
  Shape   = ElemData.Shape;
  SecHndl = ElemData.SecHndl;

  % TODO: Permute section resultants
  %      N   Mz  Vz  T  My  Vy
  % Ps = [ 1   0   0   0   0   0 ; % N
  %        0   0   0   0   0   1 ; % Vy
  %        0   0   1   0   0   0 ; % Vz
  %        0   0   0   1   0   0 ; % T
  %        0   0   0   0   1   0 ; % My
  %        0   1   0   0   0   0]; % Mz
  Ps = eye(6);

  [u, Du, DDu] = ExtrReshu(ElemState, ndf, nel);

  L = ElmLenOr(XYZ(:,[1 end]));
  jac = L/2;


  p  = zeros(ndf*nel, 1);
  kt = zeros(6*nel, 6*nel);
  Ba = zeros(6,6,ElemData.nel);
  Bb = zeros(6,6,ElemData.nel);

  %
  % Interpolate configuration variables ([2] eqn. 5.1)
  %
  xyz    = XYZ + u(1:ndm,:);
  theta  = zeros(3, nIP);
  dtheta = zeros(3, nIP);
  dx     = zeros(3, nIP);
  [~, R0] = DefGeom_3dFrm(XYZ(:,[1 end]),ElemData);

  switch ElemData.Update
    case 'Iter'
      Past  = ElemState.Pres;

      for l = 1:nIP
        for i = ElemData.nodix
          dx(:,l)     = dx(:,l)     + Shape{l,i}.Grad./jac.*xyz(:,i);
          theta(:,l)  = theta(:,l)  + Shape{l,i}.Value.*DDu(4:6,i);
          dtheta(:,l) = dtheta(:,l) + Shape{l,i}.Grad./jac.*DDu(4:6,i);
        end
      end

    otherwise
      Past  = ElemState.Past;
      if strcmp(ElemData.Update, 'Init')
        Du = u;
      end
      for l = 1:nIP
        for i = ElemData.nodix
          dx(:,l)     = dx(:,l)     + Shape{l,i}.Grad./jac.*xyz(:,i);
          theta(:,l)  = theta(:,l)  + Shape{l,i}.Value.*Du(4:6,i);
          dtheta(:,l) = dtheta(:,l) + Shape{l,i}.Grad./jac.*Du(4:6,i);
        end
      end
  end

  for l = 1:nIP

    % Update the configuration
    DR = ExpSO3(theta(:,l));
    Pres.Rotations(:,:,l) = Past.Rotations(:,:,l)*DR;
    R = Pres.Rotations(:,:,l);

    omega = DR'*Past.Curvature(:,l);
    Pres.Curvature(:,l)   = omega ...
                          + dExpSO3(theta(:,l))'*dtheta(:,l);

    % Spatial curvature
    kappa = Pres.Rotations(:,:,l)*Pres.Curvature(:,l);

    % Material strain measures
    Gamma = R'*dx(:,l) - [1; 0; 0];
    Kappa = R'*kappa;

    % Section resultants
    SecState.e = Ps*[Gamma; Kappa];
    SecResp = SecHndl('stif',l,ndm,ElemData.SecData{l},SecState);
    Ks = SecResp.ks;
    S  = SecResp.s;

    % Push-forward material resultants
    RR = blkdiag(R,R)*Ps;                    % [2], eqn 4.5c
    s = RR*S;

    %
    % Form response and tangent
    %
    Ba(:,:,:) = 0.0;
    Bb(:,:,:) = 0.0;
    for i = ElemData.nodix
      if Petrov || strcmp(ElemData.Update, 'Iter')
        Ba(:,:,i) = B_lft(Shape{l,i}, dx(:,l), jac, R, Kappa);
        if ~strcmp(ElemData.Update, 'Iter')
          Bb(:,:,i) = B_lft(Shape{l,i}, dx(:,l), jac, R, Kappa, theta(:,l), dtheta(:,l));
        else
          Bb(:,:,i) = Ba(:,:,i);
        end
      else
        Ba(:,:,i) = B_lft(Shape{l,i}, dx(:,l), jac, R, omega, theta(:,l), dtheta(:,l));
        Bb(:,:,i) = Ba(:,:,i);
      end
      p(6*(i-1)+1:6*(i)) = p(6*(i-1)+1:6*(i)) + Ba(:,:,i)*S .* wIP(l);
    end


    if strcmp(action,'stif')

      % Material tangent, [2] eqn 5.6
      % --------------------------------------------------------------

      for i = ElemData.nodix
        for j = ElemData.nodix
          I = 6*(i-1) + 1;
          J = 6*(j-1) + 1;
          kt(I:I+5, J:J+5) = kt(I:I+5, J:J+5) ...
                           + Ba(:,:,i)*Ks*Bb(:,:,j)'.*wIP(l);
        end
      end

      % Geometric tangent
      % --------------------------------------------------------------
      sn = Spin(s(1:3));
      sm = Spin(s(4:end));
      sN = Spin(S(1:3));
      sM = Spin(S(4:end));
      kg = zeros(6, 6);

      % incremental matrices
      th = theta(:,l);
      T = dExpSO3(theta(:,l));
      sNT = sN*T';
      sMT = sM*T';
      Tl = T';
      Xi = dTanSO3(theta(:,l), dtheta(:,l), 'R');

      for i = ElemData.nodix
        aa = [Shape{l,i}.Grad./jac.*eye(3)  zeros(3)
                    zeros(3)                Shape{l,i}.Value.*eye(3)
                    zeros(3)                Shape{l,i}.Grad./jac.*eye(3)];

        for j = ElemData.nodix
          kg(:,:)         = 0.0;

          if  Petrov || strcmp(ElemData.Update, 'Iter')

            kg(1:3,4:end)   = -R*sN .* Shape{l,i}.Grad./jac * Shape{l,j}.Value;
            kg(4:end,1:3)   =  R*sN .* Shape{l,j}.Grad./jac * Shape{l,i}.Value;
            kg(4:end,4:end) = ( ...
                -R*sM .* Shape{l,i}.Grad./jac .* Shape{l,j}.Value +     ...
                 Spin(dx(:,l)) * sN .* Shape{l,i}.Value .* Shape{l,j}.Value ...
                ... ((R'*dx(:,l)) * S(1:3)' -(R'*dx(:,l))'*S(1:3)*eye(3) ...
                ... + Kappa*S(4:6)' - Kappa'*S(4:6)*eye(3)).* Shape{l,i}.Value .* Shape{l,j}.Value ...
            );

          elseif false

            ab = [Shape{l,j}.Grad./jac.*eye(3)  zeros(3)
                        zeros(3)                Shape{l,j}.Value.*eye(3)
                        zeros(3)                Shape{l,j}.Grad./jac.*eye(3)];

            k2 = zeros(9);

            k2(1:3,4:6) = -R*sN*T';
%           k2(1:3,4:6) = -sn*T;
            k2(4:6,1:3) = k2(1:3,4:6)';

            k2(4:6,4:6) = T*sN*Spin(R'*dx(:,l))*T' + dTanSO3(th,sN*R'*dx(:,l), 'L') ...
                    ... + T*Spin(omega)*sM*T' + ddExpSO3(th,sM*DR'*kappa) + T'*sM*Xi + Xi'*sM*T ...
                        + T*sM*Spin(omega)*T' + dTanSO3(th,sM*omega, 'L') + T'*sM*Xi + Xi'*sM*T ...
                  ; % + H(th,dth,s(4:6));

            % Pi^T
            k2(4:6,7:9) =  T*sM*T' + dTanSO3(th, S(4:6), 'L')';
            k2(7:9,4:6) =  k2(4:6,7:9)';

            kg = aa'*k2*ab;

          else
            kg(1:3,4:end)   = -R*sNT .* Shape{l,i}.Grad./jac * Shape{l,j}.Value;
            kg(4:end,1:3)   =  T*sN*R' .* Shape{l,j}.Grad./jac * Shape{l,i}.Value;
            kg(4:end,4:end) = ( ...
              (T*Spin(R'*dx(:,l))*sNT) .*Shape{l,j}.Value.*Shape{l,i}.Value ...
              - T*sMT .*Shape{l,i}.Grad./jac.*Shape{l,j}.Value ...
              ...
              + Shape{l,i}.Value.*(dTanSO3(th, sN*R'*dx(:,l), 'L').*Shape{l,j}.Value) ...
              + Shape{l,i}.Grad./jac.*(dTanSO3(th, S(4:6), 'L').*Shape{l,j}.Value) ...
              + Shape{l,i}.Value.*(...
                    dTanSO3(th, S(4:6), 'R')'.*Shape{l,j}.Grad./jac ...
                 + ddTanSO3(th, dtheta(:,l), S(4:6)).*Shape{l,j}.Value ...
                ) ...
            );
          end % if

          kt(6*(i-1)+1:6*i, 6*(j-1)+1:6*j) = kt(6*(i-1)+1:6*i, 6*(j-1)+1:6*j) + wIP(l) .* kg;
        end % for j
      end % for i

    end % action == 'stif'
  end
  %
  % END OF GAUSS LOOP
  %

  Resp = ElemState;
  
  I = eye(3);
% RR0 = blkdiag(R0,R0,R0,R0)';
  RR0 = blkdiag(I,R0,I,R0);
  p = RR0*p;
  if strcmp(action,'stif')
    kt = RR0*kt;
    Resp.ke = kt.*jac;
  end

  Resp.p = p.*jac;
  if ~strcmp(ElemData.Update, 'Init')
    Resp.Pres = Pres;
  end
  Resp.ConvFlag = true;

end

