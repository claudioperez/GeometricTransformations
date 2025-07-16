function ElemResp = GeomTran3dFrm (action,el_no,xyz,ElemData,ElemState)
%
% GeomTran3dFrm Geometrically transformed 3d frame element.
%
%  =========================================================================================
%  function added                                                                    03-2024                                                       
%
ndf = 6;       % no. of element DOFs per node (2-node, 3d frame element)
ElemType.ndm   =       3;
ElemType.ndf   =       6;
ElemType.nsr   =       3;
ElemType.isr   = [1 5 6]; % s := [ N Mz My ];
ElemType.nqv   =       4;
ElemType.nen   =       2;
ElemType.Shear =   false;
ElemType.Release =  true;
ElemType.Type  = 'Frame';

%% report size of element arrays, check element data; else retrieve element data
switch action
  case 'size'
    arsz = ones(2,ndf);
    ElemResp = arsz;  % return size of element arrays
    return
  case 'chec'
    %% check element data; assign default values, if necessary
    if ~isfield(ElemData,'yornt'), disp('Element');disp(el_no);
                                   error('y-axis orientation missing'); end
    if ~isfield(ElemData,'rho'),      ElemData.rho    = 0;              end
    if ~isfield(ElemData,'JntOff'),   ElemData.JntOff = zeros(3,2);     end
    if ~isfield(ElemData,'LdIdx'),    ElemData.LdIdx  = 0;              end
    if ~isfield(ElemData,'LdIdy'),    ElemData.LdIdy  = 0;              end
    if ~isfield(ElemData,'LdIdz'),    ElemData.LdIdz  = 0;              end
    if ~isfield(ElemData,'rho'),      ElemData.rho    = 0;              end
    if ~isfield(ElemData,'SubDivNo'), ElemData.SubDivNo = 5;            end

    % if Geom is provided, issue a warning
    if  isfield(ElemData, 'Geom')
      warning('off','backtrace');
      warning('E:W',[
        '>> Element ', num2str(el_no,'%i'), ':  ',mfilename,' - unused field Geom="',ElemData.Geom,'".'
      ]);
      warning('on','backtrace');
    end

    % raise error if more than two nodes are supplied
    if size(xyz, 2) ~= 2
      error([
        '>> Element ', num2str(el_no,'%i'), ':  ', mfilename,' only supports 2 nodes.'
      ]);
    end

    % select default basic element type
    if ~isfield(ElemData,'BElemTyp'), ElemData.BElemTyp = 'BInel3dFrm_wEPLHNMYS'; end

    % Transformation type and data
    ElemData.ElemType = ElemType;
    ElemData.ElemType.Petrov = false;
    ElemData.ElemType.Solve = 'Init';
    ElemData.CoroData.Translate = 5;
    ElemData.ndm = 3;

    % Allow BDInel* elements to be initialized from elastic element data
    if contains(ElemData.BElemTyp, 'BDInel')
      ElemData.ElemType.Field = true;
    end

    % Check standard frame data
    ElemData = Check3dFrm(el_no, xyz, ElemData, ElemData.ElemType);

    % check for properties of basic element
    BElemTyp = str2func(ElemData.BElemTyp);
    L = ElmLenOr(xyz+ElemData.JntOff);
    ElemData = BElemTyp('chec',L,ElemData);
    ElemResp = ElemData;
    return

  case 'defo'
    %% return function handle for deformed shape
    ElemResp = str2func('DeformShape3dFrm');
    return
end

% function handle for basic element
BElemTyp = str2func(ElemData.BElemTyp);

% rigid joint offsets and y orientation
GeomData.JntOff = ElemData.JntOff;
GeomData.yornt  = ElemData.yornt;

% undeformed element length
L = ElmLenOr(xyz+GeomData.JntOff);

%% element actions
ElemResp = [];  % if not otherwise specified, ElemResp is empty

switch action
  case 'init'
    %% initialization and specification of history variables
    ElemState.Pres.Elem = BElemTyp('init', L, ElemData).Pres;

    % Initialize coordinate system
    nen = size(xyz,2);
    ElemState.Pres.Local.v    = zeros(ElemType.nqv,1);
    ElemState.Pres.Local.w    = zeros(1,6*nen);
    ElemState.Pres.Local.u    = zeros(1,6*nen);
    ElemState.Pres.Local.Du   = zeros(1,6*nen);
    ElemState.Pres.Local.DDu  = zeros(1,6*nen);
    for i=1:nen
      ElemState.Pres.Local.Triad(i).Rmat = eye(3);
    end

    ElemState.Pres.Rotations = RotationInit(xyz, ElemData);
    % set up transformation matrix from global to local coordinates
    ElemState.Pres.Geom.CorotTriad = ElemState.Pres.Rotations(:,:,1);

    ElemResp  = ElemState;
    
  case {'stif','forc'}
    %% state determination

    % Update and pull state into corotational frame
    [Pres, BasicState] = GeomTran3dFrm_Pull(ElemData.CoroData, xyz, ElemState, ElemData.ElemType);
    Pres.Local      = BasicState;
    % matrix of rigid body modes;
    %      N   Vy   Vz   T    My   Mz  | N   Vy   Vz    T   My   Mz
    av = [-1    0    0    0    0    0    1    0    0    0    0    0   % N
           0    0    0    0    0    1    0    0    0    0    0    0   % Mz
           0    0    0    0    0    0    0    0    0    0    0    1   % Mz
           0    0    0   -1    0    0    0    0    0    1    0    0   % T
           0    0    0    0    1    0    0    0    0    0    0    0   % My
           0    0    0    0    0    0    0    0    0    0    1    0]; % My

    BasicState.v    = av*Pres.Local.u(:);
    BasicState.Dv   = av*Pres.Local.Du(:);
    BasicState.DDv  = av*Pres.Local.DDu(:);
    BasicState.Pres = ElemState.Pres.Elem;
    BasicState.Past = ElemState.Past.Elem;

    BasicState = BElemTyp('stif',L,ElemData,BasicState);

    Pres.Elem = BasicState.Pres;
    ElemState.ConvFlag = BasicState.ConvFlag;

    BasicState.q = av'*BasicState.q;
    BasicState.k = av'*BasicState.k*av;

    % Push response back to current configuration in global coordinates
    [p, ke] = GeomTran3dFrm_Push(ElemData.CoroData, BasicState, ElemState, Pres.Geom, action, ElemData.ElemType);

    % Copy over response quantities
    ElemState.p        = p;
    if strcmp(action,'stif')
      ElemState.ke = ke;
    end

    % Commit update to history
    ElemState.Pres = Pres;
    if strcmp(ElemData.CoroData.RotName, 'Init')
      ElemState.Pres.Rotations = ElemState.Past.Rotations;
    end
    ElemResp = ElemState;

% ==========================================================================================
  case 'post'
    %% post-processing information - coordinates of deformed shape
    % matrix of rigid body modes
    av = [-1    0    0    0    0    0    1    0    0    0    0    0
           0    0    0    0    0    1    0    0    0    0    0    0
           0    0    0    0    0    0    0    0    0    0    0    1
           0    0    0   -1    0    0    0    0    0    1    0    0
           0    0    0    0    1    0    0    0    0    0    0    0
           0    0    0    0    0    0    0    0    0    0    1    0];
    % extract displacements from ElemState and reshape to array
    BasicState = ElemState.Pres.Local;
    BasicState.v    = av*ElemState.Pres.Local.u(:);
    BasicState.Dv   = av*ElemState.Pres.Local.Du(:);
    BasicState.DDv  = av*ElemState.Pres.Local.DDu(:);
    BasicState.Pres = ElemState.Pres.Elem;
    BasicState.Past = ElemState.Past.Elem;
    %% post-processing of basic force-deformation
    % call state determination for effective q and k under given v
    ElemPost = BElemTyp('post',L,ElemData,BasicState);
    ElemPost.v  = BasicState.v;   

    ElemResp = ElemPost;
    return
    
% ==========================================================================================
  otherwise
    %% other actions not supported
    warning('off','backtrace');
    warning('E:W',['>> Element ',num2str(el_no,'%i'),...
            ': Action "',action,'" not supported for ',mfilename,' element']);
    warning('on','backtrace');
end

end % end main function GeomTran3dFrm
