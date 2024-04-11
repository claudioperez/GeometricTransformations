%% 3d cantilever column model under biaxial tip displacements with axial load

% Clear workspace memory and close all windows
CleanStart

% input data
nen    =  2;      % number of nodes per element
ne     =  5;      % number of elements
nostep = 10;

% Basic parameters
L  =  1;
E  =  1;
A  = 50;
G  = 10;
J  =  1;
% Specify tip load
P = ExpSO3([0 0 pi/2.1])*[1; 0; 0]*-3;

% Rotation matrix for coordinate system
R = eye(3); % ExpSO3([pi/5 pi/6 -pi/8]);

%% Element Data
% specify element type
% [ElemName{1:ne}] = deal('LE3dFrm');
% [ElemName{1:ne}] = deal('GeomWrap3dFrm');
[ElemName{1:ne}] = deal('GeomTran3dFrm'); % Geometrically transformed 3d frm
% [ElemName{1:ne}] = deal('RemoCoroLE3dFrmOrig');
% [ElemName{1:ne}] = deal('GenDInel3dFrm');
% [ElemName{1:ne}] = deal('DisplShear3dFrm_wCM');

ElemData = cell(ne);
for el=1:ne
% Data used by all elements
  ElemData{el}.E  = E;
  ElemData{el}.A  = A;
  ElemData{el}.Iz = 1;
  ElemData{el}.Iy = 10;
  ElemData{el}.G  = G;
  ElemData{el}.J  = J;
  ElemData{el}.nIP = 10;
  ElemData{el}.yornt = R*[0;1;0];
  ElemData{el}.JntOff = zeros(3,2);

%% Transformation data
% Data used by both GeomWrap3dFrm and GeomTran3dFrm
  switch 5
  case 1
  % Most intuitive parameters
    ElemData{el}.CoroData.RotName  = 'Init';
  case 2
  % Nearest to FEDEASLab v5.2
    ElemData{el}.CoroData.RotName  = 'Init';
    ElemData{el}.CoroData.IsoName  =   'B1';
    ElemData{el}.CoroData.Quaternion = true;
  case 3
  % Fastest
    ElemData{el}.CoroData.RotName  = 'Iter';
    ElemData{el}.CoroData.IsoName  =   'B1';
  case 4
  % Robust
    ElemData{el}.CoroData.RotName  = 'Incr';
    ElemData{el}.CoroData.IsoName  =   'B2';
    % ElemData{el}.CoroData.Quaternion = true;
  case 5
  % Remo
    ElemData{el}.CoroData.RotName  = 'Incr';
    ElemData{el}.CoroData.IsoName  =   'C1';
  end
  ElemData{el}.CoroData.Petrov = false;

%% Element-specific data
  % Data for LE3dFrm and GenDInel3dFrm
  ElemData{el}.Geom = 'corotational';

  % Data for GeomWrap3dFrm
  ElemData{el}.ElemName          = "LE3dFrm";
  % ElemData{el}.ElemName          =  "DisplShear3dFrm_wCS";

  % Data for DisplShear3dFrm_wCS
  ElemData{el}.Update            = 'Incr';

  % Data for GeomTran3dFrm and GenDInel3dFrm
  ElemData{el}.BElemTyp          = "BInel3dFrm_wEPLHNMYS"; % "BDInel3dFrm_EBwFF"; % 

  % Data for BInel3dFrm_wEPLHNMYS
  ElemData{el}.Mp = [1e5, 1e5];
  ElemData{el}.Np =  1e5;
end

%% Model Generation
% Generate node coordinates for column of length L
nn = ne*(nen-1)+1;
XYZ( 1:nn,:) = [ linspace(0,L,nn)'  zeros(nn,1)  zeros(nn,1)]*R';

% connectivity array
CON = zeros(ne,nen);
for i=1:ne
  CON(i,:) = 1 + (i-1)*(nen-1):(i-1)*(nen-1)+nen;
end
CON = num2cell(CON, 2);

% boundary conditions
BOUN = zeros(nn,6);
BOUN( 1,1:6) = ones(1,6);

Model = Create_Model(XYZ,CON,BOUN,ElemName);

% check element data and supply default values
ElemData = Structure('chec',Model,ElemData);

% initialize solution strategy parameters
SolStrat = Initialize_SolStrat;
SolStrat.IncrStrat.Dlam0   = 1/nostep;
SolStrat.IterStrat.maxiter = 50;
SolStrat.Debug   = false;
SolStrat.Output  = 0;

%% Loading
Pe = zeros(nn,6);
Pe(nn,1:3) = R*P;

Loading = Create_Loading(Model, Pe);

% perform multi-step incremental analysis
[State, Post, SolStrat] = MultiStep(Model,ElemData,Loading,nostep,SolStrat);
PrintSolve(SolStrat);


%% Post-processing

Ufin = Post(end).U(:);
for i = 1:nn, Ufin(Model.DOF(i,:)) = blkdiag(R,R)'*Ufin(Model.DOF(i,:)); end
% Ufin = Ufin(1:Model.nf);

% For comparison against the exact planar theory.
fprintf("Tip planar displacement:\t%16.12f\n", sqrt(Ufin(Model.DOF(end,2))^2 + Ufin(Model.DOF(end,3))^2));

% if SolStrat.ConvFlag
    % display model
    Create_Window (0.5,0.5);    % open figure window
    PlotOpt.MAGF  = 1;          % magnification factor for deformed shape
    Plot_Model(Model);
    XLim = [0, 1.1*L];
    YLim = [-L/2, L/2];
    ZLim = [-L/10, L/2];
    Draw_3dAxisCross(XLim, YLim, ZLim,PlotOpt);
    % Plot_Model(Model,Ufin,PlotOpt);
    Plot_DeformedStructure(Model,ElemData,Post(end).U(1:Model.nf),Post(end),PlotOpt);
% end
