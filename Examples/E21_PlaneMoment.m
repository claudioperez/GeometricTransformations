%% Pure flexure of a 3d cantilever
%
%
CleanStart
clc

%% Parameters
% input data
ne     =  5;      % number of elements
nen    =  2;      % number of nodes per element
nostep =  1;      % number of load steps
scale  =  2;      % number of loops in spiral; scale=2 corresponds to two loops

%
L  = 1;
E  = 1;
I  = 2;
A  = 2;
G  = 1;
J  = 2;

%% Model definition (no input options)
EndMoment = -scale*E*I/L*pi*2/nostep;
Pt        = 0.0;

R = eye(3); %ExpSO3([0.0 0  pi/8]);

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

% Element type
% [ElemName{1:ne}] = deal('GeomWrap3dFrm');
[ElemName{1:ne}] = deal('DisplShear3dFrm_wCM');

Model = Create_Model(XYZ,CON,BOUN,ElemName);

ElemData = cell(ne);
for el=1:ne
  ElemData{el}.E  = E;
  ElemData{el}.A  = A;
  ElemData{el}.Iz = I;
  ElemData{el}.Iy = I;
  ElemData{el}.G  = G;
  ElemData{el}.J  = J;

  ElemData{el}.yornt = R*[0;1;0];
  ElemData{el}.JntOff = zeros(3,2);

  ElemData{el}.ElemName =  "DisplShear3dFrm_wCS"; % "Euler"; %

  ElemData{el}.Update           = 'Iter';
  ElemData{el}.CoroData.RotName = 'Iter';
  ElemData{el}.CoroData.Petrov  =  false;
  ElemData{el}.CoroData.IsoName = 'SFIN';

  ElemData{el}.CoroData.Translate = 4;

  % ElemData{el}.Geom = 'linear';

end

% check element data and supply default values
ElemData = Structure('chec',Model,ElemData);

%% Loading
Pe = zeros(nn,6);
Pe(nn,1) = 0;
Pe(nn,2) = 0;
Pe(nn,3) = Pt;
Pe(nn,4:6) = [0 0 EndMoment]';

for i = 1:nn, Pe(i,:) = blkdiag(R,R)*Pe(i,:)'; end
Loading = Create_Loading(Model,Pe);

%% Strategy and solve
% initialize solution strategy parameters
SolStrat = Initialize_SolStrat;
SolStrat.IncrStrat.Dlam0   = 1;
SolStrat.IncrStrat.LFCtrl  = 'no';
SolStrat.IterStrat.LFCtrl  = 'no';
SolStrat.IterStrat.maxiter = 40;
SolStrat.IterStrat.tol = 1e-12;
SolStrat.Debug  = true;
SolStrat.Output = 3;

% Perform multi-step incremental analysis
tic
[State, Post, SolStrat] = MultiStep(Model,ElemData,Loading,nostep,SolStrat);
toc
PrintSolve(SolStrat);

%% Post-processing

Ufin = Post(end).U(:);
for i = 1:nn, Ufin(Model.DOF(i,:)) = blkdiag(R',R')*Ufin(Model.DOF(i,:)); end
Ufin = Ufin(1:Model.nf);

% For comparison against the exact planar theory.
fprintf("Tip planar displacement:\t%g\n", sqrt(Ufin(Model.DOF(end,2))^2 + Ufin(Model.DOF(end,3))^2));

% Print torsional displacement, should be zero.
fprintf("Max torsional displacement:\t%g\n", max(Post(end).U(Model.DOF(:,6))));

% plot deformed shape of structural model
% if SolStrat.ConvFlag
    Create_Window(0.5, 0.5);
    PlotOpt.MAGF  = 1;
    Plot_Model(Model);
    XLim = [0, 1.1*L];
    YLim = [-L/2, L/2];
    ZLim = [-L/10, L/2];
    Draw_3dAxisCross(XLim, YLim, ZLim,PlotOpt);
    PlotOpt.PlNod = 'yes';
%     Plot_Model(Model,Post(end).U(1:Model.nf),PlotOpt);
    Plot_DeformedStructure(Model,ElemData,Post(end).U(1:Model.nf),Post(end),PlotOpt);
% end

% Save output for rendering
WriteModelJson('..\Rendering\json\Helix0.json', Model, ElemData, Post);

