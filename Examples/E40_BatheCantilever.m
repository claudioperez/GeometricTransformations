%% Bathe's Curved Cantilever
% 
% Clear workspace memory and close all windows
CleanStart; 
% Hide warnings to avoid cluttering table output.
warning('off','all');

% List of cases that will be run;
%        Interp  Isom    Param
Cases = { 
    ... {'Iter',   'B2', 'Iter'} ...
    ... {'Iter', 'None', 'Iter'} ...
    ... {'Iter', 'SFIN', 'Iter'} ...
    ... {'Incr', 'None', 'Incr'} ...
    ... {'Incr', 'SFIN', 'Incr'} ...
    ... {'Init', 'None', 'Init'} ...
    ... {'Init', 'SFIN', 'Init'} ...
        {'None', 'C3',   'Incr'} ...
    ... {'Iter', 'SFIN', 'Init'} ...
    ... {'Init', 'SFIN', 'Iter'} ...
    ... {'Incr', 'SFIN', 'Iter'} ...
    ... {'Iter', 'None', 'Init'}
};

% Rotate global coordinate system.
R = eye(3); %ExpSO3([pi/5 pi/6 -pi/8]);

% Loop over cases
for Config = Cases
  % Loop over different load paths
  for Path=[2 5 3]
    c = Config{1};
    % Run the analysis
    [Model, ElemData, State, Post, SolStrat] = Bathe(Path, R, c{:});
    % Print a row of the results table
    Print(Path, Model, SolStrat, Post, R, c{:});
  end
  fprintf("\n");
end

% Turn warnings back on
warning('on','all');

%%-------------------------------------------------------------------------------
function [Model, ElemData, State, Post, SolStrat] = Bathe(Path, R, Interp, Isometry, Param, Petrov)
arguments
 Path;
 R;
 Interp;
 Isometry
 Param;
 Petrov = [false false];
end
% input data
nen = 2;      % number of nodes per element
ne  = 8;      % number of elements

%
E  = 1e3;
A  = 1e4;
G  = 5e2;
I  = A/12;
J  = 1e5/(12*5);
P  = R*ExpSO3([0 0 0])*[0; 1; 0]*600;

tol = 1e-16; %12;
switch Path
  case 1
    % tol = 1e-7;
    nostep = 3;
    steps = ones(1,nostep)/nostep;
  case 2
    % Simo and Vu-Quoc 1986 and (ii) from Jelenic and Crisfield
    steps = [0.5 0.25 0.25];
  case 3
    nostep = 10;
    steps = ones(1,nostep)/nostep;
  case 4
    nostep = 6;
    steps = ones(1,nostep)/nostep;
  case 5
    % I am not aware of anyone else who has solved with these steps.
    nostep = 8;
    steps = ones(1,nostep)/nostep;
  case 6
    steps = [0.5 0.25 0.125 0.0625 0.0625];
end


% Element type
if  strcmp(Isometry, "None") && strcmp(Interp, Param)
  [ElemName{1:ne}] = deal('DisplShear3dFrm_wCM');
else
  [ElemName{1:ne}] = deal('GeomTran3dFrm');
  % [ElemName{1:ne}] = deal('GeomWrap3dFrm');
  % [ElemName{1:ne}] = deal('RemoCoroLE3dFrmOrig');
end

ElemData = cell(ne);
for el=1:ne
  ElemData{el}.E  = E;
  ElemData{el}.A  = A;
  ElemData{el}.Iz = I;
  ElemData{el}.Iy = I;
  ElemData{el}.G  = G;
  ElemData{el}.J  = J;
  ElemData{el}.yornt = R*[0;1;0];
 
  ElemData{el}.CoroData.Petrov = Petrov(2);
  ElemData{el}.CoroData.RotName = Param;
  ElemData{el}.CoroData.IsoName  = Isometry;
  ElemData{el}.Update          = Interp ;
  ElemData{el}.Petrov          = Petrov(1);
  % Element to be called by the wrapper
  ElemData{el}.ElemName        = "DisplShear3dFrm_wCS";
  ElemData{el}.BElemTyp        = "BInel3dFrm_wEPLHNMYS";

% ElemData{el}.CoroData.Translate = 4; % 3;

  % Data for BInel3dFrm_wEPLHNMYS
  ElemData{el}.Mp = [1e10, 1e10];
  ElemData{el}.Np =  1e10;

end

%% Model Generation
% Coordinates
nn = ne*(nen-1)+1; % number of nodes

rad = 100;
arc = linspace(0, pi/4, nn);
XYZ = zeros(nn,3);
for i=1:nn
  XYZ(i,:) = [rad*(sin(arc(i))) 0 rad*(1 - cos(arc(i)))]*R';
end

% Connectivity
CON = zeros(ne,nen);
for i=1:ne
  CON(i,:) = 1 + (i-1)*(nen-1):(i-1)*(nen-1)+nen;
end
CON = num2cell(CON, 2);

% Boundary
BOUN = zeros(nn,6);
BOUN( 1,1:6) = ones(1,6);

Model = Create_Model(XYZ,CON,BOUN,ElemName);

% check element data and supply default values
ElemData = Structure('chec',Model,ElemData);

%% initialize solution strategy parameters
SolStrat = Initialize_SolStrat;
SolStrat.IterStrat.maxiter = 50;
SolStrat.IterStrat.tol  = tol; % Jelenic numbers are for 1e-7
SolStrat.Debug   = true;
SolStrat.Output  = 0;
% SolStrat.TestName = 'ConvergeJelenic';
% SolStrat.Symm  = 2;
% % 
% S_InitialStep



%% Loading
Pe = zeros(nn,6);
Pe(nn,1:3) = P;

% for i = 1:nn, Pe(i,:) = blkdiag(R,R)*Pe(i,:)'; end
Loading = Create_Loading(Model,Pe);

% perform multi-step incremental analysis
SolStrat.IncrStrat.Dlam0   = steps(1);
[State, Post, SolStrat] = MultiStep(Model,ElemData,Loading,1,SolStrat);

if SolStrat.ConvFlag
  for dlam = steps(2:end)
      SolStrat.IncrStrat.Dlam0   = dlam;
      [State, Post, SolStrat] = MultiStep(Model,ElemData,Loading,1,SolStrat,State,Post);
  end
end

%% Post-processing
%
fname = ['..\Rendering\json\' ...
         'Bathe' int2str(Path) '_'      ...
         ElemData{1}.Update                   ...
         convertStringsToChars(ElemData{1}.CoroData.IsoName) ...
         ElemData{1}.CoroData.RotName      '_' ...
         int2str(ElemData{1}.CoroData.Petrov) ...
         int2str(ElemData{1}.Petrov)      '_' ...
         int2str(nen)                         ...
         '.json'];
WriteModelJson(fname, Model, ElemData, Post);

if false(1)
    % display model
    Create_Window(0.5,0.5);    % open figure window
    PlotOpt.MAGF  = 1;          % magnification factor for deformed shape
    Plot_Model(Model);
    XLim = [0, 1.1*L];
    YLim = [-L/2, L/2];
    ZLim = [-L/10, L/2];
    Draw_3dAxisCross(XLim, YLim, ZLim,PlotOpt);
    % Plot_Model(Model,Ufin,PlotOpt);
    Plot_DeformedStructure(Model,ElemData,Post(end).U(1:Model.nf),Post(end),PlotOpt);
end
end

%%--------------------------------------------------------------------------------
function Print(Path, Model, State, Post, R, Interp, Isometry, Param)

  fprintf("%s%6s %s %i ", Interp, Isometry, Param, Path);
  if ~State.ConvFlag
    fprintf(" %i \n",length(Post));
    return
  end

  Ufin = Post(end).U(:);
  for i = 1:length(Model.DOF), Ufin(Model.DOF(i,:)) = blkdiag(R,R)'*Ufin(Model.DOF(i,:)); end


  fprintf("%14.2f%4d %15.8g %15.8g %15.8g\n", ...
    State.avg_iter, State.max_iter, ...
    Ufin(Model.DOF(end,1)), Ufin(Model.DOF(end,2)), Ufin(Model.DOF(end,3)));
end
