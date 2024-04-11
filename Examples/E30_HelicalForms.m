%% Flexure of a 3d cantilever
% Notably used by [1], and generalized by [2] to include out-of-plane force.
% Investigated by [3] with corotational elements.
%
%  [1] Simo JC, Vu-Quoc L (1986) A three-dimensional finite-strain rod model
%      Part II: Computational aspects.
%      Computer Methods in Applied Mechanics and Engineering, 58(1):79–116.
%      https://doi.org/10/b8wd4z
%
%  [2] Ibrahimbegovic, A. (1997) ‘On the choice of finite rotation parameters’,
%      Computer Methods in Applied Mechanics and Engineering, 149(1–4), pp. 49–71.
%      Available at: https://doi.org/10.1016/S0045-7825(97)00059-5.
%
%  [3] Battini, J.-M. and Pacoste, C. (2002) ‘Co-rotational beam elements with
%      warping effects in instability problems’, Computer Methods in Applied
%      Mechanics and Engineering, 191(17–18), pp. 1755–1789. Available at:
%      https://doi.org/10.1016/S0045-7825(01)00352-8.
%
CleanStart
%
Case = "Spiral"; %"Perturb95"; % "Circle1.75"; %1.75"; % "Perturb02"; % "Perturb97"; % "Spiral"; % "Circle2"; %
% Helix(Case, 'Incr', 'B2', 'Incr');
Helix(Case, 'Iter', 'None', 'Iter');
Helix(Case, 'Iter', 'SFIN', 'Iter');
Helix(Case, 'Incr', 'None', 'Incr');
Helix(Case, 'Incr', 'SFIN', 'Incr');
Helix(Case, 'Init', 'None', 'Init');
%
Helix(Case, 'Iter', 'SFIN', 'Incr');
Helix(Case, 'Incr', 'SFIN', 'Iter');
Helix(Case, 'Iter', 'None', 'Incr');
Helix(Case, 'Incr', 'SFIN', 'Iter');
%
function [Model, ElemData, State, Post] = Helix(Case, Interp, Isometry, Param)
% Function to create and analyze a 3D cantilever with point moment and force.

%% Parameters
% input data
nen    = 3;       % number of element nodes

R = eye(3); % ExpSO3([0.2 .3 0.4]);

LoadTransf = "None";
Ctrl = 'no';
switch Case
  case "Circle"
    % Simo & Vu-Quoc
    ne     = 5;       % number of elements
    nostep = 1; % 20;
    loops  = 2;     % number of loops in spiral; scale=2 corresponds to two loops
    Pt = 0.0;
    Pa = 0.0;
    
    L  = 1;
    E  = 1;
    Iy = 2; % 2;
    Iz = 2; %Iy/2;
    A  = 2; % 000; %1;
    G  = 1;
    J  = 2;
    Ctrl = 'no';
    % R = R*ExpSO3([ pi/2 0 0]);

  case "Circle2"
    L  = 10;
    E  = 1e4;
    Iy = 1e-2;
    Iz = 1e-2;
    A  = 1; %
    G  = 1e4;
    J  = 1e-2;
    ne = 5 ; %10;
    nostep =  1;
    loops = 2;
    Pt = 0.000;
    Pa = 0.0;
    Ctrl = 'no';

  case {"Circle0.125", "2.5pi"}
    L  = 10;
    E  = 1e4;
    Iy = 1e-2;
    Iz = 1e-2;
    A  = 1; %
    G  = 1e4;
    J  = 1e-2;
    ne = 10;
    nostep =  1;
    loops = 0.125;
    Pt = 0.000;
    Pa = 0.0;
    Ctrl = 'no';

  case "Circle0.7"
    L  = 10;
    E  = 1e4;
    Iy = 1e-2;
    Iz = 1e-2;
    A  = 1; %
    G  = 1e4;
    J  = 1e-2;
    ne = 5;
    nostep =  1;
    loops = 0.7; %1.75;
    Pt = 0.000;
    Pa = 0.0;
    Ctrl = 'no';

  case {"Perturb95", "Perturb97", "Perturb01", "Perturb02"}
    L  = 10;
    E  = 1e4;
    Iy = 1e-2;
    Iz = 1e-2;
    A  = 1; %
    G  = 1e4;
    J  = 1e-2;
    ne = 10;
    nostep =  1;
    loops = 0.125;
    Pt = 0.0625; % 0.001;
    Pa = 0.0;
    if strcmp(Case, "Perturb97")
      Pt = 0.001;
    end
    if strcmp(Case, "Perturb02")
      nen = 3;
      if strcmp(Param, "Incr")
        LoadTransf = "Iter2Incr";
      end
    end

  case "Spiral"
    L  = 10;
    E  = 1e4;
    Iy = 1e-2;
    Iz = 1e-2;
    A  = 1; %
    G  = 1e4;
    J  = 1e-2;
    ne = 100;
    nostep =  200; % 500; %800
    loops = 10; %*.8/10;
    Pt = 50;
    Pa = 0.0;
    Ctrl = 'yes';
    if strcmp(Param, "Incr")
      LoadTransf = "Iter2Incr";
    end

  case "Spiral2"
    L  = 10;
    E  = 1e4;
    Iy = 1e-2;
    Iz = 1e-2;
    A  = 1; %
    G  = 1e4;
    J  = 1e-2;
    ne = 60;
    nostep =  200; %800
    loops = 10*0.6;
    Pt = 30;
    Pa = 0.0;
    Ctrl = 'yes';
end

EndMoment = loops*E*Iy/L*pi*2;

%% Model definition (no input options)
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
if strcmp(Isometry, "None") && strcmp(Interp, Param)
  [ElemName{1:ne}] = deal('DisplShear3dFrm_wCS');
else
  [ElemName{1:ne}] = deal('GeomWrap3dFrm');
end

Model = Create_Model(XYZ,CON,BOUN,ElemName);

QuatOpt = false;
ElemData = cell(ne);
for el=1:ne
  ElemData{el}.E  =  E;
  ElemData{el}.A  =  A;
  ElemData{el}.Iz = Iz;
  ElemData{el}.Iy = Iy;
  ElemData{el}.G  =  G;
  ElemData{el}.J  =  J;
  ElemData{el}.yornt = R*[0;0;1];

  ElemData{el}.CoroData.Petrov = false;
  ElemData{el}.CoroData.RotName = Param;
  ElemData{el}.CoroData.IsoName = Isometry;
  ElemData{el}.Update          = Interp;
  ElemData{el}.Petrov          = false;
  ElemData{el}.ElemName        = 'DisplShear3dFrm_wCS';

%% Transformation options
  ElemData{el}.Geom = 'basic';
  ElemData{el}.CoroData.Quaternion = QuatOpt;
  ElemData{el}.CoroData.Translate = 4;

end

% check element data and supply default values
ElemData = Structure('chec', Model, ElemData);

% 
% S_InitialStep

%% Strategy
% initialize solution strategy parameters
SolStrat = Initialize_SolStrat;
SolStrat.IncrStrat.Dlam0   = 1./nostep;
SolStrat.IncrStrat.LFCtrl  = Ctrl; % 'yes';
SolStrat.IterStrat.LFCtrl  = Ctrl; %'yes';
SolStrat.IterStrat.maxiter = 60;
SolStrat.IterStrat.tol = 1e-16;
% SolStrat.IterStrat.Type = 'LnSrch';
% SolStrat.IterStrat.StifUpdt = 'no';
% SolStrat.IterStrat.niter4SU = 5;
SolStrat.Debug  = false;
SolStrat.Output = 0;
% SolStrat.Symm = true;

%% Loading
Pe = zeros(nn,6);
Pe(nn,1) = Pa;
Pe(nn,3) = Pt;
Pe(nn,4:6) = [0; 0; EndMoment];

for i = 1:nn, Pe(i,:) = blkdiag(R,R)*Pe(i,:)'; end
Loading = Create_Loading(Model, Pe);
Loading.MomentDOFS   = Model.DOF(end,4:6);
Loading.MomentTransf = LoadTransf;

%% Run analysis
[State, Post, SolStrat] = MultiStep(Model,ElemData,Loading,nostep,SolStrat);

%% Post-processing
% 
Print(Model, SolStrat, Post, R, Interp, Isometry, Param)

if false && SolStrat.Debug
    figure
    hold on
    plot(arrayfun(@(i) i.Debug.IncrState.lamda*10, Post(2:end)), ...
         arrayfun(@(i) i.Debug.IncrState.DU(Model.DOF(end, 4)), Post(2:end)), '.');
    plot(arrayfun(@(i) i.Debug.IncrState.lamda*10, Post(2:end)), ...
         arrayfun(@(i) i.Debug.IncrState.DU(Model.DOF(end, 5)), Post(2:end)), 'x');
    plot(arrayfun(@(i) i.Debug.IncrState.lamda*10, Post(2:end)), ...
         arrayfun(@(i) i.Debug.IncrState.DU(Model.DOF(end, 6)), Post(2:end)), 'o');
    hold off
end

%%
if false(1)
  figure
  hold on
  plot(arrayfun(@(i) i.lamda*10, Post), ...
       arrayfun(@(i) i.U(Model.DOF(end, 3)), Post));
  
  plot(arrayfun(@(i) i.lamda*10, Post), ...
       arrayfun(@(i) i.U(Model.DOF(end, 4)), Post));
  hold off
  
  figure
  plot(arrayfun(@(i) i.lamda*10, Post), ...
       arrayfun(@(i) i.U(Model.DOF(end, 2)), Post));
  %%
  % plot deformed shape of structural model
  Create_Window(0.5, 0.5);
  PlotOpt.MAGF  = 1;
  % Plot_Model(Model);
  XLim = [   0,   L/loops]*1.1;
  YLim = [-L/2,   L/2];
  ZLim = [-L/10,  L/2];
  Draw_3dAxisCross(XLim, YLim, ZLim,PlotOpt);
  % for i=length(Post)
  %   Plot_DeformedStructure(Model,ElemData,Post(i).U(1:Model.nf),Post(i),PlotOpt);
  % end
  Plot_DeformedStructure(Model,ElemData,Post(end).U(1:Model.nf),Post(end),PlotOpt);
end

%% Save outputs
if R ~= eye(3)
  exten = '_R.json';
else
  exten = '.json';
end

if SolStrat.Debug
    fname = ['..\Rendering\json\' ...
             convertStringsToChars(Case) '_'      ...
             ElemData{1}.Update                   ...
             convertStringsToChars(ElemData{1}.CoroData.IsoName) ...
             ElemData{1}.CoroData.RotName      '_' ...
             convertStringsToChars(LoadTransf)    ...
             int2str(ElemData{1}.CoroData.Petrov) ...
             int2str(ElemData{1}.Petrov)      '_' ...
             int2str(nen)                         ...
             exten];
    WriteModelJson(fname, Model, ElemData, Post);
end

Plots.x = arrayfun(@(i) i.lamda, Post);%, ...
Plots.y = {{arrayfun(@(i) i.U(Model.DOF(end, 1)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 2)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 3)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 4)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 5)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 6)), Post)}};

fname = ['..\Plotting\json\' ...
         convertStringsToChars(Case) '_' ...
         ElemData{1}.Update                   ...
         convertStringsToChars(ElemData{1}.CoroData.IsoName) ...
         ElemData{1}.CoroData.RotName      '_' ...
         convertStringsToChars(LoadTransf)    ...
         int2str(ElemData{1}.CoroData.Petrov) ...
         int2str(ElemData{1}.Petrov)   '_nen' ...
         int2str(nen)                  '_nst' ...
         int2str(nostep)                      ...
         exten];

JsonDump(Plots, fname);
end
%%--------------------------------------------------------------------------------
function Print(Model, State, Post, R, Interp, Isometry, Param)

  fprintf("%s%6s %s ", Interp, Isometry, Param);
  if ~State.ConvFlag
    fprintf("\n");
    return
  end

  Ufin = Post(end).U(:);
  for i = 1:length(Model.DOF), Ufin(Model.DOF(i,:)) = blkdiag(R,R)'*Ufin(Model.DOF(i,:)); end
  % Ufin = Ufin(1:Model.nf);
  
  fprintf("%14.3f%4d%14.8g\t%14.8g\t%14.7g\n", ...
    State.avg_iter, State.max_iter, ...
    Ufin(Model.DOF(end,1)), Ufin(Model.DOF(end,2)), Ufin(Model.DOF(end,3)));
end

