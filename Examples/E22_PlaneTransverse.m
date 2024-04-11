%
%    Batista, M. (2016) ‘A closed-form solution for Reissner planar finite-strain
%    beam using Jacobi elliptic functions’, International Journal of Solids and
%    Structures, 87, pp. 153–166. Available at:
%    https://doi.org/10.1016/j.ijsolstr.2016.02.020.

% Clear workspace memory and close all windows
CleanStart; 
% Hide warnings to avoid cluttering table output.
warning('off','all');

% List of cases and analytical solutions
%          GAv      1 - x          y
Cases = [
        % 5e20  -0.056433236  0.301720774
           500  -0.061315658  0.317813874
            50  -0.103284917  0.465413303
            10  -0.252136606  1.167095878
             5  -0.376121399  2.104087473  ];

for c = Cases'
  GAv = c(1);
  for nen=[2 3 4]
    for param = ["Iter" "Incr"  "Init"]
      Batista(GAv, param, 'None', param, nen);
      Batista(GAv, param, 'SFIN',   param, nen);
    end
  end
  % Print analytic solution
  fprintf("                                  %14.8f %14.8f\n", c(2), c(3));
end

warning('on','all');

%%
function [Model, ElemData, State, Post] = Batista(GAv, Interp, Isometry, Param, nen)
Param = convertStringsToChars(Param);
Interp = convertStringsToChars(Interp);
%% Parameters
Case   = "Batista";
ne     = 10;       % number of elements

Ctrl = 'no';
nostep = 1; % 20;
Pt     = 10;
Pa     = 0.0;

L  =    1;
E  =   10;
Iy =    1; % 2;
Iz =    1; %Iy/2;
A  =   1e7;
G  =   GAv/A;
J  =    A;

R = eye(3); % ExpSO3([0.2 .3 0.4]);

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

ElemData = cell(ne);
for el=1:ne
  ElemData{el}.E  =  E;
  ElemData{el}.A  =  A;
  ElemData{el}.Iz = Iz;
  ElemData{el}.Iy = Iy;
  ElemData{el}.G  =  G;
  ElemData{el}.J  =  J;
  ElemData{el}.yornt = [0;0;1];

  ElemData{el}.CoroData.Petrov = false;
  ElemData{el}.CoroData.RotName = Param;
  ElemData{el}.CoroData.IsoName  = Isometry;
  ElemData{el}.Update            = Interp;
  ElemData{el}.Petrov            = false;
  ElemData{el}.ElemName          = 'DisplShear3dFrm_wCS';
  ElemData{el}.CoroData.Translate = 4;
end

% check element data and supply default values
ElemData = Structure('chec', Model, ElemData);

%% Strategy
% initialize solution strategy parameters
SolStrat = Initialize_SolStrat;
SolStrat.IncrStrat.Dlam0   = 1./nostep;
SolStrat.IncrStrat.LFCtrl  = Ctrl;
SolStrat.IterStrat.LFCtrl  = Ctrl;
SolStrat.IterStrat.maxiter = 30;
SolStrat.IterStrat.tol = 1e-16;
SolStrat.Debug  = false;
SolStrat.Output = 0;

%% Loading
Pe = zeros(nn,6);
Pe(nn,1) = Pa;
Pe(nn,2) = Pt;

for i = 1:nn, Pe(i,:) = blkdiag(R,R)*Pe(i,:)'; end
Loading = Create_Loading(Model, Pe);

%% Run analysis
[State, Post, SolStrat] = MultiStep(Model,ElemData,Loading,nostep,SolStrat);

%% Post-processing
% 
Print(Model, SolStrat, Post, R, Interp, Isometry, Param, nen);

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
if false
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
             int2str(ElemData{1}.CoroData.Petrov) ...
             int2str(ElemData{1}.Petrov)      '_' ...
             int2str(nen)                         ...
             exten];
    WriteModelJson(fname, Model, ElemData, Post);
end

%% Export data to plot with Python
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
         int2str(ElemData{1}.CoroData.Petrov) ...
         int2str(ElemData{1}.Petrov)   '_nen' ...
         int2str(nen)                  '_nst' ...
         int2str(nostep)                      ...
         exten];

JsonDump(Plots, fname);
end
%%--------------------------------------------------------------------------------
function Print(Model, State, Post, R, Interp, Isometry, Param, nen)

  fprintf("%s%6s %s ", Interp, Isometry, Param);
  if ~State.ConvFlag
    fprintf("\n");
    return
  end

  Ufin = Post(end).U(:);
  for i = 1:length(Model.DOF), Ufin(Model.DOF(i,:)) = blkdiag(R,R)'*Ufin(Model.DOF(i,:)); end
  % Ufin = Ufin(1:Model.nf);
  
  fprintf("%4d %8.2f%4d%14.8g\t%14.8g\t%14.8g\n", ...
    nen, State.avg_iter, State.max_iter, ...
    Ufin(Model.DOF(end,1)), Ufin(Model.DOF(end,2)), Ufin(Model.DOF(end,3)));
end

