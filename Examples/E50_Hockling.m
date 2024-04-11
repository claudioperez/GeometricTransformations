%% Hockling
%
CleanStart
warning('off','all');
%
% figure
% hold on
[Model, ElemData, State, Post] = Hockling('Iter', 'None', 'Iter');
[Model, ElemData, State, Post] = Hockling('Iter', 'SFIN',   'Iter');
[Model, ElemData, State, Post] = Hockling('Incr', 'None', 'Incr');
[Model, ElemData, State, Post] = Hockling('Incr', 'SFIN',   'Incr');
[Model, ElemData, State, Post] = Hockling('Init', 'SFIN',   'Init');
% plot_solution(Post, Model, 'o');
[Model, ElemData, State, Post] = Hockling('Iter', 'SFIN',   'Incr');
% plot_solution(Post, Model, 'x');
% plot_solution(Post, Model, '-');
% 
% ylabel("torque");
% xlabel("twist");
fprintf("\n");

L = 240;
% plot deformed shape of structural model
Create_Window(0.5, 0.5);
PlotOpt.MAGF  = 1;
Plot_Model(Model);
XLim = [   0, 1.1*L];
YLim = [-L/2,   L/2];
ZLim = [-L/10,  L/2];
Draw_3dAxisCross(XLim, YLim, ZLim,PlotOpt);
for i=1:5:length(Post)
  Plot_DeformedStructure(Model,ElemData,Post(i).U(1:Model.nf),Post(i),PlotOpt);
end
%
warning('on','all');
%
%%--------------------------------------------------------------------------------
function plot_solution(Post, Model, marker)
  scale = 10;
  nostep = 65;
  EI = 0.0833*71240;
  L  = 240;
  EndMoment = scale*EI/L*pi*2/nostep;
  plot(arrayfun(@(post) post.U(Model.DOF(end,4)), Post), ...
       arrayfun(@(post) post.lamda*EndMoment, Post), marker);
end

function [Model, ElemData, State, Post] = Hockling(Interp, Isometry, Param)
%% Parameters
Case = "Hockling";
ne     = 20;       % number of elements
nen    =  2;       % number of element nodes

LoadTransf = "None";

L  = 240;
E  = 71240;
G  = 27190;
I  = 0.0833;
A  = 10;
J  = 2.16;

nostep = 65;
scale  = 10.0;
R = ExpSO3([0, 0, 0.005]); % ExpSO3([0 .0050 0]); % 001

Tol = "Small"; % "Large"; % "Small";
if strcmp(Param, "Incr")
  % LoadTransf = "Iter2Incr";
end

EndMoment = scale*E*I/L*pi*2/nostep;

%% Model definition

% Coordinates
nn = ne*(nen-1)+1;
XYZ( 1:nn,:) = [ linspace(0,L,nn)'  zeros(nn,1)  zeros(nn,1)]*R';

% Connectivity
CON = zeros(ne,nen);
for i=1:ne
  CON(i,:) = 1 + (i-1)*(nen-1):(i-1)*(nen-1)+nen;
end
CON = num2cell(CON, 2);

% Boundary
BOUN = zeros(nn,6);
BOUN( 1,1:6) = ones(1,6);
BOUN(end, [2 3 5 6]) = 1;

% Element type
[ElemName{1:ne}] = deal('GeomWrap3dFrm');
% [ElemName{1:ne}] = deal('Simo3dFrm');
% [ElemName{1:ne}] = deal('DisplShear3dFrm_wCS');

Model = Create_Model(XYZ,CON,BOUN,ElemName);

ElemData = cell(ne);
for el=1:ne
  ElemData{el}.E  = E;
  ElemData{el}.A  = A;
  ElemData{el}.Iz = I;
  ElemData{el}.Iy = I;
  ElemData{el}.G  = G;
  ElemData{el}.J  = J;
  ElemData{el}.yornt = R*[0;0;1];

  ElemData{el}.CoroData.Petrov = false;
  ElemData{el}.CoroData.RotName = Param;
  ElemData{el}.CoroData.IsoName  = Isometry;
  ElemData{el}.Update          = Interp;
  ElemData{el}.Petrov          = false;
  ElemData{el}.ElemName =  'DisplShear3dFrm_wCS'; %  'Euler';  %

  ElemData{el}.CoroData.Translate = 4;
end

% check element data and supply default values
ElemData = Structure('chec', Model, ElemData);


%% Strategy
% Set solution strategy parameters
SolStrat = Initialize_SolStrat;
SolStrat.Output =  0;
SolStrat.Symm   =  0;
SolStrat.Debug  =  true;
SolStrat.IncrStrat.Dlam0   = 1;
SolStrat.IncrStrat.LFCtrl  = 'yes';
SolStrat.IterStrat.LFCtrl  = 'yes';
if strcmp(Tol, "Small")
  if strcmp(ElemData{1}.CoroData.RotName, 'Incr')
    SolStrat.IterStrat.tol     = 1e-8;  % Battini 2002 uses 1e-6
  else
    SolStrat.IterStrat.tol     = 1e-6;
  end
% else
%   SolStrat.IterStrat.tol     = 1e-16;
end
% else
%   if strcmp(ElemData{1}.CoroData.RotName, 'Incr')
%     SolStrat.IterStrat.tol     = 1e-10;
%   else
%     SolStrat.IterStrat.tol     = 1e-7;
%   end
% end

SolStrat.IterStrat.maxiter = 50;

%% Loading
Pe = zeros(nn,6);
Pe(nn,4:6) = R*[EndMoment 0 0]';
% Apply coordinate system rotation to load vector
for i = 1:nn, Pe(i,:) = blkdiag(R,R)*Pe(i,:)'; end
Loading = Create_Loading(Model, Pe);
Loading.MomentDOFS   = Model.DOF(end,4:6);
Loading.MomentTransf = LoadTransf;

%% Run analysis
[State, Post, SolStrat] = MultiStep(Model,ElemData,Loading,nostep,SolStrat);
PrintSolve(SolStrat, "");

% SolStrat.IncrStrat.Dlam0   = 2;
% [State, Post] = MultiStep(Model,ElemData,Loading,50,SolStrat,State,Post);


% SolStrat.IncrStrat.Dlam0   = 0.5;
% [State, Post] = MultiStep(Model,ElemData,Loading, 5,SolStrat,State,Post);

%
% [State, Post] = MultiStep(Model,ElemData,Loading, 10,SolStrat,State,Post);

%% Post-processing

% transform quaternion representation to angles
% if contains(ElemName{1}, 'LE3dFrm')
%   for t=1:length(Post)
%     for n=2:Model.nn
%       Post(t).U(Model.DOF(n,4:6)) = ...
%         TransfQuat(Post(t).U(Model.DOF(n,4:6)),zeros(3,1));
%     end
%   end
% end

Ufin = Post(end).U(:);
for i = 1:nn, Ufin(Model.DOF(i,:)) = blkdiag(R',R')*Ufin(Model.DOF(i,:)); end
Ufin = Ufin(1:Model.nf);



if R ~= eye(3)
  exten = '_R.json';
else
  exten = '.json';
end

if SolStrat.Debug
  % Export model and state information for rendering
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

% Export state information for plotting
Plots.x = arrayfun(@(i) i.lamda, Post);%, ...
Plots.y = {{arrayfun(@(i) i.U(Model.DOF(end, 1)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 2)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 3)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 4)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 5)), Post)} ...
           {arrayfun(@(i) i.U(Model.DOF(end, 6)), Post)}};

fname = ['..\Plotting\json\' ...
         'Hockling_' ...
         convertStringsToChars(ElemName{1}) '_' ...
         ElemData{1}.Update                   ...
         convertStringsToChars(ElemData{1}.CoroData.IsoName) ...
         ElemData{1}.CoroData.RotName      '_' ...
         int2str(ElemData{1}.CoroData.Petrov) ...
         int2str(ElemData{1}.Petrov)   '_ne'  ...
         int2str(ne)                   '_nen' ...
         int2str(nen)                         ...
         exten];
JsonDump(Plots, fname);
end
