
% L = 60;
% E = 29e3;
% I = 110;
% A = 112;
% G = 11.2e3;

L  = 240;
E  = 71240;
G  = 27190;
I  = 0.0833;
A  = 10;
J  = 2.16;

ElemData.E   = E;
ElemData.G   = G;
ElemData.J   = 2*I;
ElemData.I   = I;
%ElemData.Iy  = I;
%ElemData.Iz  = I*10;
ElemData.A   = A;
ElemData.Ay  = 3/6*A;
ElemData.Az  = 3/6*A;
ElemData.Form  = 2;
ElemData.Shear = 1;
ElemData.Goption = 'exact';
%ElemData.Geom = 'corotational';


Pref = [-1 zeros(1,5)];
%% Pin-Pin
BOUN = [1 1 1   1 1 0;
        0 1 1   1 1 0];
Pcrit = E*I*pi^2/(L*1)^2;     %   pin-pin
%Pcrit = E*I*pi^2/(L*1)^2;    %   fix-slide
%Pcrit = E*I*pi^2/(L*0.5)^2;  %   fix-fix
%Pcrit = E*I*pi^2/(L*0.7)^2;  %   fix-pin
%Pcrit = E*I*pi^2/(L*2)^2;    %   fix-free
%Pcrit = E*I*pi^2/(L*2)^2;    %   pin-slide

ne =  10;
%% Hockling
BOUN = [1 1 1   1 1 1;
        1 1 1   0 1 1];
Pcrit = 8.988*E*I/(L);
Pref = [0 0 0   1    0 0];

%%
%       D

ElemData.nen = 3;
[Load, State] = SBuckleLoad("ExactShear3dFrm", ElemData, BOUN, L, Pcrit, ne, Pref);
%[Load, State] = SBuckleLoad("LE2dFrm_wPdelta", ElemData, BOUN, L, Pcrit);
ElemData.nen = 2;
[Load, State] = SBuckleLoad("Prism3dFrm", ElemData, BOUN, L, Pcrit, ne, Pref);

%%
function [Load, State] = SBuckleLoad(ElemName, ElemOpts, BOUN, L, Pcrit, ne, Pref)
% Perform an element stability analysis. Two procedures are performed:
%   1. A linearized stability analysis
%   2. A limit analysis
%
  if nargin < 5
    ne = 1;
  end
  ProbOpts.L   = L;
  ProbOpts.nen = ElemOpts.nen;
  ProbOpts.ne  = ne;

  % Create the column
  [Model,ElemData] = CreateColumn(ElemName, ElemOpts, ProbOpts, BOUN);

  %% 1. Step: set up the linear stiffness matrix Kl
  State  = Initialize_State(Model,ElemData);
  State  = Structure('stif',Model,ElemData,State);
  Kl = full(State.Kf);

  State = Analyze(Model, ElemData,  Pref, 1);
  Kg = (full(State.Kf)-Kl)/norm(Pref);
  
  %% 4. Step: solve eigenvalue problem
  [EigVec,EigVal] = eig(Kl,Kg);
  % eigenvalues are reported in a diagonal matrix: extract diagonal and sort
  [EigVal,iseig] = sort(abs(diag(EigVal)));
  
  % Relevant eigenvalue is the smallest (lowest buckling load)
  Load = min(abs(EigVal))/Pcrit % *norm(Pref)
  ife = isfinite(EigVal);
% EigVal(ife)*norm(Pref)/Pcrit

  return

  for P = linspace(0, Pcrit, 100)

    State = Analyze(Model, ElemData,  P*Pref, 1);
    K = full(State.Kf);
    ev = min(eig(K));
    ev = det(K);
%   fprintf("%f\t%8.3g\t%8.3g\t%8.3g\n", P/Pcrit,  ev, cond(K), det(K));

    if ev <= 0.0 % || imag(ev) ~= 0
      Limit = P/Pcrit
      break;
    end
  end

  Create_Window(0.5,0.5);    % open figure window
  PlotOpt.MAGF  = 1;          % magnification factor for deformed shape
  Plot_Model(Model);
  XLim = [0, 1.1*L];
  YLim = [-L/2, L/2];
  ZLim = [-L/10, L/2];
  Draw_3dAxisCross(XLim, YLim, ZLim,PlotOpt);
  % Plot_Model(Model,Ufin,PlotOpt);
  Plot_DeformedStructure(Model,ElemData,State.U(1:Model.nf),[],PlotOpt);

end


%%
function [Model,ElemData] = CreateColumn(ElemName, ElemOpts, ProbOpts, BOUN)

  L   = ProbOpts.L;
  nen = ProbOpts.nen;
  ne  = ProbOpts.ne;
  Q = eye(3);

  if contains(ElemName, "1D")
    ne = 20;
  end

  % Note: displacement elements should use reduced integration
  if contains(ElemName, "Shear")
    ElemOpts.nIP = nen - 1;
  end

  ElemData = cell(ne);
  for el=1:ne, ElemData{el} = ElemOpts; end


  if contains(ElemName, "1D")
    ndm = 1;
    ndf = 1;
    nen = 3;

  elseif contains(ElemName, "2d")
    ndm = 2;
    ndf = 3;
    dofs = [1 2 6];

  else
    ndm = 3;
    ndf = 6;
    dofs = 1:6;

    for el=1:ne
      % ElemData{el}.nIP   = ElemOpts.nIP;
      ElemData{el}.Iz    = ElemOpts.I;
      ElemData{el}.Iy    = ElemOpts.I*10;
      ElemData{el}.yornt = Q*[0;1;0];
      ElemData{el}.GeomData.yornt = Q*[0;1;0];
      ElemData{el}.JntOff = zeros(3,2);
    end
  end

  nn = ne*(nen-1)+1;
  XYZ = [ linspace(0,L,nn)'  zeros(nn,1)  zeros(nn,1)]*Q';

  BOU = zeros(size(XYZ,1), ndf);
  BOU(1,:)   = BOUN(1,dofs);
  BOU(end,:) = BOUN(end,dofs);

  CON = zeros(ne,nen);
  for i=1:ne, CON(i,:) = 1 + (i-1)*(nen-1):(i-1)*(nen-1)+nen; end
  CON = num2cell(CON, 2);

  [EName{1:ne}] = deal(ElemName);

  Model = Create_Model(XYZ(:,1:ndm), CON, BOU, EName);

  ElemData = Structure('chec', Model, ElemData);
end

function [State, Post, SolStrat] = Analyze(Model, ElemData, P, nostep)
    Mmax = 0;
    angle = 0; %pi/2;
    nn = size(Model.XYZ,1);


    %% Loading
    switch Model.ndm
      case 1
        Pmax = P;
        Pe = zeros(nn,1);
        nostep = 30;
        % Force is applied as body force for Elastica element.
        for el=1:length(ElemData)
          ElemData{el}.Q   =  sin(angle)*Pmax/nostep;
          ElemData{el}.P   = -cos(angle)*Pmax/nostep;
          % EData{el}.ang = angle;
        end

      case 2
        Pe = zeros(nn,3);
%       Pe(end, 1:3) = ExpSO3([0 0 angle])*[Pmax; 0; Mmax];
        Pe(end, 1:3) = P;

      case 3
        Q = eye(3);
        Pe = zeros(nn,6);
        Pe(end, 1:6) = P; % + 1e-4;% Q*ExpSO3([0 0 angle])*[Pmax; 0; 0];
%       Pe(end, 1:3) = Q*ExpSO3([0 0 angle])*[Pmax; 0; 0];
%       Pe(end, 4:6) = M + 1e-8; % Q*[0; 0; Mmax];
%       Pe(end, 4:6) = Q*[0; Mmax/nostep; 0];
       %Pe(end, 4:6) =   [0; 0; Mmax/nostep];

    end
    Loading = Create_Loading(Model, Pe/nostep);
      

    SolStrat = Initialize_SolStrat;
    SolStrat.IncrStrat.Dlam0   = 1/nostep;
    SolStrat.IncrStrat.LFCtrl  = 'no';
    SolStrat.IterStrat.LFCtrl  = 'no';
    SolStrat.IterStrat.maxiter = 40;
    SolStrat.IterStrat.tol = 1e-12;
    SolStrat.Debug  = false;
    SolStrat.Output = 0;
    SolStrat.PUHist = false;


    [State, Post, SolStrat] = MultiStep(Model,ElemData,Loading,nostep,SolStrat);

end

