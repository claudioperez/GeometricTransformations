%
%  [1] Jelenić G, Crisfield MA (1999) Geometrically exact 3D beam theory:
%      implementation of a strain-invariant finite element for statics and
%      dynamics. Computer Methods in Applied Mechanics and Engineering,
%      171(1–2):141–171. https://doi.org/10/dj37b3
%
% clear all
CleanStart
L  = 1;
nn = 3;

Data.A  = 1;
Data.E  = 1;
Data.G  = 1;
Data.J  = 1;
Data.Iy = 1;
Data.Iz = 1;
Data.yornt = [0 1 0]';
Data.nodix = 1:nn;
Data.ElemName = 'DisplShear3dFrm_wCS';
Data.CoroData.Translate = 4;
Data.CoroData.IsoName = 'SFIN';

% R = ExpSO3([  0.2  1.2 -0.50]);

% xyz = [ linspace(0,L,nn)'  zeros(nn,1)  zeros(nn,1)]';

Formulations = {'Iter', 'Incr', 'Init'};

if false
for i=1:2
% clear State
% State.u(  4:6  )   = [  1.0 -0.5  0.25]';
% State.u(6+4:6+6  ) = [ -0.4  0.7  0.10]';
  Data.Update = Formulations{i};
  [Resp01, Resp02] = Jelenic('DisplShear3dFrm_wCS', Data, xyz);

  Kappa01 = (Resp01.Pres.Rotations{1}'*Resp01.Pres.Curvature(:,1))';

  Kappa02 = (Resp02.Pres.Rotations{1}'*Resp02.Pres.Curvature(:,1))';
  
  fprintf("      %f %f %f %f %f %f\n", Kappa01, Kappa02);
end
end


Triad = {'None', 'C3', 'SFIN'};
for j = 1:3
for nn = 2:4
  Data.Update = Formulations{2};
  Data.CoroData.RotName = Data.Update;
  Data.CoroData.IsoName = Triad{j};
  Data.nodix = 1:nn;
  xyz = [ linspace(0,L,nn)'  zeros(nn,1)  zeros(nn,1)]';

  [Resp01, Resp02] = Jelenic('GeomWrap3dFrm', Data, xyz);

  Kappa01 = (Resp01.Pres.Elem.Rotations{1}'*Resp01.Pres.Elem.Curvature(:,1))';

  Kappa02 = (Resp02.Pres.Elem.Rotations{1}'*Resp02.Pres.Elem.Curvature(:,1))';
  
  fprintf("%5s %2i %f %f %f %f %f %f\n", Triad{j}, nn, Kappa01, Kappa02);

end
end


function [Resp01, Resp02] = Jelenic(Elem, Data, xyz)
  State = [];
  State.u = zeros(1,length(Data.nodix)*6);

  State.u(  4:6  )   = [  1.0 -0.5  0.25]';
  State.u(end-2:end) = [ -0.4  0.7  0.10]';
  State.Du  = State.u';
  State.DDu = State.u';
  State.lamda = 1;

% Data.Update = Formulations{i};
  Data  = feval(Elem, 'chec', [], xyz, Data, State);
  State = feval(Elem, 'init', [], xyz, Data, State);

  Resp01 = feval(Elem, 'forc', [], xyz, Data, State);
% Kappa01 = (Resp.Pres.Rotations{1}'*Resp.Pres.Curvature)';

%%
  Ar = ExpSO3([0.2 1.2 -0.5]);

% clear State
  State.u(  4:6  )   = LogSO3(Ar*ExpSO3(State.u(4:6)));     %[ 1.0015 0.3468 -0.8372 ];
  if length(Data.nodix) > 2
    State.u(10:12)     = LogSO3(Ar*ExpSO3(State.u(10:12)));
  end
  if length(Data.nodix) > 3
    State.u(16:18)     = LogSO3(Ar*ExpSO3(State.u(16:18)));
  end
  State.u(end-2:end) = LogSO3(Ar*ExpSO3(State.u(end-2:end))); %[ 0.0885 1.9332 -0.0819 ];
  State.Du  = State.u';
  State.DDu = State.u';
  State.lamda = 1;
  State = feval(Elem, 'init', [], xyz, Data, State);
  Resp02  = feval(Elem, 'forc', [], xyz, Data, State);
% Kappa02 = (Resp.Pres.Rotations{1}'*Resp.Pres.Curvature)';

end
