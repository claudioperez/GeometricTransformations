function [Geom, Local] = Transform_IsoPull(CoroData, xyz, ElemState, Rotations)
% ==========================================================================================
%
%  [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric Transformations of Finite
%      Elements" Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
% ==========================================================================================
% function by Claudio M. Perez                                                          2023
% ------------------------------------------------------------------------------------------
%
ndm = 3;
% number of element nodes
nn = size(xyz,2);
%                        ndf
u = ExtrReshu(ElemState,  6,   nn);

Pres = ElemState.Pres;
Iter = ElemState.Pres;
Incr = ElemState.Past;

[Ln, e1]   = ElmLenOr(xyz(:, [1 end]) + u(1:ndm,[1 end]));

[Tr, AuxMat] = Tran3dFrm_IsoPull(Rotations, e1, CoroData.IsoName, ...
                                 ElemState.Past.Geom.CorotTriad); %, incr);

Geom.CorotLength  = Ln;
Geom.CorotTriad   = Tr;
Geom.AuxMat       = AuxMat;
Geom.nn           = nn;

Local.lamda = ElemState.lamda;
Local.u     = zeros(1, nn*6);
Local.Du    = zeros(1, nn*6);
Local.DDu   = zeros(1, nn*6);
Local.xyz   = zeros(3, nn);

for i=1:nn
  Local.Triad(i).Rmat = Tr'*Rotations(:,:,i);
end

Ga = Tran3dFrm_IsoPush(CoroData, Geom);

L  = ElmLenOr(xyz);
ic = floor(0.5*(nn + 1));
Qi = repmat({Tr},nn*2,1);
Q  = blkdiag(Qi{:});
w  = Ga*Q'*ElemState.DDu;

switch CoroData.Translate
  case 1
    Local.xyz(1,:) = linspace(0, L, nn) - L/(nn-1)*(ic-1);
    c = xyz(:,ic) + u(1:ndm,ic);
    for i=1:nn
      Local.u (6*(i-1)+1:6*i-3) = Tr'*(xyz(:,i) + u(1:ndm,i) - c)  - Local.xyz(:,i) ;%- u(1:ndm,1);
      Local.Du(6*(i-1)+1:6*i-3) = Tr'*(ElemState.Du(6*(i-1)+1:6*i-3) );
    end

  case 2
    Local.xyz = Tr'*xyz;
    for i=1:nn
      Local.u (6*(i-1)+1:6*i-3) = Tr'*u(1:ndm,i);
      Local.Du(6*(i-1)+1:6*i-3) = Tr'*ElemState.Du(6*(i-1)+1:6*i-3);
    end

  case 4
    c = xyz(:,ic) ;%+ u(1:ndm,1);
    Local.xyz = Tr'*(xyz-c);
    for i=1:nn
      v = cross(Tr*w, xyz(:,i)+u(1:ndm,i) - c); %+ ElemState.DDu(1:3);
      Local.DDu(6*(i-1)+1:6*i-3) = Tr'*(ElemState.DDu(6*(i-1)+1:6*i-3) - v);
      Local.Du(6*(i-1)+1:6*i-3)  = Pres.Local.Du(6*(i-1)+1:6*i-3) + Local.DDu(6*(i-1)+1:6*i-3);
%     Local.u(6*(i-1)+1:6*i-3)   = Pres.Local.u(6*(i-1)+1:6*i-3)  + Local.DDu(6*(i-1)+1:6*i-3);
      Local.u (6*(i-1)+1:6*i-3)  = Tr'*(xyz(:,i) + u(1:ndm,i) - c)  - Local.xyz(:,i) ; %+ u(1:ndm,1);
    end

  case {3, 5, 6}
%   if Past.Local.Du == Pres.Local.Du
%     Pres.Local.Du(:) = 0;
%   end
    c  = xyz(:,ic) + u(1:ndm,ic);
    dL = Geom.CorotLength - L;
    Local.xyz(1,:) = linspace(0, L, nn) - L/(nn-1)*(ic-1);
%   Local.xyz = Tr'*(xyz - c );
    for i=1:nn
      v = cross(Tr*w, xyz(:,i)+u(1:ndm,i) - c) + ElemState.DDu(6*(ic-1)+1:6*ic-3);
%     Local.xyz(1:ndm, i) = Local.xyz(1:ndm, i) - Pres.Local.u(6*(i-1)+1:6*i-3)';

      Local.DDu(6*(i-1)+1:6*i-3) = Tr'*(ElemState.DDu(6*(i-1)+1:6*i-3) - v);
      Local.Du(6*(i-1)+1:6*i-3) = (Tr'*Pres.Geom.CorotTriad)*Pres.Local.Du(6*(i-1)+1:6*i-3)' + Local.DDu(6*(i-1)+1:6*i-3)';
      switch CoroData.Translate
        case 3
          Local.u(6*(i-1)+1:6*i-3)  = Pres.Local.u(6*(i-1)+1:6*i-3)  + Local.DDu(6*(i-1)+1:6*i-3);
        case 6
          Local.u(6*(i-1)+1:6*i-3)  = (Tr'*Pres.Geom.CorotTriad)*Pres.Local.u(6*(i-1)+1:6*i-3)'  + Local.DDu(6*(i-1)+1:6*i-3)';
        case 5
          Local.u(6*(i-1)+1:6*i-3)  = [dL*(i-1)/(nn-1) 0 0]';
      end
    end
end

Geom.xyz = Local.xyz;

