function [Pres,Local] = GeomTran3dFrm_Pull(CoroData, xyz, ElemState, ElemType)
% ==========================================================================================
%
%  [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric Transformations of Finite
%      Elements" Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
% ==========================================================================================
% function by Claudio M. Perez                                                          2023
% ------------------------------------------------------------------------------------------
%
Pres = ElemState.Pres;

%% (a) Update Rotation Parameters
%---------------------------------------------------------------------------
Rotations      = Transform_RotPull(CoroData, xyz, ElemState);
Pres.Rotations = Rotations;

%% (b) Isometric Transformation
%---------------------------------------------------------------------------
[Geom, Local]  = Transform_IsoPull(CoroData, xyz, ElemState, Rotations);
Pres.Geom      = Geom;

%% (c) Interpolation Parameters
%---------------------------------------------------------------------------

Iter = ElemState.Pres;
Incr = ElemState.Past;
nn   = size(xyz,2);    % number of element nodes

for i=1:nn
% Local.Triad(i).Axis = LogSO3(Local.Triad(i).Rmat);
% Local.Triad(i).Incr = LogSO3(Local.Triad(i).Rmat*Incr.Local.Triad(i).Rmat');
% Local.Triad(i).Spin = LogSO3(Local.Triad(i).Rmat*Iter.Local.Triad(i).Rmat');

% Local.u  (6*(i-1)+4:6*i)  = Local.Triad(i).Axis;
% Local.Du (6*(i-1)+4:6*i)  = Local.Triad(i).Incr;
% Local.DDu(6*(i-1)+4:6*i)  = Local.Triad(i).Spin;

  Local.u  (6*(i-1)+4:6*i)  = LogSO3(Local.Triad(i).Rmat);
  Local.Du (6*(i-1)+4:6*i)  = LogSO3(Local.Triad(i).Rmat*Incr.Local.Triad(i).Rmat');
  Local.DDu(6*(i-1)+4:6*i)  = LogSO3(Local.Triad(i).Rmat*Iter.Local.Triad(i).Rmat');
end

end

%  =========================================================================================
function ve = Deformations(xyz, CoroState, Rotations)
%
  Ln  = CoroState.CorotLength;
  L   = ElmLenOr(xyz);

  %% Compute the axial deformation
  ul = Ln - L;

  %% Compute the basic rotational deformation
  Tr = CoroState.CorotTriad;

  for i=1:size(xyz,2)
    vrot(:,i) = LogSO3(Tr'*Rotations(:,:,i));
  end

  vthetaI = vrot(:,1);
  vthetaJ = vrot(:,end);

  % Deformation vector
  ve = [ul;
        vthetaI(3);
        vthetaJ(3);
        vthetaJ(1)-vthetaI(1);
        vthetaI(2);
        vthetaJ(2)];
end

