function Rotations = Transform_RotPull(TranData, xyz, ElemState)
% ==========================================================================================
%
%  [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric Transformations of Finite
%      Elements" Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
% ==========================================================================================
% function by Claudio M. Perez                                                          2023
% ------------------------------------------------------------------------------------------
%
nn = size(xyz,2); % number of element nodes

% Rotation update type
switch TranData.RotName
  case 'Iter'
    [~, ~, incr] = ExtrReshu(ElemState, 6, nn);
    Past = ElemState.Pres;
  case 'Incr'
    [~, incr, ~] = ExtrReshu(ElemState, 6, nn);
    Past = ElemState.Past;
  case 'Init'
    [incr, ~, ~] = ExtrReshu(ElemState, 6, nn);
    Past = ElemState.Past;
  otherwise
end

R0 = Past.Rotations;
if TranData.Quaternion
  for i=1:nn, Rotations(:,:,i)  = ExpSU2(Qvec2Quat(incr(4:6,i)))*R0(:,:,i); end
else
  for i=1:nn, Rotations(:,:,i)  = ExpSO3(incr(4:6,i))*R0(:,:,i); end
end
