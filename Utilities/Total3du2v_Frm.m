function [v,vthetaI,vthetaJ,RefTriads] = Total3du2v_Frm (xyz,GeomData,u)
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------


% Initialize any history variables that may be used by the underlying implementation.
GeomState = Corot3dFrm_Initialize(xyz, GeomData);

u = reshape(u,1,[])';
ElemState.u   = u;
ElemState.Du  = u;
ElemState.DDu = u;

% Form updated nodal triads
RefTriads = Corot3dFrm_Update(xyz, GeomData, GeomState, ElemState);

% Form the element deformations
[v,vthetaI,vthetaJ] = Corot3dFrm_Deformation(GeomData, RefTriads);

end
