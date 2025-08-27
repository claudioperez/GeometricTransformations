
function ElemMass = Mass_3dFrm02 (xyz,ElemData)
%MASS_3dFRM consistent mass matrix for 3d frame element
%  ELEMMASS = MASS_3dFRM (XYZ,ELEMDATA)
%  the function determines the lumped and consistent mass matrix
%  for a 3d frame element in global coordinates;
%  the array XYZ supplies the end node coordinates and ELEMDATA the element properties;
%  the field RHO is the mass density and the field JNTOFF of ELEMDATA are the rigid end offsets 
%  for a prismatic frame element the field A is the area;
%  for a frame element with variable cross-section the sub-field SECDATA of ELEMDATA
%  contains the section area in field A; SECDATA also contains the integration type in field
%  INTTYP and the number of integration points in field NIP;
%  the function returns ELEMMASS with the lumped mass matrix in field ml and
%  the consistent mass matrix in field mc

%  =========================================================================================
%  FEDEASLab - Release 5.4, July 2024
%  MATLAB Finite Elements for Design, Evaluation and Analysis of Structures
%  Professor Filip C. Filippou (filippou@berkeley.edu)
%  Department of Civil and Environmental Engineering, UC Berkeley
%  Copyright(c) 1998-2024. The Regents of the University of California. All Rights Reserved.
%  =========================================================================================
%  function added                                                                    05-2024

ndf = 6;          % no of element DOFs per node (2node, 2d frame element)

GeomData.JntOff = ElemData.JntOff;
GeomData.yornt  = ElemData.yornt;
rho    = ElemData.rho;

%% geometric transformation to global reference
% element length L and triad T of undeformed geometry
[L,Tr0] = DefGeom_3dFrm (xyz, GeomData);
% rotation matrix ar from global to local reference system
art = Tr0';
Mr  = blkdiag(art,art,art,art);


% lumped mass matrix
ml  = zeros(2*ndf,1);
% total mass for lumped mass determination
tm  = rho*sum(wIP.*A);
ml ([1:3 7:9]) = 0.5*Jac*tm.*ones(6,1);

% transform consistent mass matrix to global reference system
mc  = Mr'*mcb*Mr;

% return lumped and consistent element mass in ElemMass
ElemMass.ml = ml;
ElemMass.mc = mc;
end



function ShearPrismMass(ElemData)
end

function EulerPrismMass(ElemData)
end

function EulerCubicMass(ElemData)
%% mass determination
if ~isfield(ElemData,'IntTyp'), IntTyp = 'Gauss'; else, IntTyp = ElemData.IntTyp; end
IntTyp = str2func(IntTyp);
if ~isfield(ElemData,   'nIP'),    nIP = 4;       else, nIP    = ElemData.nIP;    end

if ~isfield(ElemData,'SecData')
  A = ones(nIP,1).*ElemData.A;
else
  Atmp = [ElemData.SecData{:}];
  A    = [Atmp.A]';             % column vector of section areas           
end
% polar moment of area (for circular sections it factors the torsional mass term by J/A) 
J = ElemData.J;

% obtain integration point locations and weights
[xIP,wIP] = IntTyp (nIP);
% element integration; set up displacement interpolation functions
Jac = 0.5.*L;
lp  = Lagrange (1,0,xIP);
hp  =  Hermite (3,0,xIP);

% determine consistent mass matrix in local system
mcb = zeros(2*ndf,2*ndf);
for i=1:nIP
  % displacement interpolation functions in local system
  Ai   = A(i);
  abxi = [ lp(1,i)      0          0          0          0             0        ;
             0      hp(1,i)        0          0          0         hp(2,i).*Jac ;
             0          0     -hp(1,i)        0      hp(2,i).*Jac      0        ;
             0          0          0      lp(1,i)        0             0        ];
  abxj = [ lp(2,i)      0          0          0          0             0        ;
             0      hp(3,i)        0          0          0         hp(4,i).*Jac ;
             0          0     -hp(3,i)        0      hp(4,i).*Jac      0        ;
             0          0          0      lp(2,i)        0             0        ];
  abx  = [ abxi abxj ];
  
  ms   = blkdiag( (rho*Ai).*eye(3) , rho*J/Ai );
  % determine contribution of section i to consistent mass matrix
  mcb  = mcb  + wIP(i) .* abx' * ms * abx;
end
% consistent mass matrix in local system
mcb = Jac.*mcb;

end
