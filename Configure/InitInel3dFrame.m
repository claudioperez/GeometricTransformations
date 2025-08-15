function ElemResp = InitInel3dFrame(ElemData, ElemType, ElemState)

ndm = ElemType.ndm;
nsr = ElemType.nsr;
nqv = ElemType.nqv;
ndf = ElemType.ndf;
nIP = ElemData.nIP;

if nargin < 3
  ElemState = [];
end

%
%% Section
%
for i=1:nIP
  % initialize secction model
  SecState = feval (ElemData.SecName,'init',i,ndm,ElemData.SecData{i});
  % section and material history variables
  ElemState.Pres.Sec{i} = SecState.Pres;
  % variables for post-processing
  ElemState.Pres.Sec{i}.e   = zeros(nsr,1);
  ElemState.Pres.Sec{i}.De  = zeros(nsr,1);
  ElemState.Pres.Sec{i}.DDe = zeros(nsr,1);
end
%
%% Basic system
%
ElemState.Pres.q = zeros(nqv,1);
%
ElemResp = ElemState;



