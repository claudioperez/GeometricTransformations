function Resp = Check3dFrm(el_no, xyz, ElemData, ElemType)
% Check element data for 3D frames. Assign default values, if necessary.
% 
% This function is used to allow all 3d Frame element to accept the same
% configuration data. This can be in one of three forms:
% 
% 1) Elastic section properties may be specified directly in ElemData
% 2) ...
% 3) ...
%
% -----------------------------------------------------------------------------
%
% General
%        yornt   = local y-axis orientation in global reference system (column vector)
%        Geom    = character variable for geometric transformation of node variables
%                  (linear, PDelta or corotational) (default=linear)
%        w       = uniform element load ( w(1) = longitudinal, w(2),w(3) = transverse in y and z, resp.)
%        rho     = mass density
%        JntOff  = rigid joint offsets in global X, Y, Z at element ends;
%                  column 1 for node i, column 2 for node j
%
% Section I
%        G       = shear modulus
%        J       = polar moment of inertia
%
% Section II
%         E       = modulus of elasticity
%         A       = cross sectional area
%         Iy      = moment of inertia about y-axis
%         Iz      = moment of inertia about z-axis
%         e0      = initial deformations ( e(1) = axial strain, e(2),e(3) = curvature about y and z, resp.)
%         Release = axial, torsional and end flexural releases in column vector (0=cont,1=hinge) (default=[0;0;0;0;0;0])
%
% Section III
%        nIP     = number of integration points
%        IntTyp  = function name for element integration
%
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

%
%% Common frame element parameters
%
  ndm = ElemType.ndm;
  ndf = ElemType.ndf;
  nsr = ElemType.nsr;
  nqv = ElemType.nqv;

  nen = size(xyz, 2);

%
%% General
%
  if ~isfield(ElemData,'w'),      ElemData.w      = zeros(6,1); end
  if ~isfield(ElemData,'LdId'),   ElemData.LdId   = zeros(6,1); end

%
%% Geometry
%

  if ~isfield(ElemData,'Geom'),    ElemData.Geom    = 'linear';     end
  if ~isfield(ElemData,'JntOff'),  ElemData.JntOff  = zeros(ndm,2); end
  if ~isfield(ElemData,'Release'), ElemData.Release = zeros(6,1);   end

  if ndm > 2
      if ~isfield(ElemData,'yornt')
        ElemData.yornt   = [0;1;0]; 
        warning('off','backtrace');
        warning('E:W',[
          '>> Element ', num2str(el_no,'%i') ...
          ': y-axis orientation for the section "yornt" missing in "ElemData". Default vector [0;1;0] used.'
        ]);
        warning('on','backtrace');
      end
      if isrow(ElemData.yornt), ElemData.yornt = ElemData.yornt'; end
      ElemData.yornt = Check3dFrmAxes(el_no, ElemData.yornt, xyz(:, [1 end]) + ElemData.JntOff);
  end

  if ~isfield(ElemData,'GeomData')
    if ndm > 2
       ElemData.GeomData.yornt  = ElemData.yornt;
    end
    ElemData.GeomData.JntOff = ElemData.JntOff;
  end

  L = ElmLenOr(xyz(:, [1 end]) + ElemData.JntOff);
  ElemData.L = L;
%
%% Geometry
%
  if ndm > 2
    ElemData.CoroData.yornt = ElemData.yornt;
  end
  if ~isfield(ElemData.CoroData, 'LogOpt')
    ElemData.CoroData.LogOpt = 'Exact';
  end
  if ~isfield(ElemData.CoroData, 'Quaternion')
    ElemData.CoroData.Quaternion = false;
  end
  if ~isfield(ElemData.CoroData, 'Translate')
    ElemData.CoroData.Translate = 4;
  end
  if ~isfield(ElemData.CoroData, 'Petrov')
    ElemData.CoroData.Petrov = false;
  end
  if ~isfield(ElemData.CoroData, 'IsoName')
    ElemData.CoroData.IsoName = 'B2';
  end
  if isfield(ElemType, 'iqv')
    ElemType.CoroData.iqv   = ElemType.iqv;
  end

%
%% Dynamics
%
  if  ~isfield(ElemData, 'rho')
    if isfield(ElemData, 'linear_mass')
      ElemData.rho = ElemData.linear_mass/ElemData.A;
    else
      ElemData.rho = 0.0;
    end
  end

%
%% Quadrature
%
% If an interpolated field is declared in ElemType, check quadrature
% configuration
  if isfield(ElemType, 'Field')
    if ~isfield(ElemData,'nIP'),    ElemData.nIP    = 4;         end
    if ~isfield(ElemData,'IntTyp'), ElemData.IntTyp = 'Lobatto'; end
    % obtain integration point locations and weights
    [xIP,wIP] = feval(ElemData.IntTyp,ElemData.nIP);
    nIP = length(xIP);

    % TODO: Which elements expect jacobian scaling?
    Jac = 0.5*L;
    xIP = Jac.*(1.+xIP);
    wIP = Jac.*wIP;
    if ~isfield(ElemData, 'xIP'), ElemData.xIP = xIP; end
    if ~isfield(ElemData, 'wIP'), ElemData.wIP = wIP; end
  end

%
%% Section
%
  % If there is an interpolated field, and a section has not been
  % specified, try setting up an elastic cross section

  if isfield(ElemType, 'Field') 

    % if element is uses integration but no section name or data is provided,
    if ~isfield(ElemData, 'SecData') && ~isfield(ElemData, 'SecName')

      % Assume elastic section properties were given in ElemData,
      % and use these to set up elastic cross sections.
      % If elastic section properties were not supplied in ElemData,
      % then the following calls to ElasticSection will fail.

      ElemData.SecType.Shear = ElemType.Shear;
      for i=1:ElemData.nIP
        % Run 'chec' on elastic section using ElemData
        ElemData.SecData{i}     = ElasticSection('chec', el_no, ndm, ElemData);
        ElemData.SecData{i}.ns  = nsr;
        ElemData.SecData{i}.nsr = nsr;
        ElemData.SecData{i}.isr = ElemType.isr;
        ElemData.SecData{i}.xIP = ElemData.xIP(i);
        ElemData.SecData{i}.wIP = ElemData.wIP(i);
      end
      rmfield(ElemData, 'SecType');

      ElemData.SecName = 'ElasticSection';
      ElemData.SecHndl = str2func(ElemData.SecName);

    else
%     if ~isfield(ElemData,'SecData') || ~isfield(ElemData,'SecData')
%       error('E:E',['>> Element ',num2str(el_no,'%i'),': Section data "SecData{}" missing in "ElemData".\n\n']);
%     end

      ldata = length(ElemData.SecData);
      if ldata < nIP
        for i = ldata:nIP
          ElemData.SecData{i} = ElemData.SecData{ldata};
        end
      end
      for i=1:nIP
        ElemData.SecData{i}.ns  = nsr;
        ElemData.SecData{i}.nsr = nsr;
        ElemData.SecData{i}.isr = ElemType.isr;
        ElemData.SecData{i} = feval(ElemData.SecName,'chec',i,ndm,ElemData.SecData{i});
        ElemData.SecData{i}.xIP = ElemData.xIP(i);
        ElemData.SecData{i}.wIP = ElemData.wIP(i);
      end
      ElemData.SecHndl = str2func(ElemData.SecName);
    end

  else
      % Elastic element; put elastic section properties right into ElemData.
      % Here the ElasticSection is used to validate any elastic cross section
      % parameters in ElemData. Note that the ElasticSection
      % function may not actually be used by the element for state determination.
      ElemData.SecType = ElemType;
      ElemData = ElasticSection('chec', el_no, ndm, ElemData);
      ElemData = rmfield(ElemData, 'SecType');

  end

  Resp = ElemData;
end

