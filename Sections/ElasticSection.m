function SecResp = ElasticSection(action,SecNo,ndm,SecData,SecState)
%
%          y
%          |
%      ---------
%      |   |   |
%  <---|---+   |
%  z   |       |
%      |______ |
%
%  ACTION = 'chec': check section property data for omissions and assign 
%                   default values
%           'init': initialize section history variables
%           'forc': report section resisting forces
%           'stif': report section stiffness matrix and resisting forces
%           'post': report post-processing information
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------
%
switch action
  case 'chec'
    if ~isfield(SecData, 'SecType')
      SecData.SecType = struct('Shear', false);
    end

    if ~isfield(SecData,'E'),  disp('Section ');disp(SecNo); error('Young modulus missing');end
    if ~isfield(SecData,'A'),  disp('Section ');disp(SecNo); error('area missing');end

    if ndm > 2 || SecData.SecType.Shear
      if ~isfield(SecData,'G'), disp('Section '); disp(SecNo); error('shear modulus G missing');end
    end

    if SecData.SecType.Shear
      if ~isfield(SecData,'G'), disp('Section '); disp(SecNo); error('shear modulus G missing');end
    end
    if ~isfield(SecData,'Ay')
      SecData.Ay = SecData.A;
    end

    if ndm == 2
      if ~isfield(SecData,'I')
        disp('Section ');disp(SecNo); error('moment of inertia I missing');
      else
        SecData.Iz = SecData.I;
      end

    else
      if ~isfield(SecData,'J'),    disp('Section ');disp(SecNo); error('polar moment of inertia J missing');end
      if ~isfield(SecData,'Iy'),   disp('Section ');disp(SecNo); error('moment of inertia Iy missing');end
      if ~isfield(SecData,'Iz'),   disp('Section ');disp(SecNo); error('moment of inertia Iz missing');end
      if ~isfield(SecData,'Iyz')
        SecData.Iyz = 0.0;
      end
      if ~isfield(SecData,'Az')
        SecData.Az = SecData.A;
      end
    end


    SecResp = SecData;
    return;

  case {'init', 'forc', 'stif', 'post'}
    isr  = SecData.isr;
    % Extract and form elastic constants.
    % TODO: Use corrected shear area
    EA   = SecData.E*SecData.A;
    GAy  = SecData.G*SecData.Ay;
    GAz  = SecData.G*SecData.Az;
    GA   = SecData.G*SecData.A;
    GJ   = SecData.G*SecData.J;                 % St. Venant warping constant
    EIy  = SecData.E*SecData.Iy;                % Centroidal MOI about y, (integral z^2)
    EIz  = SecData.E*SecData.Iz;                % Centroidal MOI about z
    EIyz = SecData.E*SecData.Iyz;               % Centroidal MOI
    GI0  = SecData.E*(SecData.Iy + SecData.Iy); % Centroidal polar MOI


    Ks = [ EA       0       0       0       0       0    ;
            0      GAy      0       0       0       0    ;
            0       0      GAz      0       0       0    ;
            0       0       0      GJ       0       0    ;
            0       0       0       0      EIy     EIyz  ;
            0       0       0       0      EIyz    EIz   ];

    if strcmp(action, 'init')
      SecState = [];
      e  = zeros(size(isr,2), 1);
    else
      e = SecState.e;
    end

    S = Ks(isr,isr)*e;

    SecResp = SecState;
    SecResp.s    = S;
    SecResp.e    = e;
    SecResp.ks   = Ks(isr,isr);
    SecResp.ks0  = Ks(isr,isr);
    SecResp.Past = [];
    SecResp.Pres = [];

    return;

  case {'disp'}
    SecResp = [];
end
