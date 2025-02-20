function [Ga, La] = Tran3dFrm_IsoPush(CoroData, LocalGeom)
%
%  - Crisfield, M.A. and Jelenić, G. (1999) ‘Objectivity of strain measures in the geometrically exact three-dimensional beam theory and its finite-element implementation’, Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences, 455(1983), pp. 1125–1147. Available at: https://doi.org/10.1098/rspa.1999.0352.
%  - Nour-Omid, B. and Rankin, C.C. (1991) ‘Finite rotation analysis and consistent linearization using projectors’, Computer Methods in Applied Mechanics and Engineering, 93(3), pp. 353–384. Available at: https://doi.org/10.1016/0045-7825(91)90248-5.
%  - Crisfield, M.A. (1990) ‘A consistent co-rotational formulation for non-linear, three-dimensional, beam-elements’, Computer Methods in Applied Mechanics and Engineering, 81(2), pp. 131–150. Available at: https://doi.org/10.1016/0045-7825(90)90106-V.
%  - Pacoste, C. and Eriksson, A. (1997) ‘Beam elements in instability problems’, Computer Methods in Applied Mechanics and Engineering, 144(1–2), pp. 163–197. Available at: https://doi.org/10.1016/S0045-7825(96)01165-6.
%  - Battini, J.-M. and Pacoste, C. (2002) ‘Co-rotational beam elements with warping effects in instability problems’, Computer Methods in Applied Mechanics and Engineering, 191(17–18), pp. 1755–1789. Available at: https://doi.org/10.1016/S0045-7825(01)00352-8.
%
%  =========================================================================================
%  function by Claudio Perez                                                         02-2023
%  -----------------------------------------------------------------------------------------

% Current length
Ln = LocalGeom.CorotLength;

% Current local reference frame
Tr = LocalGeom.CorotTriad;


% Default that tends to work well
G   = [0   0    0  , 0.5 0 0, 0  0    0  , 0.5  0  0;
       0   0   1/Ln,  0  0 0, 0  0  -1/Ln,  0   0  0;
       0 -1/Ln  0  ,  0  0 0, 0 1/Ln  0  ,  0   0  0];

if nargout > 1
  La = zeros(LocalGeom.nn*6, 3);
  for i=1:LocalGeom.nn
    La(6*(i-1)+1:6*i, :) = [
     -Spin(LocalGeom.xyz(:,i) + LocalGeom.u(1:3, i))
      eye(3) 
    ];
  end
end
%   L   = [    zeros(3)    ;
%               eye(3)     ;
%           0     0      0 ;  % \
%           0     0     Ln ;  % | = -Spin([Ln 0 0])
%           0   -Ln      0 ;  % /
%               eye(3)     ];


switch CoroData.IsoName
  case 'None'
    Ga = zeros(6*LocalGeom.nn, 3)';
    return

  case {'R1', 'R2', 'KM2'}
    % Nour-Omid and Rankin (1991)

    t  = Tr'*LocalGeom.AuxMat(:,2);

    if strcmp(CoroData.IsoName, 'R1')
      n = 0.0;
    else
      n =  t(1)/t(2);
    end

    %     |             |          |              |          |
    G   = [0   0   n/Ln , 1  -n  0 , 0   0  -n/Ln , 0   0  0 ;
           0   0   1/Ln , 0   0  0 , 0   0  -1/Ln , 0   0  0 ;
           0 -1/Ln  0   , 0   0  0 , 0  1/Ln  0   , 0   0  0 ];

  case {'E1'}
    % Pacoste, C. and Eriksson, A. (1997)
    t  = Tr'*LocalGeom.AuxMat(:,1);

    n   =  t(1)/t(2);

    G = [
      0   0   n/Ln,  0.5  -n/2   0, 0   0   -n/Ln,  0.5   -n/2   0 
      0   0   1/Ln,   0     0    0, 0   0   -1/Ln,   0      0    0 
      0 -1/Ln   0 ,   0     0    0, 0  1/Ln   0  ,   0      0    0
    ];

  case {'B1'}
    % Battini and Pacoste (2002)
    %
    %     |            |         |            |         |
    G   = [0   0    0    0.5 0 0   0  0    0    0.5  0  0;
           0   0   1/Ln   0  0 0   0  0  -1/Ln   0   0  0;
           0 -1/Ln  0     0  0 0   0 1/Ln  0     0   0  0];

  case {'B2'}
    % Battini, J.-M. and Pacoste, C. (2002)
    % 
    t  = Tr'*LocalGeom.AuxMat(:,1);
    t1 = Tr'*LocalGeom.AuxMat(:,2);
    t2 = Tr'*LocalGeom.AuxMat(:,3);

    n   =  t(1)/t(2);
    n11 = t1(1)/t(2);
    n12 = t1(2)/t(2);
    n21 = t2(1)/t(2);
    n22 = t2(2)/t(2);

    G = [
      0   0   n/Ln, n12/2  -n11/2  0, 0   0   -n/Ln, n22/2 -n21/2  0 
      0   0   1/Ln,   0       0    0, 0   0   -1/Ln,   0      0    0 
      0 -1/Ln   0 ,   0       0    0, 0  1/Ln   0  ,   0      0    0
    ];


  case {'SFIN'}
    % Perez and Filippou (2024), Crisfield and Jelenić (1999)
    Phi = LocalGeom.AuxMat;
    phi = norm(Phi(:,1));
    G(:) = 0;

    nn = LocalGeom.nn;
    I = floor(0.5*(nn + 1));
    J = floor(0.5*(nn + 2));

    Ga = zeros(3, LocalGeom.nn*6);
    if phi > 1e-10
      c = 1/phi*tan(0.25*phi);
      Ga(1:3, 6*(I-1)+4:6*I) = (0.5*(eye(3) + c*Spin(Phi(:,1))));
      Ga(1:3, 6*(J-1)+4:6*J) = (0.5*(eye(3) + c*Spin(Phi(:,2))));
    else
      Ga(1:3, 6*(I-1)+4:6*I) = (0.5*(eye(3) + 0.25*Spin(Phi(:,1))));
      Ga(1:3, 6*(J-1)+4:6*J) = (0.5*(eye(3) + 0.25*Spin(Phi(:,2))));
    end
    return

  case { 'C3', 'C4'} % {'C1', 'C2', 'lerp'}
    % Crisfield, M.A. (1990)
  
    Tr  = LocalGeom.CorotTriad;
    if strcmp(CoroData.IsoName, 'C3')
      E   = eye(3);
      R   = Tr'*LocalGeom.AuxMat;
    else
      E  = Tr;
      R  = LocalGeom.AuxMat;
    end

    e1  = E(:,1);
    e2  = E(:,2);
    e3  = E(:,3);

    r1  = R(:,1);
    r2  = R(:,2);
    r3  = R(:,3);

    % truss effect
    at  = (1/Ln)*(eye(3) - e1*e1');

    de3 = crisfield_variation(r3,r1,e1,at);
    de2 = crisfield_variation(r2,r1,e1,at);


    G   =  [-e2'*de3;
             e1'*de3;
            -e1'*de2];

end % switch

Ga = zeros(3, LocalGeom.nn*6);
Ga(:, 1:6) = G(:, 1:6);
Ga(:, end-6:end) = G(:, end-6:end);

for i=1:LocalGeom.nn
  switch i
    case 1
      Ga(:,6*(i-1)+1:6*i) = G(:,1:6);
    case LocalGeom.nn
      Ga(:,6*(i-1)+1:6*i) = G(:,end-5:end);
  end
end

end % function

function L = crisfield_variation(ri, r1, e1, A)
%
% Equations 51 and 52 from [1]
%
% [1] Crisfield, M.A. (1990) ‘A consistent co-rotational formulation for
%     non-linear, three-dimensional, beam-elements’, Computer Methods in Applied
%     Mechanics and Engineering, 81(2), pp. 131–150. Available at:
%     https://doi.org/10.1016/0045-7825(90)90106-V.
%
  L1  = ri'*e1 * A/2 + A*ri*(e1 + r1)'/2;
  Sri = Spin(ri);
  L2  = Sri/2 - (ri'*e1)*Spin(r1)/4 - Sri*e1*(e1 + r1)'/4;

  L = [L1' L2' -L1' L2'];

end
