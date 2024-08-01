function [a0, a1, a2, a3, b1, b2, b3, c1, c2, c3] = GibSO3(vec)
%
% Compute coefficients of the Rodrigues formula and their derivatives.
%
%                        [1,2]     | [3]
%                        ----------+-------
%                        a1        |
%                        a2        | c4
%                        a3        | c5
%
%                        b1 = -c0  | c1
%                        b2        | c2 
%                        b3        | c3
%
%
%  [1] Perez, C. M., and Filippou F. C. (2024) "On Nonlinear Geometric 
%      Transformations of Finite Elements" 
%      Int. J. Numer. Meth. Engrg. 2024 (Expected)
%
%  [2] Ritto-Corrêa, M. and Camotim, D. (2002) "On the differentiation of the
%      Rodrigues formula and its significance for the vector-like parameterization
%      of Reissner-Simo beam theory"
%      Int. J. Numer. Meth. Engrg., 55(9), pp.
%      1005–1032. Available at: https://doi.org/10.1002/nme.532.
%
%  [3] Ibrahimbegović, A. and Mikdad, M.A. (1998) ‘Finite rotations in dynamics of
%      beams and implicit time‐stepping’, 41, pp. 781–814.
%
%  [4] Pfister, F. (1998) ‘Bernoulli Numbers and Rotational Kinematics’,
%      Journal of Applied Mechanics, 65(3), pp. 758–763. 
%      Available at: https://doi.org/10.1115/1.2789120.
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------
%
  angle2 =  vec'*vec;  % = angle^2;

% if angle  <= 1e-12
  if angle2  <= 1e-07 % 1e-07

    a0    =   0.0;
    a1    =   1.0       - angle2*(1.0/6.0   - angle2*(1.0/120.0  - angle2/5040.0));
    a2    =   0.5       - angle2*(1.0/24.0  - angle2*(1.0/720.0  - angle2/40320.0));
    a3    =   1.0/6.0   - angle2/(1.0/120.0 - angle2/(1.0/5040.0 - angle2/362880.0));
    b1    = - 1.0/3.0   + angle2*(1.0/30.0  - angle2*(1.0/840.0  - angle2/45360.0));
    b2    = - 1.0/12.0  + angle2*(1.0/180.0 - angle2*(1.0/6720.0 - angle2/453600.0));
    b3    = - 1.0/60.0  + angle2*(1.0/1260.0 - angle2*(1.0/60480 - angle2/4989600.0));
    c1    = b3   - b2; 
    c2    = 1.0/90.0    - angle2*(1.0/1680.0 - angle2*(1.0/75600.0 - angle2/5987520.0)); 
    c3    = 1.0/630.0   - angle2*(1.0/15120.0 - angle2*(1.0/831600.0 - angle2/77837760.0));

  else

    angle  = norm(vec);
%   angle  = sqrt(angle2);
    sn     = sin(angle);
    cs     = cos(angle);
    angle3 = angle*angle2; %.^3; %
    angle4 = angle*angle3; %.^4; %
    angle5 = angle*angle4; %.^5; %
    angle6 = angle*angle5; %.^6; %
                                                    
    a0   = cs;
    a1   = sn / angle;   
    a2   = ( 1.0 - a0 ) / angle2;    
%   a3   = ( 1.0 - a1 ) / angle2;
    a3   = (angle - sn)/(angle3);

    if nargout > 4
      b1   = ( angle*cs - sn)/angle3;
      b2   = ( angle*sn - 2 + 2*cs)/angle4;    
      b3   = ( 3*sn - 2*angle -  angle*cs )/angle5;    

      c1 = (3*sn - angle^2*sn - 3*angle*cs)/(angle5);
      c2 = (8 - 8*cs - 5*angle*sn + angle^2*cs)/(angle5*angle);
      c3 = (8*angle + 7*angle*cs + angle2*sn - 15*sn)/(angle5*angle^2);
    end
%     c1 = (3.0*sn - 2.0*angle - angle*cs)/angle5 - (angle*sn + 2.0*cs - 2.0)/angle4;
%     c2 = (angle*cs - sn)/angle5 -  4.0*(angle*sn + 2.0*cs - 2.0)/angle6;
%     c3 = ((angle*sn + 2.0*cs - 2.0)/angle4 -  5.0*((3.0*sn - 2.0*angle -  angle*cs)/angle5))/angle2;
  end
end
