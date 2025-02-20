function [ag, kg] = CorotRDS_Tangent(CoroState, GeomData, q, v)
% ==========================================================================================
% The following matrices transform the local element variables
% (q,v) to the form used by Crisfield (1990) and Remo's original
% implementation.
%
%            |       |   |       |
%          N  Miz Mjz  T  Miy Mjy
    ac = [ 0   0   0  -1   0   0 ;  % T
           0   1   0   0   0   0 ;
           0   0   0   0  -1   0 ;
           0   0   0   1   0   0 ;  % T
           0   0   1   0   0   0 ;
           0   0   0   0   0  -1 ;  % 
           1   0   0   0   0   0 ]; % N

%            |           |           |
%          N  Ti  Mi1 Mi2 Tj  Mj1 Mj2
    P  = [ 0   1   0   0   0   0   0 ;  %  T
           0   0   0   1   0   0   0 ;  %  Mi2
           0   0  -1   0   0   0   0 ;  % -Mi1
           0   0   0   0   1   0   0 ;
           0   0   0   0   0   0   1 ;
           0   0   0   0   0  -1   0 ;
           1   0   0   0   0   0   0 ];

    [~, vi, vj] = Corot3dFrm_Deformation(GeomData, RefTriads);

    ul = P*[v(1); vi; vj];

%   ul = ac*v;
    MN = ac*q;

    vtheta = ul(1:6);

    % Extract variables defining the current geometry
    Ln = RefTriads.Ln;

    e1 = CoroState.CorotTriad(1:3,1);
    e2 = CoroState.CorotTriad(1:3,2);
    e3 = CoroState.CorotTriad(1:3,3);


    t1 = RefTriads.Ti(1:3,1);
    t2 = RefTriads.Ti(1:3,2);
    t3 = RefTriads.Ti(1:3,3);

    u1 = RefTriads.Tj(1:3,1);
    u2 = RefTriads.Tj(1:3,2);
    u3 = RefTriads.Tj(1:3,3);

    r1 = RefTriads.Rbar(1:3,1);
    r2 = RefTriads.Rbar(1:3,2);
    r3 = RefTriads.Rbar(1:3,3);

    %
    % Form the static matrix
    %
    A = (1/Ln)*(eye(3) - e1*e1');

    Lr2 = getLmatrix(r2,r1,e1,A);
    Lr3 = getLmatrix(r3,r1,e1,A);

    O = zeros(3,3);

    AA= [ A
          O
         -A
          O ];

    O = zeros(3,1);
    %       d1                a                 d2   b
    h1 = [  O', (-Spin(t3)*e2 + Spin(t2)*e3)',  O',  O']';
    h2 = [  O', (-Spin(t2)*e1 + Spin(t1)*e2)',  O',  O']';
    h3 = [  O', (-Spin(t3)*e1 + Spin(t1)*e3)',  O',  O']';

    Fbar(:,1) =         Lr3*t2 - Lr2*t3 + h1;
    Fbar(:,2) = Lr2*t1  +AA*t2          + h2;
    Fbar(:,3) = Lr3*t1           +AA*t3 + h3;

    %       d1  a   d2                 b
    h4 = [  O', O',  O', (-Spin(u3)*e2 + Spin(u2)*e3)']';
    h5 = [  O', O',  O', (-Spin(u2)*e1 + Spin(u1)*e2)']';
    h6 = [  O', O',  O', (-Spin(u3)*e1 + Spin(u1)*e3)']';

    Fbar(:,4) =         Lr3*u2 - Lr2*u3 + h4;
    Fbar(:,5) = Lr2*u1  +AA*u2          + h5;
    Fbar(:,6) = Lr3*u1           +AA*u3 + h6;


    F = zeros(12,7);
    for i = 1:6
      F(:,i) = Fbar(:,i)/(2*cos(vtheta(i)));
    end
    F(:,7) = [-e1' O' e1' O']';

    ar = RefTriads.Tr';
    ar = blkdiag(ar,ar,ar,ar);
    ag = ac'*F'; %*ar;

    %
    % Geometric tangent components
    %
    if nargout > 1
      m = MN(1:6)./(2*cos(vtheta));

      % Ksigma1 -------------------------------
      % Equation C.6
      Ks1_11 = MN(7)*A;
      Ks1_33 =  Ks1_11;
      Ks1_13 = -Ks1_11;
      Ks1_31 = -Ks1_11;

      O = zeros(3);

      Ks1 = [Ks1_11  O  Ks1_13  O;
                O    O     O    O;
             Ks1_31  O  Ks1_33  O;
                O    O     O    O];

      % Ksigma3 -------------------------------
      Kbar2 = -Lr2*(m(4)*Spin(t3) + m(2)*Spin(t1)) + ...
               Lr3*(m(4)*Spin(t2) - m(3)*Spin(t1)) ;

      Kbar4 =  Lr2*(m(4)*Spin(u3) - m(5)*Spin(u1)) - ...
               Lr3*(m(4)*Spin(u2) + m(6)*Spin(u1));

      O = zeros(12,3);

      Ks3 = [O Kbar2 O Kbar4];

      % Ksigma4 -------------------------------
      Ks4_22 =  m(4)*( Spin(e2)*Spin(t3) - Spin(e3)*Spin(t2)) + ...
                m(2)*(-Spin(e1)*Spin(t2) + Spin(e2)*Spin(t1)) + ...
                m(3)*(-Spin(e1)*Spin(t3) + Spin(e3)*Spin(t1));

      Ks4_44 = -m(4)*( Spin(e2)*Spin(u3) - Spin(e3)*Spin(u2)) + ...
                m(5)*(-Spin(e1)*Spin(u2) + Spin(e2)*Spin(u1)) + ...
                m(6)*(-Spin(e1)*Spin(u3) + Spin(e3)*Spin(u1));

      O = zeros(3);
      Ks4 = [   O    O     O    O;
                O  Ks4_22  O    O;
                O    O     O    O;
                O    O     O  Ks4_44];


      % Ksigma5 -------------------------------
      Ks5_12 = -(m(2)*A*Spin(t2) + m(3)*A*Spin(t3));
      Ks5_32 = -Ks5_12;

      Ks5_14 = -(m(5)*A*Spin(u2) + m(6)*A*Spin(u3));
      Ks5_34 = -Ks5_14;

      Ks5_21 =  Ks5_12';
      Ks5_23 = -Ks5_21;

      Ks5_41 =  Ks5_14';
      Ks5_43 = -Ks5_41;

      v5 = (1/Ln)*(m(2)*t2 + m(3)*t3 + m(5)*u2 + m(6)*u3);

      Ks5_11 = A*v5*e1' + e1*v5'*A + (e1'*v5)*A;
      Ks5_33 =  Ks5_11;
      Ks5_13 = -Ks5_11;
      Ks5_31 = -Ks5_11;

      Ks5 = [Ks5_11 Ks5_12 Ks5_13 Ks5_14;
             Ks5_21   O    Ks5_23   O   ;
             Ks5_31 Ks5_32 Ks5_33 Ks5_34;
             Ks5_41   O    Ks5_43   O   ];

      % Ksigma -------------------------------

      Ks2r2t3_u3 = Ks2(r2, t3-u3, r1, e1, A, Ln);
      Ks2r3u2_t2 = Ks2(r3, u2-t2, r1, e1, A, Ln);
      Ks2r2t1    = Ks2(r2,    t1, r1, e1, A, Ln);
      Ks2r3t1    = Ks2(r3,    t1, r1, e1, A, Ln);
      Ks2r2u1    = Ks2(r2,    u1, r1, e1, A, Ln);
      Ks2r3u1    = Ks2(r3,    u1, r1, e1, A, Ln);

      F6 = F(:,1:6);

      kg = Ks1 + ...
           Ks3 + Ks3' + Ks4 + Ks5 + ...
           F6 * diag (MN(1:6).* tan(vtheta))*F6' + ...
           m(4)*(Ks2r2t3_u3 + Ks2r3u2_t2) + ...
           m(2)*Ks2r2t1 + m(3)*Ks2r3t1 + m(5)*Ks2r2u1 + m(6)*Ks2r3u1;
    end
end

%%  function getLmatrix --------------------------------------------------------------------
function L = getLmatrix (ri, r1, e1, A)

  L1  = ri'*e1 * A/2 + A*ri*(e1 + r1)'/2;
  Sri = Spin(ri);
  L2  = Sri/2 - ri'*e1*Spin(r1)/4 - Sri*e1*(e1 + r1)'/4;

  L = [L1' L2' -L1' L2']';

end
  
%%  ----- function Ks2 ---------------------------------------------------------------------
function Ksigma2 = Ks2(ri, z, r1, e1, A, Ln)
  
  U = (-1/2)*A*z*ri'*A + ri'*e1*A*z*e1'/(2*Ln)+...
      z'*(e1+r1)*A*ri*e1'/(2*Ln);
    
  K11 = U + U' + ri'*e1*(2*(e1'*z)+z'*r1)*A/(2*Ln);
  K13 = -K11;
  K31 = -K11;
  K33 =  K11;
  
  Sri = Spin(ri);
  Sr1 = Spin(r1);
    
  K12 = 0.25*(-A*z*e1'*Sri - A*ri*z'*Sr1 - z'*(e1+r1)*A*Sri);
  K14 =  K12;
  K32 = -K12;
  K34 = -K12;
    
  K21 =  K12';
  K41 =  K21;
  K23 = -K21;
  K43 = -K21;
  
  Sz = Spin(z);
  K22 = (1/8)*((-ri'*e1)*Sz*Sr1 + Sr1*z*e1'*Sri + ...
        Sri*e1*z'*Sr1 - (e1+r1)'*z*Spin(e1)*Sri + 2*Sz*Sri);
    
  K24 = K22;
  K42 = K22;
  K44 = K22;
  
  Ksigma2 = [K11 K12 K13 K14;
             K21 K22 K23 K24;
             K31 K32 K33 K34;
             K41 K42 K43 K44];

end
