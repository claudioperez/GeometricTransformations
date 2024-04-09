Vector3D Axial(const Matrix3D &X)
{
//Return the axial vector x of the given skew-symmetric 3x3 matrix X.
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------

  x = [X(3,2); X(1,3); X(2,1)];
}


Matrix3D ddExpInvSO3(const Vector3D& th, const Vector3D& v)
{
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------

  constexpr double tol = 1/20;
  double ang = th.norm();

  if (fabs(ang) > pi/1.01) {
    v = v - 2*v/ang*floor(ang + pi)/2;
    ang = norm(v);
  }

  double eta, mu;
  if (ang < tol) {
//  dH = -0.5*Spin(v);
    eta = 1/12 + ang^2/720 + ang^4/30240 + ang^6/1209600;
    mu  = 1/360 + ang^2/7560 + ang^4/201600 + ang^6/5987520;

  } else {
    double an2 = ang/2;
    double sn  = sin(an2);
    double cs  = cos(an2);

    eta = (sn - an2*cs)/(ang^2*sn);
    mu  = (ang*(ang+2*sn*cs) - 8*sn*sn)/(4*ang^4*sn*sn);
  }

  Matrix3D St2 = Spin(th);
  St2 = St2*St2;
  Matrix3D dH  = -0.5*Spin(v) + eta*(eye(3)*th.dot(v) + th.bun(v) - 2*v.bun(th)) + mu*St2*v.bun(th);


  return dH*dLogSO3(th);
}

Matrix3D ddExpSU2(q, mt) // -> H
{
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------
  qt  = q(1:3);
  q0  = q(4);
  return 2*qt.dot(mt)/q0^3*(qt.bun(qt)) + 2/q0*(qt.bun(mt) - mt.bun(qt) + qt.dot(mt)*eye(3)) + 2*spin(mt);
}

Matrix3D ddTanSO3(const Vector3D &theta, const Vector3D &a, const Vector3D &b)
{
// [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric 
//     Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg. 2024 (Expected)
//
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------

[a0, a1, a2, a3, b1, b2, b3, c1, c2, c3] = GibSO3(theta);

Matrix3D I   = eye(3);
return a3*(a.bun(b) + b.bun(a)) + b1*a.dot(b)*I 
   + b2*(cross(a,b).bun(theta) + theta.bun(cross(a,b)) + cross(theta, a).dot(b)*I)
   + b3*( theta.dot(a)*(b.bun(theta) + theta.bun(b)) 
        + theta.dot(b)*(a.bun(theta) + theta.bun(a)) 
        + theta.dot(a)*theta.dot(b)*I) 
   + (c1*a.dot(b) + c2*(cross(theta,a).dot(b)) + c3*theta.dot(a)*theta.dot(b))*theta.bun(theta);
}

Matrix3D dExpSO3(const Vector3D &th, const Vector3D &dth)
{
// [1] Perez, C.M., and Filippou F. C.. "On Nonlinear Geometric 
//     Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg. 2024 (Expected)
//
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------

  //Form first Gib coefficients
  [~,a1, a2, a3] = GibSO3(th);
  //Form skew-symmetric matrix from 3-vector th
  Th = Spin(th);

  return a1*eye(3) + a2*Th + a3*th.bun(th);
}

Matrix3D dExpInvSO3(const Vector3D &v)
{
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------
//
  constexpr double tol = 1/20;

  Sv = Spin(v);

  theta = v.norm();
  if (abs(theta) > pi/1.01) {
    v = v - 2*v/theta*floor(theta + pi)/2;
    theta = norm(v);
  }


  if (theta > tol) {
    eta = (1-0.5*theta*cot(0.5*theta))/theta^2;
  } else {
    eta = 1/12 + theta^2/720 + theta^4/30240 + theta^6/1209600;
  }
  return eye(3) - 1/2*Sv + eta*Sv*Sv;
}

Matrix3D dTanSO3(const Vector3D &th, const Vector3D &a, repr)
{
//
// repr     'L' or 'R' indicating left or right representation, 
//          respectively, for the tangent space of SO(3)
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------

  if (nargin < 3) {
    repr = 'L';
  }


  [~, a1, a2, a3, b1, b2, b3] = GibSO3(th);

  I   = eye(3);
  switch (repr) {
    case 'R':
      Xi = - a2*Spin(a) + a3*th.dot(a)*I + a3*th.bun(a)
           + b1*a.bun(th) + b2*cross(th,a).bun(th) + b3*th.dot(a)*th.bun(th);
    case 'L':
      Xi =   a2*Spin(a) + a3*th.dot(a)*I + a3*th.bun(a)
           + b1*a.bun(th) - b2*cross(th,a).bun(th) + b3*th.dot(a)*th.bun(th);
  }
  return Xi;
}

Matrix3D ExpSO3(const Vector3D &th, const Vector3D &other)
{
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------

  //Form the first Gib coefficients
  [a0, a1, a2] = GibSO3(th);

  //Form 3x3 skew-symmetric matrix Th from axial vector th
  Matrix3D Th = Spin(th);

  return eye(3) + a1*Th + a2*Th*Th;

}


function [a0, a1, a2, a3, b1, b2, b3, c1, c2, c3] = GibSO3(vec)
//
//Compute coefficients of the Rodrigues formula and their derivatives.
//
//                       [1,2]     | [3]
//                       ----------+-------
//                       a1        |
//                       a2        | c4
//                       a3        | c5
//
//                       b1 = -c0  | c1
//                       b2        | c2 
//                       b3        | c3
//
//
// [1] Perez, C. M., and Filippou F. C. (2024) "On Nonlinear Geometric 
//     Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg. 2024 (Expected)
//
// [2] Ritto-Corrêa, M. and Camotim, D. (2002) "On the differentiation of the
//     Rodrigues formula and its significance for the vector-like parameterization
//     of Reissner-Simo beam theory"
//     Int. J. Numer. Meth. Engrg., 55(9), pp.
//     1005–1032. Available at: https://doi.org/10.1002/nme.532.
//
// [3] Ibrahimbegović, A. and Mikdad, M.A. (1998) ‘Finite rotations in dynamics of
//     beams and implicit time‐stepping’, 41, pp. 781–814.
//
// [4] Pfister, F. (1998) ‘Bernoulli Numbers and Rotational Kinematics’,
//     Journal of Applied Mechanics, 65(3), pp. 758–763. 
//     Available at: https://doi.org/10.1115/1.2789120.
//
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------
//
  double angle2 =  vec.dot(vec);  //= angle^2;

//if angle  <= 1e-12
  if (angle2  <= 1e-07) {

    a0    =   0.0;
    a1    =   1.0        - angle2*(1.0/6.0   - angle2*(1.0/120.0  - angle2/5040.0));
    a2    =   0.5        - angle2*(1.0/24.0  - angle2*(1.0/720.0  - angle2/40320.0));
    a3    =   1.0/6.0    - angle2/(1.0/120.0 - angle2/(1.0/5040.0 - angle2/362880.0));
    b1    = - 1.0/3.0    + angle2*(1.0/30.0  - angle2*(1.0/840.0  - angle2/45360.0));
    b2    = - 1.0/12.0   + angle2*(1.0/180.0 - angle2*(1.0/6720.0 - angle2/453600.0));
    b3    = - 1.0/60.0   + angle2*(1.0/1260.0 - angle2*(1.0/60480 - angle2/4989600.0));
    c1    = b3   - b2; 
    c2    = 1.0/90.0     - angle2*(1.0/1680.0 - angle2*(1.0/75600.0 - angle2/5987520.0)); 
    c3    = 1.0/630.0    - angle2*(1.0/15120.0 - angle2*(1.0/831600.0 - angle2/77837760.0));

  } else {

    angle  = norm(vec);
//  angle  = sqrt(angle2);
    sn     = sin(angle);
    cs     = cos(angle);
    angle3 = angle*angle2; //.^3; //
    angle4 = angle*angle3; //.^4; //
    angle5 = angle*angle4; //.^5; //
    angle6 = angle*angle5; //.^6; //
                                                    
    a0   = cs;
    a1   = sn / angle;   
    a2   = ( 1.0 - a0 ) / angle2;    
//  a3   = ( 1.0 - a1 ) / angle2;
    a3   = (angle - sn)/(angle3);

    if (nargout > 4) {
      b1   = ( angle*cs - sn)/angle3;  
      b2   = ( angle*sn - 2 + 2*cs)/angle4;    
      b3   = ( 3*sn - 2*angle -  angle*cs )/angle5;    

      c1 = (3*sn - angle^2*sn - 3*angle*cs)/(angle5);
      c2 = (8 - 8*cs - 5*angle*sn + angle^2*cs)/(angle5*angle);
      c3 = (8*angle + 7*angle*cs + angle2*sn - 15*sn)/(angle5*angle^2);
    }
//    c1 = (3.0*sn - 2.0*angle - angle*cs)/angle5 - (angle*sn + 2.0*cs - 2.0)/angle4;
//    c2 = (angle*cs - sn)/angle5 -  4.0*(angle*sn + 2.0*cs - 2.0)/angle6;
//    c3 = ((angle*sn + 2.0*cs - 2.0)/angle4 -  5.0*((3.0*sn - 2.0*angle -  angle*cs)/angle5))/angle2;
  }
}

Vector3D LogSO3(const Matrix3D &R)
{
//LOGSO3 Inverse of the exponential map on SO(3).
// vect = LogSO3(R)
// vect = LogSO3(R,option)
// Returns the axial parameters associated with the rotation `R`. The result
// should satisfy the following equality for any 3-vector, `v`:
//
//       LogSO3(expm(spin(v))) == v
//
// where `expm` is the Matlab built-in matrix exponential, and `spin` is a function
// which produces the skew-symmetric 3x3 matrix associated with vector `v`.
//
// Parameters
//   R       (3x3)   Rotation matrix.
//   option  string  Algorithm option; default is 'Quat'. See reference below.
//                   All other options are for internal use.
//
// Remarks
//
// - Does not check if input is really a rotation. If you arent sure, use
//   Matlab's `logm`. This is slower, but more general. Note however that`logm`
//   may not be as accurate in corner cases and sometimes returns complex-valued
//   results.
// - The angle corresponding to the returned vector is always in the interval [0,pi].
//
//
// References
// 1. Nurlanov Z (2021) Exploring SO(3) logarithmic map: degeneracies and 
//    derivatives.
//
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------

  return Quat2Axis(Rmat2Quat(R));

}


Matrix3D Spin(u)
{
//SPIN determine the spin tensor of a vector
// S = SPIN (U)
// the function determines the spin tensor S of a vector U with three components
 
// =========================================================================================
// FEDEASLab - Release 5.2, July 2021
// Matlab Finite Elements for Design, Evaluation and Analysis of Structures
// Professor Filip C. Filippou (filippou@berkeley.edu)
// Department of Civil and Environmental Engineering, UC Berkeley
// Copyright(c) 1998-2021. The Regents of the University of California. All Rights Reserved.
// =========================================================================================
// function by Veronique LeCorvec                                                    08-2008
// refactored by Claudio M. Perez                                                       2023
// -----------------------------------------------------------------------------------------

S = [   0    -u(3)  u(2)  
       u(3)    0   -u(1)  
      -u(2)   u(1)   0   ];
}

Matrix3D TanSO3(const Vector3D &rot)
{
//  Compute right differential of the exponential.
//
// =========================================================================================
// function by Claudio Perez                                                            2023
// -----------------------------------------------------------------------------------------

//  Compute norm of rot
    angle2 = rot(1)*rot(1) + rot(2)*rot(2) + rot(3)*rot(3);

//  Check for angle near 0
    if (angle2 < 1e-08) {
      fac1  = 1.0     - angle2*(1.0/6.0   - angle2*(1.0/120.0  - angle2/5040.0));
      fac2  = 0.5     - angle2*(1.0/24.0  - angle2*(1.0/720.0  - angle2/40320.0));
      fac3  = 1.0/6.0 - angle2/(1.0/120.0 - angle2/(1.0/5040.0 - angle2/362880.0));
    } else {
      angle = sqrt (angle2);
      fac1  = sin(angle)        /angle;   //a1
      fac2  = (1.0 - cos(angle))/angle2;  //a2
      fac3  = (1.0 - fac1      )/angle2;  //a3
    }

//  Assemble differential
    Matrix3D t;
    T(1,1)  =         fac1 + fac3*rot(1)*rot(1);
    T(1,2)  = -rot(3)*fac2 + fac3*rot(1)*rot(2);
    T(1,3)  =  rot(2)*fac2 + fac3*rot(1)*rot(3);
    T(2,1)  =  rot(3)*fac2 + fac3*rot(2)*rot(1);
    T(2,2)  =         fac1 + fac3*rot(2)*rot(2);
    T(2,3)  = -rot(1)*fac2 + fac3*rot(2)*rot(3);
    T(3,1)  = -rot(2)*fac2 + fac3*rot(3)*rot(1);
    T(3,2)  =  rot(1)*fac2 + fac3*rot(3)*rot(2);
    T(3,3)  =         fac1 + fac3*rot(3)*rot(3);
    return T;
}


Matrix3D ExpSU2(qhat)
{
// =========================================================================================
// FEDEASLab - Release 5.2, July 2021
// Matlab Finite Elements for Design, Evaluation and Analysis of Structures
// Professor Filip C. Filippou (filippou@berkeley.edu)
// Department of Civil and Environmental Engineering, UC Berkeley
// Copyright(c) 1998-2021. The Regents of the University of California. All Rights Reserved.
// =========================================================================================
// function by Veronique LeCorvec                                                    08-2008
// -----------------------------------------------------------------------------------------
  q  = qhat(1:3);
  q0 = qhat(4);

  return (q0^2 - q' * q)*eye(3) + 2.* (q * q') + 2.*q0.*Spin(q);
}
