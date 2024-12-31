function dH = ddExpInvSO3(th, v)
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

  v   = reshape(v,  [3,1]);
  th  = reshape(th, [3,1]);

  tol = 1/20;
  ang = norm(th);
  if abs(ang) > pi/1.01
    v = v - 2*v/ang*floor(ang + pi)/2;
    ang = norm(v);
  end

  if ang < tol
%   dH = -0.5*Spin(v);
    eta = 1/12 + ang^2/720 + ang^4/30240 + ang^6/1209600;
    mu = 1/360 + ang^2/7560 + ang^4/201600 + ang^6/5987520;

  else
    an2 = ang/2;
    sn  = sin(an2);
    cs  = cos(an2);

    eta = (sn - an2*cs)/(ang^2*sn);
    mu  = (ang*(ang+2*sn*cs) - 8*sn*sn)/(4*ang^4*sn*sn);
  end

  St2 = Spin(th);
  St2 = St2*St2;
  dH  = -0.5*Spin(v) + eta*(eye(3)*(th'*v) + (th*v') - 2*(v*th')) + mu*St2*(v*th');


  dH = dH*dLogSO3(th);
end
