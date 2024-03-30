function Xi = ddExpSO3(th, a)
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------
  th = reshape(th, [3 1]);
  a  = reshape(a , [3 1]);

  ang = norm(th);
  if abs(ang) > pi/1.01
    th = th - 2*th/ang*floor(ang + pi)/2;
    ang = norm(th);
  end

  if abs(ang) < 1e-12 %eps(0.0)
    Xi = -0.5*Spin(a);
    return
  end
  cs  = cos(ang);
  sn  = sin(ang);
  c1  = (ang*cs-sn)/ang^3;
  c2  = (ang*sn+2*cs-2)/ang^4;
  c3  = (3*sn - 2*ang - ang*cs)/ang^5;
  c4  = (1-cs)/ang^2;
  c5  = (ang-sn)/ang^3;

  I = eye(3);
  Xi = c1*(a*th') - c2*(cross(th,a)*th') + c3*(th'*a)*(th*th')  ...
     + c4*spin(a) + c5*(th'*a*I + th*a');
  return

  %-----------------------------
  n  = th/ang;
  b1 = (sin(ang/2)/(ang/2))^2;
  c1 = -(sn/ang - b1);
  c2 = 0.5*b1;
  c3 = (cs - sn/ang)/ang;
  c4 = (1-sn/ang)/ang;
  Xi = c1 * cross(n,a)*n' + c2*spin(a) + c3*(a*n' - (n'*a)*n*n') + c4*(n*a' - 2*(n'*a)*n*n' + (n'*a)*eye(3));
end
