function H = dLogSO3(v, option)
%  Formerly `dExpInvSO3`
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------
%
  tol = 1/20;

  % Ensure v is a column vector
  v  = reshape(v, [3 1]);
  Sv = Spin(v);

  theta = norm(v);
  if abs(theta) > pi/1.01
    v = v - 2*v/theta*floor(theta + pi)/2;
    theta = norm(v);
  end

%%
  if theta > tol
    eta = (1-0.5*theta*cot(0.5*theta))/theta^2;
  else
    eta = 1/12 + theta^2/720 + theta^4/30240 + theta^6/1209600;
  end

  H = eye(3) - 1/2*Sv + eta*Sv*Sv;
end

