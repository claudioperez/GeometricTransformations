function H = dLogC91(th, Lm)
  H = trace(Lm)*eye(3) - Lm;
  for i=1:3
    H(i,:) = H(i,:)./(2*cos(th(i)));
  end
end
