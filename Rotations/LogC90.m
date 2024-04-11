function th = LogC91(R)
  th = asin(Axial(skew(R)));
end

function S = skew(M)
  S = 0.5*(M - M');
end
