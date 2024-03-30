function x = Axial(X)
% Return the axial vector x of the given skew-symmetric 3x3 matrix X.
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

  x = [X(3,2); X(1,3); X(2,1)];
end

