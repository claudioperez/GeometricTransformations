function GDT = TanGLn(a, b, c, N, Hat, CkT, zk, opTransp, opInvers)
%
%  [1] Todesco, J. and Brüls, O. (2023) ‘Highly accurate differentiation of the 
%      exponential map and its tangent operator’, Mechanism and Machine Theory, 
%      190, p. 105451. 
%      Available at: https://doi.org/10.1016/j.mechmachtheory.2023.105451.
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------
%
bc  = [b, c];
nm  = size(bc,2);
GDT = zeros(zk, zk);
GDa = zeros(zk, zk, 2^nm);
GDa(:, :, 1) = eye(zk, zk);
fact = 1;
N0   = 1;
if nm == 0
    N0 = 0;
end
mi = 0;
if size(c,2) > 1
    mi = size(c,2) - 1;
end
for i = N0:(N + nm - 1)
    for e = nm:1
        for v = nchoosek(1:nm, e)'
        % for f = nchoosek(nm,e):1
            % v = subSet(nm, e, f);
            % v = V(f, )
            n = 0;
            m = 0;
            for g = 1:e
                if v(g) <= size(b,2) + mi % do
                    n = n + 1;
                else
                    m = m + 1;
                end
            end % for g
            idx = 1 + sumPower2(nm - v);
            GDa(:, :, idx) = Hat(a) * GDa(:, :, idx);
            for j = [1:n  n+2:n+m]
                idx2 = 1 + sumPower2(nm - v([1:j - 1, j + 1:e]));
                GDa(:, :, idx) = GDa(:, :, idx) ...
                               + Hat( bc(:, v(j) )) * GDa(:, :, idx2);
            end % for j
            if m > 0
                idx2 = 1 + sumPower2(nm - v([1:n, n + 2:e]));
                if opTransp == 1 % do
                    GDa(:, :, idx) = GDa(:, :, idx) ...
                                   + CkT(GDa(:, :, idx2) * bc(:, v(n + 1)));
                else
                    GDa(:, :, idx) = GDa(:, :, idx) ...
                                   - Hat(GDa(:, :, idx2) * bc(:, v(n + 1)));
                end
            end
        end % for
    end % for

    GDa(:, :, 1) = Hat(a) * GDa(:, :, 1);
    if opInvers == true
        GDT = GDT + (-1)^i * GDa(:, :, end) * bernoulli(i)/fact;
        fact = fact * (i + 1);
    else
        fact = fact * (i + 1);
        GDT = GDT + (-1)^i * GDa(:, :, end)/fact;
    end % if
end % for
end % function

%%
function s = sumPower2(z)
  s = 0;
  for e = z
    s = s + 2^e;
  end
end
