%% From Veronique's 'FullAnalysisWF.m'
    % torsional cosntant WF

    J = 1/3*((d-2*tf)*tw^3+2*tf^3*bf)
    G = E/(2*(1+nu))
    P*L/(G*J)

    T   = P;
    Iw  = tf*(d-tf)^2*bf^3/24;
    alp = sqrt(G*J/(E*Iw));
    A   =  T/(G*J*alp)*sinh(alp*L)/cosh(alp*L);
    B   = -T/(G*J*alp);
    C   = -A;
    theta = @(x) A*cosh(alp*x) + B*sinh(alp*x) + T/(G*J)*x + C;
    theta(L)
    T*L/(G*J)+C

%%
  function [A,I,ks] = WfConstants(bf,d,tf,tw,nu)
    h1 = d - 2*tf;

    I = 1/12*tw*h1^3 + 2*(1/12*bf*tf^3 + bf*tf*(d+h1)^2/4^2);

    m = 2*bf*tf/(d*tw);
    n = bf/d;

    ksdiv = (12 + 72*m + 150*m^2 + 90*m^3) + nu*(11 + 66*m + 135*m^2 + 90*m^3) + 30*n^2*(m + m^2) + 5*nu*n^2*(8*m + 9*m^2);
    ks = 10*(1+nu)*(1+3*m)^2/ksdiv;

    A = 2*bf*tf + h1*tw;
  end

%%
  function [A,I,ks] = WfConstants2(bf1,bf2,d,tf1,tf2,tw)
    % assume tf1 = tf2 at first
    h1 = d - (tf1 + tf2);

    I = 2/12*tw*h1^3 + (1/12*bf1*tf1^3  +  bf1*tf1*(d + h1)^2/4^2) + (1/12*bf2*tf2^3  +  bf2*tf2*(d + h1)^2/4^2);
    I = 1/12*tw*h1^3 + ( bf1*tf1*(d + h1)^2/4^2) + (bf2*tf2*(d + h1)^2/4^2);

    A = bf1*tf1 + bf2*tf2 + h1*tw;

    Aweb = h1*tw;
    ks = Aweb/A;
  end

%%
  function ShapeData = RectConstant(b,d)

    [yfib,zfib,wfib] = RectPatch2Fiber([d, -b; -d, b]/2,IntTyp,nyfib,nzfib);

    switch Shear
    case "Cons"     % 3D ConsShear
      Shy = 5/6;    % Shear coefficient
      Shz = 5/6;
      cy = sqrt(Shy);
      cz = sqrt(Shz);
      as = @(yi,zi) [ 1  -yi    0     0   zi     0   ;
                      0    0   cy   -zi    0     0   ;
                      0    0    0    yi    0    cz  ];
    case "Para"
      % 2D ParaShear
      as = [1 -yi     0;
            0   0   5/4*(1-4*(yi^2/d^2))];
    end

    J = b^3*d/3;
    for i=1:500
       lam = (2*(i-1) + 1)*pi/b*d/2;
       J = J-64*b^4/(pi^5)*tanh(lam)/(2*(i-1) + 1)^5;
    end
  end

%%
function [A,I,ks] = BoxwFConstants(b,bf,d,tf,tw)
    % assume tf1 = tf2 at first
    h = d-(2*tf);
    bf1 = b+2*bf;
    bf2 = b;

    I = 1/12*2*tw*h^3 + (1/12*bf1*tf^3 + bf1*tf*(d+h)^2/4^2) ...
      +(1/12*bf2*tf^3 + bf2*tf*(d+h)^2/4^2);

    %I = 1/12*2*tw*h^3+( bf1*tf*(d+h)^2/4^2)+(bf2*tf*(d+h)^2/4^2);

    A = bf1*tf + bf2*tf + h*2*tw;

    % change for the centroid
    dy = 2*bf*tf*(d/2-tf/2);
    I = I-(dy^2)/A;

    Aweb = h*2*tw;

    ks = Aweb/A;
      %% BoxwCsShear
      cz = sqrt(Shz);
      cy = sqrt(Shy);
      function as = as_matrix(yi, zi)
        if m<(nyfl*nzfl+1) || m>(nyfl*nzfl+nyw*nzw*2) 
          % flanges
          as = [1  -yi   0    0  zi   0 ;
                0    0   0  -zi   0   0 ;
                0    0   0   yi   0  cz];
        else 
          % webs
          as = [1  -yi   0    0  zi   0;
                0    0  cy  -zi   0   0;
                0    0   0   yi   0   0];
        end
      end
      %% ParaShear
      % Shear correction factor
      m = 2*b*tf/d/(2*tw);   alpha = m;
      n = b/d;
      beta = (1+3*alpha)*(2/3)/((1+2*alpha)^2-2/3*(1+2*alpha)+1/5);

      if( m<(nyfl*nzfl+1) || m>(nyfl*nzfl+nyw*nzw*2)) % flanges
        as = [1 -yi  0     0 zi 0;
             %0   0  0     -0*zi 0 0;
              0   0     beta*(2*alpha)*(zi/b)  yi 0 1];
      else % webs
        as = [1 -yi 0 0 zi 0;
              0     0    beta*((1+2*alpha)-(2*yi/d)^2) -zi 0 0];
             %0     0    0  0*yi 0 0];
      end

end

