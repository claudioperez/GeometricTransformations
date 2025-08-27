function mass = Mass3dFrm(density,L,area,torsion,extra)
% Mass3dFrm(leng, dens,  area, tors)
% Mass3dFrm(leng, dens, @area, tors)
% Mass3dFrm(leng, dens, @area, tors, order)
% Mass3dFrm(leng, dens,  area, tors, 'flag')
%
% opt:
%    <int> : Consistent, Gauss integration
%    'hrz' : HRZ lumping
%

  if isa(area,'function_handle')
    if nargin < 5
      mass = mass_integral(area,dens,L,6)
    else
      mass = mass_integral(area,density,L,extra)
    end
    return
  else
    if nargin < 5
      mass = mass_sparse(L,dens)
      return
    else
      switch extra
      case 'hrz'
        mass = mass_hrz(L,dens)
      end
    end
  end
end


function mass = mass_sparse(rho,L,tors)

  m = rho*L/420.0;
  L2 = L*L;

    [ 1, 1,  m*140.0;
      7, 7,  m*140.0;
      4, 4,  m*(Jx/A)*140.0;
     10,10,  m*(Jx/A)*140.0;
      1, 7,  m*70.0;
      7, 1,  m*70.0;
      4,10,  m*(Jx/A)*70.0;
     10, 4,  m*(Jx/A)*70.0;

      2,12, -m*13.0*L;
      3,11,  m*13.0*L;
     11, 3,  m*13.0*L;
     12, 2, -m*13.0*L;

      2, 2,  m*156.0;
      3, 3,  m*156.0;
      8, 8,  m*156.0;
      9, 9,  m*156.0;

      2, 8,  m*54.0;
      3, 9,  m*54.0;
      8, 2,  m*54.0;
      9, 3,  m*54.0;

      5, 5,  m*4.0*L*L;
     11,11,  m*4.0*L*L;
      6, 6,  m*4.0*L*L;
     12,12,  m*4.0*L*L;
     12, 6, -m*3.0*L*L;
      5,11, -m*3.0*L*L;
     11, 5, -m*3.0*L*L;
      6,12, -m*3.0*L*L;

      3, 5, -m*22.0*L;
      2, 6,  m*22.0*L;
      6, 2,  m*22.0*L;
      5, 3, -m*22.0*L;

      9,11, -mass(3,5);
      5, 9, -mass(3,11);
      8,12, -mass(2,6);
      6, 8, -mass(2,12);
     11, 9, -mass(3,5);
      9, 5, -mass(3,11);
     12, 8, -mass(2,6);
      8, 6, -mass(2,12)];
end

%% ---------------------------------
function mass = mass_integral(area,dens,L,order)

  [xi,w] = Gauss(order);

  f = @(xi,L) area(L*xi) * rho * N(xi,L) * N(xi,L)';

  mass = 0.0;
  for i=1:order, 
    mass = mass + w(i)*xi(i)*f(xi(i),L); 
  end
end

%% ---------------------------------
function shape = N(xi,L)

shape = ones(12,1);

shapes = struct(
    'vert', @(xi,L) [
         (1. - 3.*xi**2 + 2.*xi**3),
         (3.*xi**2 - 2*xi**3)
    ],
    'tran', @(xi,L) [
        (@(xi,L) 1. - 3.*xi**2 + 2.*xi**3),
        (@(xi,L) 3.*xi**2 - 2*xi**3)
    ],
    'sect', @(xi,L) [
        (@(xi,L) 1 - xi),
        (@(xi,L) xi)
    ],
    'plan', @(xi,L) [
        (@(xi,L) L*(xi - 2.*xi**2 + xi**3)),
        (@(xi,L) L*(xi**3 - xi**2))
    ],
    'elev', @(xi,L) [
        (@(xi,L) L*(xi - 2.*xi**2 + xi**3)),
        (@(xi,L) L*(xi**3 - xi**2))
    ],
);
end
