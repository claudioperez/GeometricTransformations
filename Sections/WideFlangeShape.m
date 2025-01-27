function SecData = WideFlangeShape(SecData)

if ~isfield(SecData,'d'),error('section depth missing');end
if ~isfield(SecData,'tw'),error('web thickness missing');end
if ~isfield(SecData,'bf'),error('flange width missing');end
if ~isfield(SecData,'tf'),error('flange thickness missing');end

if ~isfield(SecData,'nft')
  warning('no of flange integration points in y missing, 10 IP assumed');
  SecData.nft = 10;
end
if ~isfield(SecData,'nfl')
  warning('no of flange integration points in z missing, 10 IP assumed');
  SecData.nfl = 10;
end
if ~isfield(SecData,'nwl')
  warning('no of web integration points in y missing, 10 IP assumed');
  SecData.nwl = 10;
end
if ~isfield(SecData,'nwt')
  warning('no of web integration points in z missing, 10 IP assumed');
  SecData.nwt = 10;
end

if ~isfield(SecData,'IntTyp')
  warning('integration type missing, midpoint assumed');
  SecData.IntTyp = 'Midpoint';
end

if ~isfield(SecData,'FlangeShear')
  SecData.FlangeShear = 'no';
end

SecData.SType     = 'I-Shp';      % identify shape for plotting
SecData.ShearDist = 'VLC';

[yfib,zfib,wfib,~] = Create_IPMesh4WFShape(SecData);

% TODO: Homogeneous
if ~isfield(SecData,'MatName'), error('material type missing');end
if iscell(SecData.MatName), MatName = SecData.MatName{1}; else, MatName = SecData.MatName; end
if iscell(SecData.MatData), MatData = SecData.MatData{1}; else, MatData = SecData.MatData; end

MatHndl = str2func(MatName);
SecData.MatData = MatHndl('chec', 1, MatData);

for i = 1:length(yfib)
  SecData.Fibers{i}.MatName = MatName;
  SecData.Fibers{i}.MatData = SecData.MatData;
  SecData.Fibers{i}.A = wfib(i);
  SecData.Fibers{i}.y = yfib(i);
  SecData.Fibers{i}.z = zfib(i);
end

d       = SecData.d;       % total section depth
tw      = SecData.tw;      % web thickness
bf      = SecData.bf;      % flange width
tf      = SecData.tf;      % flange thickness
nyfl    = SecData.nft;     % no of integration points in flange in y
nzfl    = SecData.nfl;     % no of integration points in flange in z
nyw     = SecData.nwl;     % no of integration points in web in y
nzw     = SecData.nwt;     % no of integration points in web in z

% switch SecData.ShearDist
%   case "ParaWarp"
%     alph = 2*bf*tf/(d*tw);
%     beta = 3/2*((1+2*alph)^2-2/3*(1+2*alph)+1/5)/(1+3*alph);
% 
%     function as = as_matrix(m, yi, zi)
%       if (abs(yi)<(d/2-tf)) %% in the web
%           as = [1 -yi   0    0    zi  0     0; % -yi*zi neglect the web part
%                 0   0   1/beta*((1+2*alph)-4*yi^2/d^2) -2*zi    0   0      0];
%       else %%in the flanges
%           as = [1 -yi  0   0    zi 0  sign(yi)*(abs(yi)-(d-tf))*zi;
%                 0   0  0   2*sign(yi)*(abs(yi)-1/2*(d-tf))   0 5/4*(1-4*(zi^2/d^2)) 0 ];
%       end
%     end
% 
%   otherwise
%     % calculate parameters for variation of shear
%   %   n = bf/d;
%     alpha = 2*bf*tf/d/(2*tw);
%     beta = (1+3*alpha)*(2/3)/((1+2*alpha)^2-2/3*(1+2*alpha)+1/5);
%     
%     function as = as_matrix(m, yi, zi)
%       if m<(nyfl*nzfl+1) || m>(nyfl*nzfl+nyw*nzw*2)  % flanges
%         psi = beta*(alpha)*(abs(zi)/bf-1/2)*sign(zi);
%         as = [1 -yi   0    0  zi  0;
%               0   0  psi  yi   0  1];
%       else % web
%         psi = beta*((1+2*alpha)-(2*yi/d)^2);
%         as = [1 -yi   0    0  zi  0;
%               0   0  psi -zi   0  0];
%       end
%   end
% end
    alpha = 2*bf*tf/d/(2*tw);
    beta = (1+3*alpha)*(2/3)/((1+2*alpha)^2-2/3*(1+2*alpha)+1/5);
    
function as = as_matrix(m, yi, zi)
  if m<(nyfl*nzfl+1) || m>(nyfl*nzfl+nyw*nzw*2)  % flanges
    psi = beta*(alpha)*(abs(zi)/bf-1/2)*sign(zi);
    as = [1 -yi   0    0  zi  0;
          0   0  psi  yi   0  1];
  else % web
    psi = beta*((1+2*alpha)-(2*yi/d)^2);
    as = [1 -yi   0    0  zi  0;
          0   0  psi -zi   0  0];
  end
end

if length(SecData.MatData.irs) == 1
  SecData.as_matrix = @(m, yi, zi) [1 -yi  zi ];
else
  SecData.as_matrix = @as_matrix;
end

end % function
