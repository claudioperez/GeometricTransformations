function ElemName = ElementDispatch(Options)

ndm = 3;

Defaults01 = {
   'Shear',          false 
   'Strain',        'Small'  % Exact, delta
   'Field'          'Basic'  % Direct, None, Force, Displ, Mixed
};

Defaults02 = {
  {'IntType',     'Lobatto'}
};

% Apply defaults
for i=1:length(Defaults01)
  if ~isfield(Options, Defaults01{i,1})
    Options.(Defaults01{i,1}) = Defaults01{i,2};
  end
end


ElemName = Options.Field;

if Options.Shear
  ElemName = [ElemName 'Shear' num2str(ndm) 'dFrm'];
else
  ElemName = [ElemName 'Euler' num2str(ndm) 'dFrm'];
end

if strcmp(Options.Strain, 'Large')
  ElemName = [ElemName '_wCS'];
else
  ElemName = [ElemName '_wSS'];
end



