Elem.Type = {'Frame'
             'Cable'
             'Truss'
             'Shell'
             'Plane'};

%% - Shear
%----------------------------------------------------------
Elem.Shear  = {false, true};
            % { Section.Av }


%% - Strain
%----------------------------------------------------------
%               _wSS     _wCS
Elem.Strain = {'Small', 'Large', 0, 1, 2};


%% - Interpolation
%----------------------------------------------------------
Elem.Interpolation.Field  = {'Basic'   % default
                             'Displ'
                             'Mixed'
                             'Force'};

Elem.Interpolation.IntTyp = {'Lobatto'
                               '...'  };

%% - Rotation
%----------------------------------------------------------
Elem.Rotation.Update = {'Incr'
                        'Iter'
                        'Init'};

Elem.Rotation.LogOpt = {'C1', 'C2', 'Exact'}


%% Transform { Rotation }
%----------------------------------------------------------
Elem.CoroData.PullBack = {'B1', 'C1', 'C2'};
Elem.CoroData.Update = {Rotation.LogOpt};
Elem.CoroData.LogOpt = {Rotation.LogOpt}

