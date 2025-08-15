function [State,Post,SolStrat] = MultiStep(Model,ElemData,Loading,nostep,SolStrat, ...
                                           State,Post)
%% MULTISTEP : function for multi-step incremental analysis including initialization
%
%  =========================================================================================


%% Initialize state
% This will call 'init' for each element
if nargin < 6 % State not supplied
  % TODO: evalc is used to suppress output from Initialize_State
  if SolStrat.Output > 1
    State = Initialize_State(Model,ElemData); 
  else
    [~,State] = evalc('Initialize_State(Model,ElemData)');
  end
end

if nargin < 7 % Post not supplied
  [State,SolStrat] = Initialize(Model,ElemData,Loading,State,SolStrat);

  pc = 1;    % initialize plot counter
  State    = Structure('stif',Model,ElemData,State);
  Post(pc) = Structure('post',Model,ElemData,State);
  if SolStrat.Debug
    Post(pc).Debug.IterState = [];
  end
  if SolStrat.PUHist
    Post(pc).Punrm = [0];
    Post(pc).Uf = zeros(Model.nf, 1);
    % Post(pc).PUHist = struct("Punrm", 0); %State.PUHist;
  end
else
  pc = length(Post);
end

% TODO(cmp): Change to use the names already used by FEDEASLab
if ~isfield(SolStrat, 'tot_iter'), SolStrat.tot_iter = 0; end
if ~isfield(SolStrat, 'max_iter'), SolStrat.max_iter = 0; end

if SolStrat.Output == -1
  progbar = waitbar(0);
end

%% multi-step incremental analysis
tic;
for k = 1:nostep
  switch SolStrat.Output
    case -1
      waitbar(k/nostep, progbar);

    case 0

    otherwise
      % fprintf("\nLoad Step %f\n", k);
  end
  
%% Increment
  [State, SolStrat] = Increment(Model,ElemData,Loading,State,SolStrat);


%% Iterate
  if State.ConvFlag
    [State,SolStrat] = Iterate(Model,ElemData,Loading,State,SolStrat);
  else
    SolStrat.ConvFlag = false;
  end

%% Post
% if k == 1
  if ~isfield(SolStrat, 'avg_iter'),
%   SolStrat.avg_iter = 0; 
    avg_iter = SolStrat.Hist.no_iter;
  else
    n = State.Time;
    avg_iter = ((n-1)*SolStrat.avg_iter + SolStrat.Hist.no_iter)/(n);
  end

  SolStrat.avg_iter = avg_iter;
  SolStrat.tot_iter = SolStrat.tot_iter + SolStrat.Hist.no_iter;
  SolStrat.max_iter = max(SolStrat.max_iter, SolStrat.Hist.no_iter);

  if SolStrat.ConvFlag
    State = Update_State(Model,ElemData,State);
    % GoodState = State;
    % update counter and store results
    pc = pc+1;
    post = Structure('post',Model,ElemData,State);
    if SolStrat.Debug
      post.Debug = State.Debug;
    end
    % if SolStrat.PUHist
    %   post.PUHist = State.PUHist;
    % end
    Post(pc) = post;
  else
    if SolStrat.Output
      fprintf('\nNo convergence at time %.2f, step %d\n',State.Time,k);
    end
    if exist('GoodState','var')
%     State = GoodState;
    end
    break
  end
end

SolStrat.tot_time = toc;
if SolStrat.Output
  fprintf("\n");
end

