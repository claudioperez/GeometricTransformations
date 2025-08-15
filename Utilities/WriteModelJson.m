function Output = WriteModelJson(Filename, Model, ElemData, Post)
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

  ElemName='Beam';

  Output.StructuralAnalysisModel.properties.sections{1} = struct(...
      'name', 1, ...
      'bounding_polygon', [
              6.0,  0.0  
              6.0,  4.0  
             -6.0,  4.0  
             -6.0,  0.0  
             -2.0,  0.0  
             -2.0, -8.0  
              2.0, -8.0  
              2.0,  0.0 ]/100);

  for i = 1:Model.nn
    Output.StructuralAnalysisModel.geometry.nodes{i}.name = num2str(i);
    Output.StructuralAnalysisModel.geometry.nodes{i}.crd  = Model.XYZ(i,:);
  end

  for i = 1:Model.ne
    % Output.nodes{i}.crd = containers.Map('crd', Model.CON{i})
    Output.StructuralAnalysisModel.geometry.elements{i}.name = num2str(i);
    Output.StructuralAnalysisModel.geometry.elements{i}.type = ElemName;
    if Model.ndm == 3
      Output.StructuralAnalysisModel.geometry.elements{i}.yvec = ElemData{i}.yornt;
      Output.StructuralAnalysisModel.geometry.elements{i}.sections = {1};
    end
    if isfield(ElemData{i}, "A")
      Output.StructuralAnalysisModel.geometry.elements{i}.area = ElemData{i}.A;
    end

    con = Model.CON{i};
    for j=1:length(con)
      Output.StructuralAnalysisModel.geometry.elements{i}.nodes{j} = num2str(con(j));
    end
  end

  Output.IterationHistory  = [];
% Output.ConvergedHistory  = [];
% Output.Displacements  = containers.Map();

  for i = 1:length(Post)

    Output.ConvergedHistory(i) = State2Mapping(Post(i), Model.DOF);
    if isfield(Post(i), "Debug") && length(Post(i).Debug.IterState) > 0
%     Output.ConvergedHistory(i) = State2Mapping(Post(i).Debug.IncrState, Model.DOF);
      Output.IterationHistory = [Output.IterationHistory State2Mapping(Post(i).Debug.IncrState, Model.DOF)];
      for j = 1:length(Post(i).Debug.IterState)
        Output.IterationHistory = [Output.IterationHistory State2Mapping(Post(i).Debug.IterState(j), Model.DOF)];
      end
    end

  end

  fid = fopen(Filename,'wt');
  fprintf(fid, jsonencode(Output,PrettyPrint=true));
  fclose(fid);
end

%%
function Resp = State2Mapping(State, DOF)
    Resp.U        = containers.Map();
    Resp.DU       = containers.Map();
    Resp.DDU      = containers.Map();
    Resp.Time     = State.Time;
%   Resp.ConvFlag = State.ConvFlag;
    for j = 1:size(DOF,1)
      Resp.U(num2str(j))   = State.U(DOF(j,:));
      if isfield(State, "DU"),  Resp.DU(num2str(j))   = State.DU(DOF(j,:));  end
      if isfield(State, "DDU"), Resp.DDU(num2str(j))  = State.DDU(DOF(j,:)); end
    end
end
