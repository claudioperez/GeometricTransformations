function ElemState = CreateState_3dFrm(v)
  ElemState.u   = [v(1) zeros(1,5); zeros(1,6)];
  ElemState.Du  = [v(1) zeros(1,5); zeros(1,6)];
  ElemState.DDu = [v(1) zeros(1,5); zeros(1,6)];
end

