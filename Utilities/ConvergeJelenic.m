function norm = TestJelenic(Model,State)
  DuDu = 0;
  uu = 0;
  for i=1:Model.nn
    DuDu = DuDu + dot(State.DDU(Model.DOF(i,1:3)), State.DDU(Model.DOF(i,1:3)));
    uu   = uu   + dot(State.U(Model.DOF(i,1:3)), State.U(Model.DOF(i,1:3)));
  end
  du = sqrt(DuDu)/sqrt(uu); % /uu);
  norm = du;

% Df = 0.0;
% for i=1:Model.nn
% end
% df = sqrt(Df);

% norm = max(df, du);
end
