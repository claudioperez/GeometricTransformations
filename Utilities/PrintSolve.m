function PrintSolve(State, ending)
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

if nargin < 2
  ending = "";
end
fprintf("\nTime (s)   Avg. Iters  Max Iters\n");
fprintf("%f\t%f\t%f%s", State.tot_time, State.avg_iter, State.max_iter, ending);
if nargin < 2
  fprintf("\n");
end
end
