function JsonDump(Output, Filename)
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

  fid = fopen(Filename,'wt');
  fprintf(fid, jsonencode(Output,PrettyPrint=true));
  fclose(fid);
end
