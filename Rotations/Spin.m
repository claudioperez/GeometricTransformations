function S = Spin(u)
%SPIN determine the spin tensor of a vector
%  S = SPIN (U)
%  the function determines the spin tensor S of a vector U with three components
 
%  =========================================================================================
%  FEDEASLab - Release 5.2, July 2021
%  Matlab Finite Elements for Design, Evaluation and Analysis of Structures
%  Professor Filip C. Filippou (filippou@berkeley.edu)
%  Department of Civil and Environmental Engineering, UC Berkeley
%  Copyright(c) 1998-2021. The Regents of the University of California. All Rights Reserved.
%  =========================================================================================
%  function by Veronique LeCorvec                                                    08-2008
%  refactored by Claudio M. Perez                                                       2023
%  -----------------------------------------------------------------------------------------

S = [   0    -u(3)  u(2)  
       u(3)    0   -u(1)  
      -u(2)   u(1)   0   ];
