function vect = LogSO3(R)
%LOGSO3 Inverse of the exponential map on SO(3).
%  vect = LogSO3(R)
%  vect = LogSO3(R,option)
%  Returns the axial parameters associated with the rotation `R`. The result
%  should satisfy the following equality for any 3-vector, `v`:
%
%        LogSO3(expm(spin(v))) == v
%
%  where `expm` is the Matlab built-in matrix exponential, and `spin` is a function
%  which produces the skew-symmetric 3x3 matrix associated with vector `v`.
%
%  Parameters
%    R       (3x3)   Rotation matrix.
%    option  string  Algorithm option; default is 'Quat'. See reference below.
%                    All other options are for internal use.
%
%  Remarks
%
%  - Does not check if input is really a rotation. If you arent sure, use
%    Matlab's `logm`. This is slower, but more general. Note however that`logm`
%    may not be as accurate in corner cases and sometimes returns complex-valued
%    results.
%  - The angle corresponding to the returned vector is always in the interval [0,pi].
%
%
%  References
%  1. Nurlanov Z (2021) Exploring SO(3) logarithmic map: degeneracies and 
%     derivatives.
%
%  =========================================================================================
%  function by Claudio Perez                                                            2023
%  -----------------------------------------------------------------------------------------

vect = Quat2Axis(Rmat2Quat(R));

