using Nemo, LinearAlgebra, DelimitedFiles

using Revise
using SymbolicCRN

M = [63   90   36  -45  -18
     -7  -10   -4    5    2
     28   40   16  -20   -8
     -56  -80  -32   40   16
     7   10    4   -5   -2]

K = [2  5  0  0  0  0
     0  0  1  0  0  1
     0  0  0  5  1  0
     0  7  2  4  0  0
     7  0  0  0  2  5]

R, vrs = PolynomialRing(Nemo.ZZ, vcat(["k$i" for i in 1:5], ["x$i" for i in 1:4]));
F = [vrs[1] * vrs[2] * vrs[6] - vrs[8]; vrs[3] * vrs[9] + 2 * vrs[7]];

@test Jacobian(R, F, 6:9) == [vrs[1] * vrs[2]  0  -1   0
                              0  2   0  vrs[3]]
@test Jacobian(R, F) == [vrs[2] * vrs[6]  vrs[1] * vrs[6]   0  0  0  vrs[1] * vrs[2]  0  -1   0
                         0      0  vrs[9]  0  0      0  2   0  vrs[3]]
