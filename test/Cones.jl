
using Nemo, LinearAlgebra, DelimitedFiles

# using Revise
using SymbolicCRN

M = [63 90 36 -45 -18
    -7 -10 -4 5 2
    28 40 16 -20 -8
    -56 -80 -32 40 16
    7 10 4 -5 -2]

K = [2 5 0 0 0 0
    0 0 1 0 0 1
    0 0 0 5 1 0
    0 7 2 4 0 0
    7 0 0 0 2 5]

R, vrs = PolynomialRing(Nemo.ZZ, vcat(["k$i" for i in 1:5], ["x$i" for i in 1:4]));

@test all(integermultiple(raysof(SymbolicCRN.cone_positiveorthant(5))) .== LinearAlgebra.I(5))
@test integermultiple(raysof(cone_positivenullspace([0 0; 0 0]))) == [1 0; 0 1]
@test raysof(cone_positivenullspace([1 -1; 2 -2])) == transpose([1 1])
Knew = raysof(cone_positivenullspace(M))
@test all(M * Knew .== 0)
@test all(map(>=(0), Knew))
@test K == integermultiple(Knew)
p = prod([vrs[i] - i for i in 1:7])
@test integermultiple(raysof(outernormalcone(p, 5))) == [
    0 1 0 0 0 0 0
    1 0 0 0 0 0 0
    0 0 1 0 0 0 0
    0 0 0 1 0 0 0
    0 0 0 0 0 0 -1
    0 0 0 0 1 0 0
    0 0 0 0 0 1 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
]

# TODO;
# verticesof(polytope) = permutedims(Rational.(polymake_dehomog(polytope.VERTICES)))
