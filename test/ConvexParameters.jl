
using Nemo, Catalyst

# using Revise
using SymbolicCRN

## Complete M1
rn = @reaction_network begin
    k1,  X2 → X1
    k2,  X1 → X2
    k3,  X3 → X5
    k4,  X5 → X4
    k5,  X4 → X6
    k6,  X4 → X3
    k7,  X6 → X5
    k8,  X1+X5 → X2+X3
    k9,  X1+X6 → X2+X4
end k1 k2 k3 k4 k5 k6 k7 k8 k9
Catalyst.reorder_states!(rn, [2,1,3,5,4,6])

N = netstoichmat(rn)
E = raysof(cone_positivenullspace(N))
W = conservationlaws(N)
Y = substoichmat(rn)

lnum = size(E, 2) ## 1:6
hnum = size(N, 1) ## 7:14

# S, vars = PolynomialRing(Nemo.QQ, vcat(["l[$i]" for i in 1:lnum], ["h[$i]" for i in 1:hnum]))
L, ls = PolynomialRing(Nemo.QQ, ["ls[$i]" for i in 1:lnum])
H, hs = PolynomialRing(L, ["hs[$i]" for i in 1:hnum])

## Checking the order
# m = vars[1]
# for i in 2:14 global m *= vars[i]^(i) end
# (vars, collect(exponent_vectors(m)))
Jhl = SymbolicCRN.Jacobian_cp(H, L, N, Y, E);

# Two methods to compute it.
@test ldc_charJacobian_cp(H, L, N, Y, E, W) == ldc_charJacobian_cp(H, L, N, Y, E)

# Test one method.
@test ldc_charJacobian_cp(H, L, N, Y, E, W) ==
    (ls[1]^2*ls[2]*ls[3] + ls[1]^2*ls[2]*ls[4] + ls[1]^2*ls[2]*ls[5] + ls[1]^2*ls[3]*ls[5] +
    ls[1]^2*ls[4]*ls[5] + ls[1]^2*ls[5]^2 + ls[1]*ls[2]^2*ls[3] + ls[1]*ls[2]^2*ls[4] +
    ls[1]*ls[2]^2*ls[5] + ls[1]*ls[2]*ls[3]*ls[4] + ls[1]*ls[2]*ls[3]*ls[5] + ls[1]*ls[2]*ls[4]^2 +
    2*ls[1]*ls[2]*ls[4]*ls[5] + ls[1]*ls[2]*ls[5]^2 + ls[1]*ls[3]*ls[4]*ls[5] + ls[1]*ls[4]^2*ls[5] +
    ls[1]*ls[4]*ls[5]^2 + ls[2]^2*ls[3]*ls[4] + ls[2]^2*ls[4]^2 + 2*ls[2]^2*ls[4]*ls[5] +
    ls[2]*ls[3]*ls[4]*ls[5] + ls[2]*ls[4]^2*ls[5] + 2*ls[2]*ls[4]*ls[5]^2)*hs[1]*hs[3]*hs[4]*hs[5] +
    (ls[1]^2*ls[2]*ls[3] + ls[1]^2*ls[2]*ls[4] + ls[1]^2*ls[2]*ls[5] + ls[1]^2*ls[3]*ls[5] +
    ls[1]^2*ls[4]*ls[5] + ls[1]^2*ls[5]^2 + ls[1]*ls[2]^2*ls[3] + ls[1]*ls[2]^2*ls[4] +
    ls[1]*ls[2]^2*ls[5] + ls[1]*ls[2]*ls[3]*ls[4] + ls[1]*ls[2]*ls[3]*ls[5] + ls[1]*ls[2]*ls[4]^2 +
    2*ls[1]*ls[2]*ls[4]*ls[5] + ls[1]*ls[2]*ls[5]^2 + ls[1]*ls[3]*ls[4]*ls[5] + ls[1]*ls[4]^2*ls[5] +
    ls[1]*ls[4]*ls[5]^2 + ls[2]^2*ls[3]*ls[4] + ls[2]^2*ls[4]*ls[5] + ls[2]*ls[3]*ls[4]*ls[5] +
    ls[2]*ls[4]*ls[5]^2)*hs[1]*hs[3]*hs[4]*hs[6] + (ls[1]^2*ls[2]*ls[3] + ls[1]^2*ls[2]*ls[4] +
    ls[1]^2*ls[2]*ls[5] + ls[1]^2*ls[3]*ls[5] + ls[1]^2*ls[4]*ls[5] + ls[1]^2*ls[5]^2 +
    ls[1]*ls[2]^2*ls[3] + ls[1]*ls[2]^2*ls[4] + ls[1]*ls[2]^2*ls[5] + ls[1]*ls[2]*ls[3]*ls[4] +
    ls[1]*ls[2]*ls[3]*ls[5] + ls[1]*ls[2]*ls[4]*ls[5] + ls[1]*ls[2]*ls[5]^2 + ls[1]*ls[3]*ls[4]*ls[5] +
    ls[2]^2*ls[3]*ls[4] + ls[2]^2*ls[4]*ls[5] + ls[2]*ls[3]*ls[4]*ls[5] +
    ls[2]*ls[4]*ls[5]^2)*hs[1]*hs[3]*hs[5]*hs[6] + (ls[1]^2*ls[2]*ls[3] + ls[1]^2*ls[2]*ls[4] +
    ls[1]^2*ls[2]*ls[5] + ls[1]^2*ls[3]*ls[5] + ls[1]*ls[2]^2*ls[3] + ls[1]*ls[2]^2*ls[4] +
    ls[1]*ls[2]^2*ls[5] + ls[1]*ls[2]*ls[3]*ls[4] + ls[1]*ls[2]*ls[3]*ls[5] + ls[1]*ls[3]*ls[4]*ls[5] +
    ls[2]^2*ls[3]*ls[4] + ls[2]^2*ls[4]*ls[5] + ls[2]*ls[3]*ls[4]*ls[5] -
    ls[2]*ls[4]^2*ls[5])*hs[1]*hs[4]*hs[5]*hs[6] + (ls[1]^2*ls[2]*ls[3] + ls[1]^2*ls[2]*ls[4] +
    ls[1]^2*ls[2]*ls[5] + ls[1]^2*ls[3]*ls[5] + ls[1]^2*ls[4]*ls[5] + ls[1]^2*ls[5]^2 +
    ls[1]*ls[2]^2*ls[3] + ls[1]*ls[2]^2*ls[4] + ls[1]*ls[2]^2*ls[5] + ls[1]*ls[2]*ls[3]*ls[4] +
    ls[1]*ls[2]*ls[3]*ls[5] + ls[1]*ls[2]*ls[4]^2 + 2*ls[1]*ls[2]*ls[4]*ls[5] + ls[1]*ls[2]*ls[5]^2 +
    ls[1]*ls[3]*ls[4]*ls[5] + ls[1]*ls[4]^2*ls[5] + ls[1]*ls[4]*ls[5]^2 + ls[2]^2*ls[3]*ls[4] +
    ls[2]^2*ls[4]^2 + ls[2]^2*ls[4]*ls[5] + ls[2]*ls[3]*ls[4]*ls[5] + ls[2]*ls[4]^2*ls[5] +
    ls[2]*ls[4]*ls[5]^2)*hs[2]*hs[3]*hs[4]*hs[5] + (ls[1]^2*ls[2]*ls[3] + ls[1]^2*ls[2]*ls[4] +
    ls[1]^2*ls[2]*ls[5] + ls[1]^2*ls[3]*ls[5] + ls[1]^2*ls[4]*ls[5] + ls[1]^2*ls[5]^2 +
    ls[1]*ls[2]^2*ls[3] + ls[1]*ls[2]^2*ls[4] + ls[1]*ls[2]^2*ls[5] + ls[1]*ls[2]*ls[3]*ls[4] +
    ls[1]*ls[2]*ls[3]*ls[5] + ls[1]*ls[2]*ls[4]^2 + 2*ls[1]*ls[2]*ls[4]*ls[5] + ls[1]*ls[2]*ls[5]^2 +
    ls[1]*ls[3]*ls[4]*ls[5] + ls[1]*ls[4]^2*ls[5] + ls[1]*ls[4]*ls[5]^2 + ls[2]^2*ls[3]*ls[4] +
    ls[2]^2*ls[4]^2 + ls[2]^2*ls[4]*ls[5] + ls[2]*ls[3]*ls[4]*ls[5] + ls[2]*ls[4]^2*ls[5] +
    ls[2]*ls[4]*ls[5]^2)*hs[2]*hs[3]*hs[4]*hs[6] + (ls[1]^2*ls[2]*ls[3] + ls[1]^2*ls[2]*ls[4] +
    ls[1]^2*ls[2]*ls[5] + ls[1]^2*ls[3]*ls[5] + ls[1]^2*ls[4]*ls[5] + ls[1]^2*ls[5]^2 +
    ls[1]*ls[2]^2*ls[3] + ls[1]*ls[2]^2*ls[4] + ls[1]*ls[2]^2*ls[5] + ls[1]*ls[2]*ls[3]*ls[4] +
    ls[1]*ls[2]*ls[3]*ls[5] + ls[1]*ls[2]*ls[4]^2 + 2*ls[1]*ls[2]*ls[4]*ls[5] + ls[1]*ls[2]*ls[5]^2 +
    ls[1]*ls[3]*ls[4]*ls[5] + ls[1]*ls[4]^2*ls[5] + ls[1]*ls[4]*ls[5]^2 + ls[2]^2*ls[3]*ls[4] +
    ls[2]^2*ls[4]^2 + ls[2]^2*ls[4]*ls[5] + ls[2]*ls[3]*ls[4]*ls[5] + ls[2]*ls[4]^2*ls[5] +
    ls[2]*ls[4]*ls[5]^2)*hs[2]*hs[3]*hs[5]*hs[6] + (ls[1]^2*ls[2]*ls[3] + ls[1]^2*ls[2]*ls[4] +
    ls[1]^2*ls[2]*ls[5] + ls[1]^2*ls[3]*ls[5] + ls[1]^2*ls[4]*ls[5] + ls[1]^2*ls[5]^2 +
    ls[1]*ls[2]^2*ls[3] + ls[1]*ls[2]^2*ls[4] + ls[1]*ls[2]^2*ls[5] + ls[1]*ls[2]*ls[3]*ls[4] +
    ls[1]*ls[2]*ls[3]*ls[5] + ls[1]*ls[2]*ls[4]^2 + 2*ls[1]*ls[2]*ls[4]*ls[5] + ls[1]*ls[2]*ls[5]^2 +
    ls[1]*ls[3]*ls[4]*ls[5] + ls[1]*ls[4]^2*ls[5] + ls[1]*ls[4]*ls[5]^2 + ls[2]^2*ls[3]*ls[4] +
    ls[2]^2*ls[4]^2 + ls[2]^2*ls[4]*ls[5] + ls[2]*ls[3]*ls[4]*ls[5] + ls[2]*ls[4]^2*ls[5] +
    ls[2]*ls[4]*ls[5]^2)*hs[2]*hs[4]*hs[5]*hs[6]


# J= Jhl;
# mindeg=size(J, 1) - rank(J)
# unique_signs= x -> unique(sign.(coeffs(x))))
# P, x = PolynomialRing(base_ring(J), "x")
# p = divexact(charpoly(P, J), x^mindeg);
# H = Hurwitzdeterminants(p)(1)

# signs_stabilitymatrices(filename,
#                         Jhl,
#                         unique_signs=D -> unique(sign.(vcat(collect.(coeffs.(coeffs(D)))...))))

# filename = "M1.txt"
## M1 no intermediates
