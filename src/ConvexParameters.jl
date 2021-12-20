

"""
$(SIGNATURES)

Returns Nemo's `MatElem` paramatrising the Jacobian of the dynamical system associated to a biochemical reaction network. We use "convex parameters", see
```julia
rn = @reaction_network begin
    k1,  X2 → X1
    k2,  X1 → X2
end k1 k2

N = netstoichmat(rn)
Y = substoichmat(rn)

E = raysof(cone_positivenullspace(N))
W = conservationlaws(N)
```
"""
function Jacobian_cp(N::MatElem, Y::MatElem, E::AbstractMatrix, ls, hs)
    ## Diagonal with coordinates of E
    D = diagonal(base_ring(N), E * ls)
    ## Square matrix with h's diagonal
    H = diagonal(base_ring(N), hs)
    ## Return the jacobian under convex parameters
    return N * D * Y * H
end

function Jacobian_cp(S,
                     N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix;
                     ls=gens(S)[1:size(E, 2)],
                     hs=gens(S)[.+(size(E, 2), 1:size(N, 1))])
    return Jacobian_cp(matrix(S, N),
                       matrix(S, permutedims(Y)),
                       E,
                       ls,
hs)
end

Jacobian_cp(H, L,
            N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix) =
                Jacobian_cp(H,
                            N, Y, E;
                            ls=gens(L), hs=gens(H))


## Computing the coeffitient of the least degree term (ldc) of the characteristic polynomial of the Jacobian, parametrised using convex parameters.
##
## Computing the whole char polynomial
function ldc_charJacobian_cp(S,
                             N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix;
                             ls=gens(S)[1:size(E, 2)],
                             hs=gens(S)[.+(size(E, 2), 1:size(N, 1))],
                             mindeg=size(N, 1) - rank(matrix(FlintZZ, N)))
    return coeff_charpoly(Jacobian_cp(S, N, Y, E; ls=ls, hs=hs), mindeg)
end

function ldc_charJacobian_cp(H, L,
                             N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix;
                             mindeg=size(N, 1) - rank(matrix(FlintZZ, N)))
    return ldc_charJacobian_cp(H,
                               N, Y, E;
                               ls=gens(L), hs=gens(H),
                               mindeg=mindeg)
end

## Computing the coeffitient of the least degree term (ldc) of the characteristic polynomial of the Jacobian, parametrised using convex parameters.
##
## If we know W, we can use it, simply add it as a last argument!
function ldc_charJacobian_cp(S,
                             N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix, W::AbstractMatrix;
                             ls=gens(S)[1:size(E, 2)],
                             hs=gens(S)[.+(size(E, 2), 1:size(N, 1))],
                             mindeg=size(N, 1) - rank(matrix(FlintZZ, N)))
    J = Jacobian_cp(S, N, Y, E; ls=ls, hs=hs)
    for (i, j) in enumerate(findpivotsof(W)(1))
        # Change row j-th of J by i-th of W
        # Nemo do not allow broadcasting.
        for k in 1:size(J, 2)
            J[j,k] = W[i,k]
        end
    end
    return det(J)
end

function ldc_charJacobian_cp(H, L,
                             N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix, W::AbstractMatrix;
                             mindeg=size(N, 1) - rank(matrix(FlintZZ, N)))
    return ldc_charJacobian_cp(H,
                               N, Y, E, W;
                               ls=gens(L), hs=gens(H),
                               mindeg=mindeg)
end

function from_lh_to_kx(point, E, Y, nls=size(E, 2))
    xs = inv.(point[nls + 1:end])
    e = E * point[begin:nls]
    x_to_y = [prod(xs.^r) for r in eachcol(Y)]
    ks = e .* inv.(x_to_y)
    return ks, xs
end
