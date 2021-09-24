
function Jacobian_convexparameters(N::MatElem, Y::MatElem, E::AbstractMatrix, ls, hs)
    ## Diagonal with coordinates of E
    D = Diagonal(base_ring(N), E * ls)
    ## Square matrix with h's diagonal
    H = Diagonal(base_ring(N), hs)
    ## Return the jacobian under convex parameters
    return N * D * Y * H
end

function Jacobian_convexparameters(S,
                                   N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix;
                                   ls=gens(S)[1:size(E, 2)],
                                   hs=gens(S)[.+(size(E, 2), 1:size(N, 1))])
    return Jacobian_convexparameters(matrix(S, N),
                                     matrix(S, LA.transpose(Y)),
                                     E,
                                     ls,
                                     hs)
end

Jacobian_convexparameters(H, L,
                          N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix) =
                              Jacobian_convexparameters(H,
                                                        N, Y, E;
                                                        ls=gens(L),
                                                        hs=gens(H))

function coeff_charpoly(M::MatElem, ncoeff)
    P, x = PolynomialRing(base_ring(M), "x")
    # q = coeff(divexact(charpoly(P,Jhl),x^mindeg),0);
    return coeff(charpoly(P, M), ncoeff)
end

function det_Jacobian_convexparameters(S,
                                       N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix;
                                       ls=gens(S)[1:size(E, 2)],
                                       hs=gens(S)[.+(size(E, 2), 1:size(N, 1))],
                                       mindeg=size(N, 1) - rank(matrix(FlintZZ, N)))
    return coeff_charpoly(Jacobian_convexparameters(S, N, Y, E; ls=ls, hs=hs), mindeg)
end

function det_Jacobian_convexparameters(H, L,
                                       N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix;
                                       mindeg=size(N, 1) - rank(matrix(FlintZZ, N)))
    return coeff_charpoly(Jacobian_convexparameters(H, L, N, Y, E), mindeg)
end
