module SymbolicCRN

using Nemo # matrix, FlintIntegerRing, nullspace
using AbstractAlgebra
import Polymake # polytope.cone,intersection
import LinearAlgebra: dot, transpose
import PolynomialRoots


export dissect, Diagonal, partialdervativeof, Jacobian
include("NemoUtils.jl")


export filter_isreal, realpositiveroots, collect_realpositiveroots, pRoots_qPossitive
# using LinearAlgebra: dot
# using PolynomialRoots: roots
include("FindRoots.jl")


export cone_positiveorthant, cone_vectorspace, cone_positivenullspace,
    findnegativepoint, rays_outernormalcone, outernormalcone
## Nemo: matrix, FlintZZ, nullspace
## Polymake: polytope:cone,intersection
include("Cones.jl")


export negvertices, posvertices
include("PosNegVertices.jl")


export evenodd_coeffs, evenodd_polys, mix_evenodd, principaldeterminantsof,
    normalize_Hurwitz, Hurwitzmatrix, Hurwitzmatrix_sylvester, Hurwitzdeterminants,
    Hurwitz_odd_determinants, Hurwitz_subresultants, Hurwitzdeterminants_subresultant,
    signs_stabilitymatrices
include("Stability.jl")


export CharpolyCoeff, JacobianConvexparameters, JacobianDeterminantConvexparameters
incudel("ConvexParameters.jl")

# using Catalyst
# include("Catalyst.jl")

end
