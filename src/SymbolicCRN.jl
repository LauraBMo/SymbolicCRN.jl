module SymbolicCRN

using Nemo # matrix, FlintIntegerRing, nullspace
import AbstractAlgebra as AA

using LinearAlgebra: dot # old: transpose, identity

## TODO maybe use CSV.jl
using DelimitedFiles: writedlm, readdlm
using Polymake # polytope.cone,intersection
import PolynomialRoots

using DocStringExtensions: SIGNATURES, TYPEDEF
# TODO use this instead of @computing
# using ProgressMeter

export hasproperty, findpivotsof, dropzeroslices, integermultiple!, integermultiple
include("Utils.jl")

export dissect, diagonal, partialdervativeof, Jacobian
include("NemoUtils.jl")

export PolyPolyt, Newtonpolytope, vertex_point_map, save_polypolyt, negvertices
include("PolyPolyt.jl")

export save_cones, cone_positivenullspace, raysof, verticesof, outernormalcone, outercones_negvertices
## Nemo: matrix, FlintZZ, nullspace
## Polymake: polytope:cone,intersection
include("Cones.jl")

export tpoly, collect_realpositiveroots, H3root_qpositive
# using LinearAlgebra: dot
# using PolynomialRoots: roots
include("FindRoots.jl")

export evenodd_coeffs, evenodd_polys, mix_evenodd, principaldeterminantsof,
    normalize_Hurwitz, Hurwitzmatrix, Hurwitzmatrix_sylvester, Hurwitzdeterminants,
    Hurwitz_odd_determinants, Hurwitz_subresultants, Hurwitzdeterminants_subresultant,
    signs_stabilitymatrices
include("Stability.jl")

## Use Polymake.jl to compute E = raysof(cone_positivenullspace(N))
## where N = netstoichmat(rn)
export Jacobian_cp, ldc_charJacobian_cp, from_lh_to_kx
include("ConvexParameters.jl")

export toMaple_values, toMaple_matrix
include("Maple.jl")

using Requires

function __init__()
    @require HomotopyContinuation = "f213a82b-91d6-5c5d-acf7-10f1c761b327" begin
        import .HomotopyContinuation as HC
        export solve_all
        include("HomotopyContinuation.jl")
    end
end

# using Catalyst
# include("Catalyst.jl")

end
