module SymbolicCRN

using Nemo # matrix, FlintIntegerRing, nullspace
using AbstractAlgebra
import Polymake # polytope.cone,intersection
import LinearAlgebra: dot, transpose
import PolynomialRoots


export filter_isreal, realpositiveroots, collect_realpositiveroots, pRoots_qPossitive
# using LinearAlgebra: dot
# using PolynomialRoots: roots
include("FindRoots.jl")

end
