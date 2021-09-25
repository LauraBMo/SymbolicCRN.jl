
## Compute negative vertex of a polynomial where,
## a vertex of a polynomial is a term corresponding to a vertex of its Newton polytope,'
## and the sign of a vertex (a term) is the sign of its coefficient.
##
## Here, Polymake computes/gives the correspondence between vertices and exponent_vectors
## see NegativeVertex2.jl for an alternative approach, where is C++ libflint who computes/gives
## such a correspondence.

Base.@kwdef mutable struct PolyPolyt{T}
    p::MPolyElem{T}
    pointconfiguration::Union{Polymake.BigObjectAllocated,Nothing} = nothing
end

set_pointconfiguration!(pp) = pp.pointconfiguration = pointconfiguration(pp.p)

# function set_pointconfiguration!(pp)
#     pp.pointconfiguration = pointconfiguration(pp.p)
#     return pp
# end

function Newtonpolytope(pp::PolyPolyt)
    pp.pointconfiguration === nothing && set_pointconfiguration!(pp)
    return pp.pointconfiguration.CONVEX_HULL
end

function vertex_point_map(pp::PolyPolyt)
    # polytope = Newtonpolytope(pp)
    pp.pointconfiguration === nothing && set_pointconfiguration!(pp)
    hasproperty(pp.pointconfiguration, "VERTEX_POINT_MAP") || print("\n -- Computing vertices of Newton polytope -- \n\n")
    return Array(pp.pointconfiguration.VERTEX_POINT_MAP)
end

## It may be useful to "hard-push" this property to a Polymake.BigObject
## But it feels too hackie, let's try working without it:
## The functions who need VERTEX_POINT_MAP accept it as a parameter,
## so I can use the one saved from older sessions.
## See negvertices, posvertices, pRoot_qPositive,...
# function set_vertex_poss!(pp::PolyPolyt, vertex_point_map)
#     pp.pointconfiguration === nothing && set_pointconfiguration!(pp)
#     pp.pointconfiguration.VERTEX_POINT_MAP = vertex_point_map
#     return pp
# end

function save_vertex_point_map(file::String, pp::PolyPolyt)
    open(file, "w") do io
        DF.writedlm(io, vertex_point_map(pp))
    end
end

negvertices(pp::PolyPolyt, vp_map=vertex_point_map(pp)) = [i for i in (vp_map .+ 1) if coeff(pp.p, i) < zero(base_ring(pp.p))]
negvertices(p) = negvertices(PolyPolyt(;p=p))
negvertices(p, vp_map) = negvertices(PolyPolyt(;p=p), vp_map)
negvertices(p, file::String) = negvertices(p, DF.readdlm(file, Int))

posvertices(pp::PolyPolyt, vp_map=vertex_point_map(pp)) = [i for i in (vp_map .+ 1) if coeff(pp.p, i) > zero(base_ring(pp.p))]
posvertices(p) = negvertices(PolyPolyt(;p=p))
posvertices(p, vp_map) = negvertices(PolyPolyt(;p=p), vp_map)
posvertices(p, file::String) = negvertices(p, DF.readdlm(file, Int))

# hasindependentterm(p) = coeff(p, zeros(Int, nvars(parent(p)))) != base_ring(p)(0)

# function findall_verticesNewtonpolytope(poly, polytope)
#     vertices_mat = Newtonpolytope_getvertices(polytope)
#     vertices_mon = build_monomial(parent(poly)).(eachcol(vertices_mat))
#     poss = [searchsortedfirst(collect(monomials(poly)), v, rev=true) for v in vertices_mon]
#     ## If we can assume that vertices are also ordered,
#     ## something like this might be more efficient.
#     # return searchsortedfirsts(collect(monomials(p.p)), vertices; rev=true)
#     return length(poly) < last(poss) ? poss[begin:(end-1)] : poss
#     # ## See test for the check:
#     # R, (x,y) = PolynomialRing(Nemo.ZZ, ["x", "y"], ordering=:lex)
#     # @test findall_verticesNewtonpolytope(x*y^2+y-x^2*y+x*y) == [1,2,4]
# end

# findall_verticesNewtonpolytope(poly) = findall_verticesNewtonpolytope(poly,Newtonpolytope(poly))

# # function searchsortedfirsts(a, b; kargs...)
# #     out = Vector{Int}(undef, length(b))
# #     f = (x, i) -> i-1+searchsortedfirst(a[(begin+i):end], x; kargs...)
# #     out[1] = f(b[1], 0)
# #     for i in 2:length(b)
# #         out[i] = f(b[i], out[i-1])
# #     end
# #     return out
# # end

# function Newtonpolytope_getvertices(polytope)
#     out = polytope.VERTICES[:,(begin+1):end]
#     out = Polymake.@convert_to Matrix{Integer} out
#     return Int.(transpose(Matrix(out)))
# end

# build_monomial(R, exponent) = R([one(base_ring(R))], [exponent])
# # build_context = MPolyBuildCtx(R)
# # ## Important: use one(base_ring(R))!!!
# # push_term!(build_context, one(base_ring(R)), exponent)
# # return finish(build_context)

# build_monomial(R) = exponent -> build_monomial(R, exponent)
