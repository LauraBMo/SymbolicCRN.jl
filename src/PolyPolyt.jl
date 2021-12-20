
## Compute negative vertex of a polynomial where,
## a vertex of a polynomial is a term corresponding to a vertex of its Newton polytope,'
## and the sign of a vertex (a term) is the sign of its coefficient.
##
## Here, Polymake computes the correspondence between vertices and exponent_vectors.
## see NegativeVertex.old.jl for an alternative approach, where is julia and C++ libflint who computes/gives
## such a correspondence.

const ver_sufix = "_vertices.cvs"
const map_sufix = "_ver_pt_map.cvs"

"""
$(TYPEDEF)

A polynomial `p`, its Newton polytope `Newtonpolyt` (as Polymake big object, with the vertices computed) and the map vertices -> points `vertex_point_map` (vector storing the position of vertex in polynomial `p`).

Methods to create a `PolyPolyt`:
```julia
PolyPolyt(p::MPolyElem) # Requires computing the vertices of the Newton polytope.
PolyPolyt(p::MPolyElem, vertices, ver_pt_map) # Loads vertices and ver-to-pt map.
PolyPolyt(p::MPolyElem, name) # Loads them from files generated by save_polypolyt(name::String, pp)
```
"""
struct PolyPolyt{T}
    p::MPolyElem{T}
    Newtonpolyt::Polymake.BigObject
    vertex_point_map::Vector{Int}
end

"""
$(SIGNATURES)

Returns the Newton polytope of a `PolyPolyt`.
"""
Newtonpolytope(pp::PolyPolyt) = pp.Newtonpolyt

"""
$(SIGNATURES)

Returns the ver-to-pt map of a `PolyPolyt`.
"""
vertex_point_map(pp::PolyPolyt) = pp.vertex_point_map

"""
$(SIGNATURES)

Given a vertex `i` of a `PolyPolyt`, returns the position of `i` in ver-to-pt map.
That is, `vertex_point_map(pp)[vertex_index(pp, i)] == i`.
Or more interestingly,
```julia
verticesof(Newtonpolytope(pp))[:, vertex_index(pp, i)] == collect(exponent_vectors(pp.p))[:,i]
```
"""
vertex_index(pp::PolyPolyt, i::Int) = searchsortedfirst(vertex_point_map(pp), i - 1)

function Base.convert(::Type{PolyPolyt}, p)
    pt_config = Polymake.polytope.PointConfiguration(POINTS = polymake_homog(permutedims(exponents_matrix(p))))
    @computing ver_pt_map = Vector(pt_config.VERTEX_POINT_MAP) "vertices of Newton polytope"
    return PolyPolyt(p,
        Polymake.polytope.Polytope(VERTICES = pt_config.CONVEX_HULL.VERTICES),
        ver_pt_map)
end

PolyPolyt(p::MPolyElem) = convert(PolyPolyt, p)

PolyPolyt(p::MPolyElem, vertices, ver_pt_map) =
    PolyPolyt(p,
        Polymake.polytope.Polytope(VERTICES = vertices),
        vec(ver_pt_map))

PolyPolyt(p::MPolyElem, name::AbstractString) =
    PolyPolyt(p,
        readdlm(name * ver_sufix, Int),
        readdlm(name * map_sufix, Int))

# TODO save and load the polynomial as well
# PolyPolyt(p::AbstractString, vertices::AbstractString, ver_pt_map::AbstractString) where T =
#     PolyPolyt(...,
#               Polymake.polytope.Polytope(VERTICES=DF.readdlm(vertices, Int)),
#               readdlm(ver_pt_map, Int))

"""
$(SIGNATURES)

Save vertices of `pp::PolyPolyt` in file `name * ver_sufix` and ver-to-pt map in `name * map_sufix`
"""
function save_polypolyt(name::AbstractString, pp::PolyPolyt)
    open(name * ver_sufix, "w") do io
        writedlm(io, Newtonpolytope(pp).VERTICES)
    end
    open(name * map_sufix, "w") do io
        writedlm(io, vertex_point_map(pp))
    end
end

negvertices(pp::PolyPolyt, isneg_fun = _isnegative) =
    filter(i -> isneg_fun(coeff(pp.p, i)), vertex_point_map(pp) .+ 1)
negvertices(p, name::AbstractString, isneg_fun = _isnegative) = negvertices(PolyPolyt(p, name), isneg_fun = isneg_fun)
negvertices(p, isneg_fun = _isnegative) = negvertices(PolyPolyt(p), isneg_fun = isneg_fun)

# posvertices(pp::PolyPolyt, ispos_fun=_ispositive) =
#     filter(i -> ispos_fun(coeff(pp.p, i)), vertex_point_map(pp) .+ 1)
# posvertices(p, name::AbstractString, isneg_fun=_isnegative) = posvertices(PolyPolyt(p, name),isneg_fun=isneg_fun)
# posvertices(p, isneg_fun=_isnegative) = posvertices(PolyPolyt(p),isneg_fun=isneg_fun)
