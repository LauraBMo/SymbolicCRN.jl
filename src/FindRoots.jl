
function tcoeffs(coeffs, texps)
    maxdeg = maximum(texps)
    mindeg = minimum(texps)
    ## If mindeg>0 we divide by t^(mindeg), otherwise we multiply by it.
    ## Both correspond to shift by -mindeg!
    polyn = zeros(Int, maxdeg - mindeg + 1)
    for (c, i) in zip(coeffs, texps)
        polyn[i - mindeg + 1] += c
    end
    return polyn
end

"""

Structure of a tpolynomial T. Vector p of coefficients of T, with v[1] != 0 (nonzero independent term) and integer mindeg idicating the smallest degree of T.
mindeg can be negative, we have, with n = length(v)-1,

     T = x^mindeg * sum( v .* [1, x,..., x^n])

"""
struct tpoly{T}
    p::Vector{T}
    mindeg::Int
end

function tpoly(exps, coeffs, tdir)
    ## Compute exponent of t (to know the minimum in case it is negative).
    texps = [LA.dot(exp, tdir) for exp in exps]
    poly = tcoeffs(coeffs, texps)
    return tpoly{eltype(poly)}(poly, minimum(texps))
end
tpoly(p::MPolyElem, tdir) = tpoly(exponent_vectors(p), coeffs(p), tdir)

Base.evalpoly(z::Number, p::tpoly) = z^(p.mindeg) * Base.evalpoly(z, p.p)
Base.evalpoly(p::tpoly) = z -> Base.evalpoly(z, p)
Base.evalpoly(A::AbstractArray, p::tpoly) = Base.evalpoly(p).(A)

isrealof(rtol) = z -> abs(imag(z)) < rtol #
filter_isreal(A, rtol) = real.(filter(isrealof(rtol), A))
Nemo.ispositive(x) = x > zero(x)
function realpositiveroots(polyn, rtol::Real=1e-7)
    # poly = evalexponent(p, texp)
    roots = PolynomialRoots.roots(BigFloat.(polyn))
    return filter(ispositive, filter_isreal(roots, rtol))
end

function collect_realpositiveroots(p, tdir; rtol::Real=1e-7)
    poly = tpoly(p, tdir)
    troots = realpositiveroots(poly.p, rtol)
    # return tpointof(texp).(troots)
    return troots, poly
end

function pRoots_qPossitive(p, q; rtol=1e-7, nattemps::Int=10, bound::Int=50, p_vertex_point_map=nothing)
    if p_vertex_point_map === nothing
        return pRoots_qPossitive(PolyPolyt(;p=p), PolyPolyt(;p=q); rtol=rtol, nattemps=nattemps, bound=bound)
    end
    if isa(p_vertex_point_map, AbstractString) && isfile(p_vertex_point_map)
        return pRoots_qPossitive(PolyPolyt(;p=p), PolyPolyt(;p=q);
                                 rtol=rtol, nattemps=nattemps, bound=bound,
                                 pp_vertex_point_map=DF.readdlm(p_vertex_point_map, Int))
    end
    return pRoots_qPossitive(PolyPolyt(;p=p), PolyPolyt(;p=q);
                             rtol=rtol, nattemps=nattemps, bound=bound,
                             pp_vertex_point_map=p_vertex_point_map)
end

## Find roots of p for which q is possitive
function pRoots_qPossitive(pp::PolyPolyt, pq::PolyPolyt;
                           rtol=1e-7, nattemps::Int=10, bound::Int=50,
                           pp_vertex_point_map=vertex_point_map(pp))
    pp_outercones_negvertices = []
    for j in negvertices(pp, pp_vertex_point_map)
        print("Computing Outer j=$j ... ")
        push!(pp_outercones_negvertices, outernormalcone(pp, j))
        print("Computed\n")
        print("Computing rays j=$j ... ")
    end
    for i in posvertices(pq)
        print("Computing Outer i=$i ... ")
        icone = outernormalcone(pq, i)
        print("Computed\n")
        for (j, jcone) in enumerate(pp_outercones_negvertices)
            print("Computing rays j=$j ... ")
            rays = raysof(Polymake.polytope.intersection(icone, jcone))
            if !(isempty(rays))
                print("Cones i: $i and j: $j with nontrivial intersection\n")
                tdir = integermultiple(linearcombination(rays))
                print("Computing real positive roots\n")
                proots, tpolyp = collect_realpositiveroots(pp.p, tdir; rtol=rtol)
                qvals = evalpoly(proots, tpoly(pq.p, tdir))
                r = findfirst(ispositive, qvals)
                if r === nothing
                    j = 1;
                    while r === nothing && j < nattemps
                        j += 1
                        tdir = integermultiple(linearcombination(rays, bound))
                        print("Computing real positive roots $j\n")
                        proots, tpolyp = collect_realpositiveroots(pp.p, tdir; rtol=rtol)
                        qvals = evalpoly(proots, tpoly(pq.p, tdir))
                        r = findfirst(ispositive, qvals)
                    end
                end
                if !(r === nothing)
                    printfound(i, j, tdir, proots, qvals)
                end
            end
        end
    end
end

function printfound(i, j, tdir, proots, qvals)
    print("==============================================\n")
    print("==============================================\n")
    print("==============  Points found  ================\n")
    print("Negative vertex of p: $j\n")
    print("Positive vertex of q: $i\n\n")
    print("Exponent belonging to both cones\n")
    print("tdir $(tdir)\n")
    for i in findall(ispositive, qvals)
        print("proot: $(proots[i]), qval: $qvals[i]\n")
    end
end
