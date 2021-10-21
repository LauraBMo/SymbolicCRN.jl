
function tcoeffs(coeffs, texps)
    maxdeg = maximum(texps)
    mindeg = minimum(texps)
    ## If mindeg>0 we divide by t^(mindeg), otherwise we multiply by it.
    ## Both correspond to shift by -mindeg!
    polyn = zeros(eltype(coeffs), maxdeg - mindeg + 1)
    for (c, i) in zip(coeffs, texps)
        polyn[i - mindeg + 1] += c
    end
    return polyn
end

"""

Structure of a tpolynomial P. Vector p of coefficients of P, with v[1] != 0 (nonzero independent term) and integer mindeg indicating the smallest degree of P.
mindeg can be negative, we have, with n = length(v)-1,

     P = x^mindeg * sum( v .* [1, x,..., x^n])

"""
struct tpoly{T}
    p::Vector{T}
    mindeg::Int
end

Base.show(io::IO, tp::tpoly) = begin
    show(io, tp.p)
    show(io, "\nMindeg: ")
    show(io, tp.mindeg)
end

# Base.show(io::IO, ::MIME"text/plain", p::HurwitzMatrix) =
#     print(io, "HurwitzMatrix for poly (coeffs $(base_ring(p.p))) :\n", p)


function tpoly(exps, coeffs, tdir)
    ## Compute exponent of t (to know the minimum in case it is negative).
    texps = [LA.dot(exp, tdir) for exp in exps]
    poly = tcoeffs(coeffs, texps)
    return tpoly(poly, minimum(texps))
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

function pRoots_qPossitive(p, q;
                           rtol=1e-7, nattemps::Int=10, bound::Int=50,
                           p_name::Union{Nothing,AbstractString}=nothing,
                           q_name::Union{Nothing,AbstractString}=nothing)
    return pRoots_qPossitive(PolyPolyt(p, p_name), PolyPolyt(q, q_name);
                             rtol=rtol, nattemps=nattemps, bound=bound)
end

function pRoots_qPossitive(cones::Vector, q;
                           rtol=1e-7, nattemps::Int=10, bound::Int=50,
                           q_name::Union{Nothing,AbstractString}=nothing)
    return pRoots_qPossitive(cones, PolyPolyt(q, q_name);
                             rtol=rtol, nattemps=nattemps, bound=bound)
end

vertex_index(pp::PolyPolyt, i::Int) = searchsortedfirst(vertex_point_map(pp), i - 1)

function outercones_negvertices(pp::PolyPolyt)
    cones = []
    for j in negvertices(pp)
        print("Computing outer cone j=$j, vertex=$(vertex_index(pp, j))...")
        push!(cones, outernormalcone(pp, vertex_index(pp, j)))
        print("  Computed\n")
    end
    return cones
end

pRoots_qPossitive(pp::PolyPolyt, pq::PolyPolyt;
                  rtol=1e-7, nattemps::Int=10, bound::Int=50) =
                      pRoots_qPossitive(pp, outercones_negvertices(pp), pq;
                                        rtol=rtol, nattemps=nattemps, bound=bound)

## Find roots of p for which q is possitive
function pRoots_qPossitive(pp::PolyPolyt, pp_outercones::Vector, pq::PolyPolyt;
                           rtol=1e-7, nattemps::Int=10, bound::Int=50)
    for i in posvertices(pq)
        print("Computing cone i=$i...")
        icone = outernormalcone(pq, vertex_index(pq, i))
        print("   Computed\n")
        for (j, jcone) in enumerate(pp_outercones)
            print("Computing rays of intersection...")
            rays = raysof(Polymake.polytope.intersection(icone, jcone))
            print("   Computed\n")
            if !(isempty(rays))
                print("Cones i: $i and j: $j with nontrivial intersection.\n")
                tdir = integermultiple(linearcombination(rays))
                print("Computing real positive roots, 1st tiemj.\n")
                println(tdir)
                proots, tpolyp = collect_realpositiveroots(pp.p, tdir; rtol=rtol)
                print("Proots: ")
                println(proots)
                qvals = evalpoly(proots, tpoly(pq.p, tdir))
                r = findfirst(>(0), qvals)
                print("Qvals: ")
                println(qvals)
                # println(minimum(qvals), maximum(qvals))
                if isnothing(r)
                    j = 1;
                    while isnothing(r) && j < nattemps
                        j += 1
                        tdir = integermultiple(linearcombination(rays, 1:bound))
                        print("Computing real positive roots, $j-th tiem.\n")
                        println(tdir)
                        proots, tpolyp = collect_realpositiveroots(pp.p, tdir; rtol=rtol)
                        print("Proots: ")
                        println(proots)
                        qvals = evalpoly(proots, tpoly(pq.p, tdir))
                        r = findfirst(>(0), qvals)
                        print("Qvals: ")
                        println(qvals)
                    end
                end
                if isnothing(r)
                    print("No positive root found for cones i=$i, j=$j\n\n")
                else
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
