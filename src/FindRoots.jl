
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

struct tpoly{T}
    p::Vector{T}
    mindeg::Int
end

function tpoly(exps, coeffs, tdir)
    ## Compute exponent of t (to know the minimum in case it is negative).
    texps = [dot(exp, tdir) for exp in exps]
    poly = tcoeffs(coeffs, texps)
    return tpoly{eltype(poly)}(poly, minimum(texps))
end
tpoly(p::MPolyElem, tdir) = tpoly(exponent_vectors(p), coeffs(p), tdir)

Base.evalpoly(z::Number, p::tpoly) = z^(p.mindeg)*Base.evalpoly(z,p.p)
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

## Find roots of p for which q is possitive
function pRoots_qPossitive(pp::PolyPolyt, pq::PolyPolyt; rtol=1e-7, nattemps::Int=10, bound::Int=50)
    for i in posvertices(pq)
        icone = outernormalcone(pq, i)
        for (j, k) in enumerate(negvertices(pp))
            jcone = outernormalcone(pp, j)
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
