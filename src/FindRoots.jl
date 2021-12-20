
## Given a vector tdir=[v1,...,vn] and a point x₀,
## create a univariate fewnomial tpoly(t) (it may have negative exponents)
## from a MPolyElem P(x_1,...,x_n), where
## tpoly(t) = P(x₀[1]t^v1,...,x₀[n]t^vn) for all t∈R

"""
$(SIGNATURES)

Structure of a tpolynomial P. Vector p of coefficients of P, with v[1] != 0 (nonzero independent term) and integer mindeg indicating the smallest degree of P.
mindeg can be negative, we have, with n = length(v)-1,
```julia
     poly = t^(P.mindeg) * sum( P.p .* [1, t,..., t^n])
```
"""
struct tpoly{T}
    # See get_texponents
    p::Vector{T}
    mindeg::Int
end

Base.show(io::IO, tp::tpoly) = begin
    print(io, tp.p)
    print(io, "\nMindeg: ")
    print(io, tp.mindeg)
end

# Given a multivariate polynomial p (given as a list of coeffs and exponents)
# and a vector tdir=[v_1,...,v_n], returns the tpoly
# (univariate rational function a_n*t^n+...a_0+a_(-1)*t^(-1)+... given as the list of coeffs)
# resulting of evaluates p at [t^(v_1),...,t^(v_n)].
# We "normalize" such rational function to obtain an actual polynomial with no root at t=0.
# We save the mindeg or least_deg of the original relational function in order to recover it when needed.
function tpoly(coeffs, exponents, tdir)
    ## Pre-compute exponent of t (to know the max and min).
    tpoly_exps = [dot(exp, tdir) for exp in exponents]
    (deg, least_deg) = (maximum(tpoly_exps), minimum(tpoly_exps))
    ## If least_deg>0, we divide by t^(least_deg), so we remove the trivial root t=0.
    ## Otherwise, we multiply by t^(least_deg), so we have a polynomial (with the same roots).
    ## Nicely, both cases correspond to simply shifting the exponents by -least_deg.!
    tpoly_exps = tpoly_exps .- (least_deg - 1)
    ## Pre-allocate the coefficients (most of them are zero).
    ## The new polynomial will have "shifted degree" equal to deg - least_deg.
    tpoly_coeffs = zeros(eltype(coeffs), deg - least_deg + 1)
    for (c, e) in zip(coeffs, tpoly_exps)
        tpoly_coeffs[e] += c
    end
    return tpoly(tpoly_coeffs, least_deg)
end

# TODO: use BigFloat (use string representation to pass from Nemo.RealField to BigFloat without losing precision).
function prepare_coeffs(p, x₀ = nothing, k₀ = nothing; n = 64)
    RR = RealField(n)
    if isnothing(x₀)
        if isnothing(k₀)
            # evaluates p=a_Ix^I at [t^(v_1),...,t^(v_n)], that is a_I*t^(<v,I>)
            # So, the coeffs are simply: a_I (converted to Julia numbers).
            return [convert(Float64, RR(c)) for c in coeffs(p)]
        end
        # We have "p=p_k(x)" and we evaluate p at x=[t^(v_1),...,t^(v_n)] and k=k₀.
        # So, the coeffs are: evaluate(a_I, k₀)
        k_real = RR.(k₀)
        return [convert(Float64, evaluate(c, k_real)) for c in coeffs(p)]
    end
    if isnothing(k₀)
        # evaluates p at [x_1*t^(v_1),...,x_n*t^(v_n)], that is a_I*(x₀)^I*t^(<v,I>)
        # So, the coeffs are: evaluate(a_I*x^I, x₀) where a_I*x^I is the I-th term.
        x_real = RR.(x₀)
        return [convert(Float64, evaluate(t, x_real)) for t in terms(p)]
    end
    # We have "p=p_k(x)" and we evaluate p at x=[x_1*t^(v_1),...,x_n*t^(v_n)] and k=k₀.
    # So, the coeffs are: evaluate(a_I(k), k₀)*evaluate(x^I,x₀).
    # or evaluate(evaluate(a_I(k)*x^I, x₀), k₀).
    k_real = RR.(k₀)
    x_real = RR.(x₀)
    return [convert(Float64, evaluate(evaluate(t, x_real), k_real)) for t in terms(p)]
end



# Add more as needed
# tpoly(p::MPolyElem{fmpq}, tdir, x₀=nothing) = tpoly(exponent_vectors(p), Rational.(coeffs(p)), tdir, x₀)
# tpoly(p::MPolyElem{fmpz}, tdir, x₀=nothing) = tpoly(exponent_vectors(p), Int.(coeffs(p)), tdir, x₀)
tpoly(p::MPolyElem, tdir, x₀ = nothing, k₀ = nothing; n = 64) =
    tpoly(prepare_coeffs(p, x₀, k₀; n = n), exponent_vectors(p), tdir)
tpoly(pp::PolyPolyt, tdir, x₀ = nothing, k₀ = nothing; n = 64) = tpoly(pp.p, tdir, x₀, k₀; n = n)

## Evaluate a tpoly. Via Base.evalpoly function.
Base.evalpoly(z::Number, p::tpoly) = z^(p.mindeg) * Base.evalpoly(z, p.p)
Base.evalpoly(p::tpoly) = z -> Base.evalpoly(z, p)
Base.evalpoly(A::AbstractArray, p::tpoly) = Base.evalpoly(p).(A)

# Given an MPolyElem p (with n variables) and a Vector{Int} tdir=[v_1,...,v_n],
# evaluates p at [t^(v_1),...,t^(v_n)], obtaining a univariate polynomial, and
# computes all the real and positive roots on t of such univariate polynomial.
# When a Vector{Number} x₀=[x_1,...,x_n] is given, evaluates p at [x_1*t^(v_1),...,x_n*t^(v_n)].
# When a Vector{Number} k₀=[k_1,...,k_m] is given, we assume that the coefficients of p are polynomials
# on m variables. So, we have "p=p_k(x)" and we evaluate p at x=[x_1*t^(v_1),...,x_n*t^(v_n)] and k=k₀.
function collect_realpositiveroots(p, tdir, x₀ = nothing, k₀ = nothing; n = 64, rtol::Real = 1e-7)
    poly = tpoly(p, tdir, x₀, k₀; n = n)
    # print("Roots of a polynomial: deg $(length(poly.p) - 1), mindeg: $(poly.mindeg)\n ")
    out = @computing PolynomialRoots.roots(BigFloat.(poly.p)) "roots of a tpoly"
    return filter_realpositive(out), poly
end

function from_t_to_xs(t₀::Number, tdir, x₀ = nothing, RR = RealField(64))
    out = .^(RR(t₀), tdir)
    if !(isnothing(x₀))
        out = x₀ .* out
    end
    return out
end
from_t_to_xs(tdir::AbstractVector, x₀ = nothing, RR = RealField(64)) = t₀ -> from_t_to_xs(t₀, tdir, x₀ = nothing, RR = RealField(64))

function H3root_qpositive(
    H3::MPolyElem, q::MPolyElem, H3tdir::Vector{Int},
    E, Y,
    x₀ = nothing, k₀ = nothing;
    n = 64, rtol::Real = 1e-7,
    nls = size(E, 2), nhs = size(Y, 1))

    H3roots, tH3 = collect_realpositiveroots(H3, H3tdir, x₀, k₀; n = n, rtol = rtol)
    # print("Roots found at t=$(H3roots)\n")
    RR = RealField(n)
    points = from_t_to_xs(H3tdir, x₀, RR).(H3roots)

    RRS, _ = PolynomialRing(RR, vcat(["l$i" for i in 1:nls], ["h$i" for i in 1:nhs]))
    RRq = change_base_ring(RR, q, parent = RRS)
    # CH3 = change_base_ring(CC, H3, parent=CS);
    #
    vals = [evaluate(RRq, pt) for pt in points]
    outhl = []
    out = []
    # 1-element Vector{acb}:
    #  [0.00043319620972865426 +/- 5.53e-21]
    for (v, p) in zip(vals, points)
        if ispositive(real(v))
            ks, xs = from_lh_to_kx(p, E, Y)
            Ks = [convert(Float64, real(k)) for k in ks]
            Xs = [convert(Float64, real(x)) for x in xs]
            push!(out, vcat(Ks, Xs))
            push!(outhl, p)
        end
    end
    return out, outhl
end

# function pRoots_qPossitive(p, q;
#                            rtol=1e-7, nattemps::Int=10, bound::Int=50,
#                            p_name::Union{Nothing,AbstractString}=nothing,
#                            q_name::Union{Nothing,AbstractString}=nothing)
#     return pRoots_qPossitive(PolyPolyt(p, p_name), PolyPolyt(q, q_name);
#                              rtol=rtol, nattemps=nattemps, bound=bound)
# end

# function pRoots_qPossitive(cones::Vector, q;
#                            rtol=1e-7, nattemps::Int=10, bound::Int=50,
#                            q_name::Union{Nothing,AbstractString}=nothing)
#     return pRoots_qPossitive(cones, PolyPolyt(q, q_name);
#                              rtol=rtol, nattemps=nattemps, bound=bound)
# end

# pRoots_qPossitive(pp::PolyPolyt, pq::PolyPolyt;
#                   rtol=1e-7, nattemps::Int=10, bound::Int=50) =
#                       pRoots_qPossitive(pp, outercones_negvertices(pp), pq;
#                                         rtol=rtol, nattemps=nattemps, bound=bound)

# ## Find roots of p for which q is possitive
# function pRoots_qPossitive(pp::PolyPolyt, pp_outercones::Vector, pq::PolyPolyt;
#                            rtol=1e-7, nattemps::Int=10, bound::Int=50)
#     for i in posvertices(pq)
#         print("Computing cone i=$i...")
#         icone = outernormalcone(pq, vertex_index(pq, i))
#         print("   Computed\n")
#         for (j, jcone) in enumerate(pp_outercones)
#             print("Computing rays of intersection...")
#             rays = raysof(Polymake.polytope.intersection(icone, jcone))
#             print("   Computed\n")
#             if !(isempty(rays))
#                 print("Cones i: $i and j: $j with nontrivial intersection.\n")
#                 tdir = integermultiple(linearcombination(rays))
#                 print("Computing real positive roots, 1st tiemj.\n")
#                 println(tdir)
#                 proots, tpolyp = collect_realpositiveroots(pp.p, tdir; rtol=rtol)
#                 print("Proots: ")
#                 println(proots)
#                 qvals = evalpoly(proots, tpoly(pq.p, tdir))
#                 r = findfirst(>(0), qvals)
#                 print("Qvals: ")
#                 println(qvals)
#                 # println(minimum(qvals), maximum(qvals))
#                 if isnothing(r)
#                     j = 1;
#                     while isnothing(r) && j < nattemps
#                         j += 1
#                         tdir = integermultiple(linearcombination(rays, 1:bound))
#                         print("Computing real positive roots, $j-th tiem.\n")
#                         println(tdir)
#                         proots, tpolyp = collect_realpositiveroots(pp.p, tdir; rtol=rtol)
#                         print("Proots: ")
#                         println(proots)
#                         qvals = evalpoly(proots, tpoly(pq.p, tdir))
#                         r = findfirst(>(0), qvals)
#                         print("Qvals: ")
#                         println(qvals)
#                     end
#                 end
#                 if isnothing(r)
#                     print("No positive root found for cones i=$i, j=$j\n\n")
#                 else
#                     printfound(i, j, tdir, proots, qvals)
#                 end
#             end
#         end
#     end
# end

# function printfound(i, j, tdir, proots, qvals)
#     print("==============================================\n")
#     print("==============================================\n")
#     print("==============  Points found  ================\n")
#     print("Negative vertex of p: $j\n")
#     print("Positive vertex of q: $i\n\n")
#     print("Exponent belonging to both cones\n")
#     print("tdir $(tdir)\n")
#     for i in findall(_ispositive, qvals)
#         print("proot: $(proots[i]), qval: $qvals[i]\n")
#     end
# end
