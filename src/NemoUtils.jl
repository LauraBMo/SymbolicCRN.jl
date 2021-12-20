
##
## A few handy methods missing in Nemo so far.

"""
$(SIGNATURES)

Iterator for the coefficients and exponent vectors of the given polynomial.

"""
function dissect(p::MPolyElem)
    return zip(coeffs(p), exponent_vectors(p))
end

get_coeffs(p, I=0:length(p)) = [coeff(p, i) for i in I]

# isneg_coeff(p, i) = isless(coeff(p, i), zero(base_ring(p)))
# ispos_coeff(p, i) = isless(zero(base_ring(p)), coeff(p, i))


_ispositive(x) = isless(zero(x), x)
_isnegative(x) = isless(x, zero(x))

# Add generic implementations (I need them for broadcasting)
AA.leading_coefficient(i::Number) = i

normalize(p) = divexact(p, content(p))

"""
$(SIGNATURES)

Given a Julia vector `V` of entries, construct the corresponding AbstractAlgebra.jl one-column matrix over the given ring `R`, assuming all the entries can be coerced into `R`.
"""
function Nemo.matrix(R, V::AbstractVector)
    return Nemo.matrix(R, reshape(V, (:, 1)))
end


"""
$(SIGNATURES)

Given a Julia vector `V` of entries, construct the corresponding AbstractAlgebra.jl diagonal matrix over the given ring `R`, assuming all the entries can be coerced into `R`.
"""
function diagonal(R, V::AbstractVector)
    n = size(V, 1)
    D = Nemo.zero_matrix(R, n, n)
    for i in 1:n D[i,i] = V[i] end
    return D
end

partialdervativeof(i::Integer) = p -> Nemo.derivative(p, i)

"""
$(SIGNATURES)

Returns the Jacobian matrix of a one-column matrix of polynomials `F` with respect to the generators of `R` indexed by `vars`. When `vars` is omitted all the generators of `R` are used.

# Examples
```julia
julia> using Nemo

julia> R, vars = PolynomialRing(ZZ, vcat(["k\$i" for i in 1:5], ["x\$i" for i in 1:4]));

julia> F = [vars[1]*vars[2]*vars[6]-vars[8]; vars[3]*vars[9]+2*vars[7]]
[k1*k2*x1-x3]
[ k3*x4+2*x2]

julia> Jacobian(R, F, 6:9)
[k1*k2  0  -1   0]
[    0  2   0  k3]

julia> Jacobian(R, F)
[k2*x1  k1*x1   0  0  0  k1*k2  0  -1   0]
[    0      0  x4  0  0      0  2   0  k3]
```
"""
function Jacobian(R::MPolyRing, F::AbstractVector, vars=1:nvars(R))
    rows = size(F, 1)
    cols = length(vars)
    J = zeros(R, rows, cols)
    for j in 1:cols
        J[:,j] = partialdervativeof(vars[j]).(F)
    end
    return J
end

"""
$(SIGNATURES)

Returns the `i`-th coefficient of the characteristic polynomial of `M`.
"""
function coeff_charpoly(M::MatElem, i)
    P, _ = PolynomialRing(base_ring(M), "x")
    return coeff(charpoly(P, M), i)
end
# # Examples
# ```jldoctest; setup = :(using SymbolicCRN, Nemo)
# julia> using Nemo
