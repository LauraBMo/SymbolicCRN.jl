

## (De)compose a polynomial into even and odd degrees
function evenodd_coeffs(p)
    l = length(p); ## degree + 1
    ## For the zero polynomial, 'odd' would be a 0-element array
    if l < 1
        even = [0]
        odd  = [0]
    else
        even = get_coeffs(p, 0:2:l)
        odd  = get_coeffs(p, 1:2:l)
        # even = reverse([coeff(p, i) for i in 0:2:l])
        # odd  = reverse([coeff(p, i) for i in 1:2:l])
    end
    return even, odd
end

evenodd_polys(p) = parent(p).(evenodd_coeffs(p))

function mix_evenodd(even, odd)
    n = length(even) + length(odd)
    out = Vector{eltype(even)}(undef, n)
    out[1:2:end] .= even
    out[2:2:end] .= odd
    return out
end

mix_evenodd(pe::PolyElem, po::PolyElem) = mix_evenodd(get_coeffs(pe), get_coeffs(po))

principaldeterminantsof(M) = i -> det(M[1:i,1:i])

# Structure for the Hurwitz matrix of a polynomial.
# Inherits matrix structure from MatElem general type.
struct HurwitzMatrix{T} <: MatElem{T}
    p::PolyElem{T}
    # HurwitzMatrix(p::PolyElem{T}) where T = new{T}(parent(p)(normalize_Hurwitz(p)))
end

AA.nrows(p::HurwitzMatrix) = degree(p.p)
AA.ncols(p::HurwitzMatrix) = AA.nrows(p)
AA.length(p::HurwitzMatrix) = AA.nrows(p)^2
AA.isempty(p::HurwitzMatrix) = false
AA.base_ring(p::HurwitzMatrix) = base_ring(p.p)

Base.show(io::IO, p::HurwitzMatrix) = print(io, p.p)
Base.show(io::IO, ::MIME"text/plain", p::HurwitzMatrix) =
    print(io, "HurwitzMatrix for poly (coeffs $(base_ring(p.p))) :\n", p)

Base.convert(::Type{AbstractMatrix}, p::HurwitzMatrix) = Matrix(p)

Base.IndexStyle(::Type{<:HurwitzMatrix}) = IndexCartesian()
function Base.getindex(p::HurwitzMatrix, i::Int, j::Int)
    n = degree(p.p)
    # For i > deg(p), we already have `coeff(p,i)=0`
    return -1 < (n - 2 * i + j) ? coeff(p.p, n - 2 * i + j) : zero(base_ring(p.p))
end


# Functions computing Hurwitz something.
function normalize_Hurwitz(p)
    even, odd = normalize.(evenodd_polys(p))
    return parent(p)(mix_evenodd(even, odd))
end

Hurwitzmatrix(p) = HurwitzMatrix(p)[1:degree(p), 1:degree(p)]
Hurwitzmatrix_sylvester(p) = sylvester_matrix(evenodd_polys(p)...)

Hurwitzdeterminants(p) = principaldeterminantsof(HurwitzMatrix(p))

Hurwitz_odd_indices(p) = (degree(p) - 1):-2:2

Hurwitz_odd_determinants(p) = reverse(Hurwitzdeterminants(p).(Hurwitz_odd_indices(p)))

Hurwitz_subresultants(p) = Hurwitz_subresultants_signs(p) .* subresultants_ducos(evenodd_polys(p)...)
Hurwitzdeterminants_subresultant(p) = AA.lead.(Hurwitz_subresultants(p))

function Hurwitz_subresultants_signs(p)
    n = length(Hurwitz_odd_indices(p))
    out = ones(Int, n)
    goes_neg = iseven(degree(p)) ? [1,2] : [2,3]
    for i in 1:n
        if i % 4 in goes_neg
            out[i] = -1
        end
    end
    return out
end

## Numerical results form test_subresultant_ducos
## Deg | length | signs
##  22 |     10 | [-1, -1, 1, 1, -1, -1, 1, 1, -1, -1]
##  20 |      9 | [-1, -1, 1, 1, -1, -1, 1, 1, -1]
##  18 |      8 | [-1, -1, 1, 1, -1, -1, 1, 1]
##
##  21 |     10 | [1, -1, -1, 1, 1, -1, -1, 1, 1, -1]
##  19 |      9 | [1, -1, -1, 1, 1, -1, -1, 1, 1]
##  17 |      8 | [1, -1, -1, 1, 1, -1, -1, 1]

#############################################
#############################################
################################### Stability

## See https://sci-hub.se/https://doi.org/10.1017/S0305004100049872
## Barnett 1971, Proceedings of the Cambridge Philosophical Society, 70(02), 269.
##
##

# Check the signs of the Hurwitz determinants for the Ruth-Hurwitz criteria
# getsigns(D) = unique(sign.(reshape(hcat, coeffs.(coeffs(D)))))
function signs_stabilitymatrices(io,
                                 J;
                                 mindeg::Integer=size(J, 1) - rank(J),
                                 unique_signs::Function=x -> unique(sign.(coeffs(x))))
    P, x = PolynomialRing(base_ring(J), "x")
    p = divexact(charpoly(P, J), x^mindeg)
    L = []
    # println(p)
    range = Hurwitz_odd_indices(p)
    # println(collect(range))
    write(io, "==============================================\n")
    write(io, "==============================================\n")
    write(io, "$(length(range)) expressions will be studied\n")
    for i in reverse(range)
        # This poly belongs to base_ring(J), which is defined in the scope calling this function
        H = Hurwitzdeterminants(p)(i)
        write(io, "\n-- $i th determinant --\n")
        signs = unique_signs(H)
        write(io, "$(signs)\n")
        if -1 in signs
            push!(L, H)
        end
    end
    write(io, "==============================================\n")
    write(io, "==============================================\n")
    return L
end

function signs_stabilitymatrices(J;
                                 mindeg::Integer=size(J, 1) - rank(J),
                                 unique_signs::Function=x -> unique(sign.(coeffs(x))))
    signs_stabilitymatrices(IOBuffer(), J, mindeg=mindeg, unique_signs=unique_signs)
end

function signs_stabilitymatrices(file::String,
                                 J;
                                 netname=nothing,
                                 mode="a", # append
                                 mindeg::Integer=size(J, 1) - rank(J),
                                 unique_signs::Function=x -> unique(sign.(coeffs(x))))
    open(file, mode) do io
        if !(netname === nothing)
            write(io, "==============================================\n")
            write(io, "==============================================\n")
            write(io, "$netname\n")
            write(io, "==============================================\n")
            write(io, "==============================================\n")
        end
        signs_stabilitymatrices(io, J, mindeg=mindeg, unique_signs=unique_signs)
    end
end
