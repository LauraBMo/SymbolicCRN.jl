
const PRINT_COMPUTING = true

macro computing(ex, msgs...)
    msg_body = isempty(msgs) ? ex : msgs[1]
    msg = string(msg_body)
    if PRINT_COMPUTING
        return quote
            print("Computing " * $msg * "...")
            local time = @elapsed out = $(esc(ex))
            print("   Finished! (time: $(time)s)\n")
            out
        end
    end
    return esc(ex)
end

exponents_matrix(p) = reduce(hcat, exponent_vectors(p))

###############################################################################
#                             Array manipulations                             #
###############################################################################

anynonzero(A) = any((!iszero).(A))

function findfirstnonzero(A)
    return findfirst((!iszero), A)
end

function findpivotsof(W)
    # return findfirstnonzero.(eachrow(W)) ## it works just for matrices, code below is more generic.
    # return mapslices(findfirstnonzero, W, dims=collect(2:ndims(W)))
    return n -> findfirstnonzero.(eachslice(W, dims = n))
end

function nonzeroslicesof(A)
    # return n->vec(mapslices(anynonzero, M, dims=deleteat!(collect(1:ndims(M)),n)))
    return n -> anynonzero.(eachslice(A, dims = n))
end

function dropzeroslices(A)
    indices = nonzeroslicesof(A).(1:ndims(A))
    return A[indices...]
end

"""
$(SIGNATURES)

Given an array of rational numbers `A`, returns a multiple `λA` whose entries are all integers, where `λ=abs(lcm(denominator.(A)))` is the minimal integer with this property.
"""
function integermultiple end

function integermultiple!(A; dims)
    for a in eachslice(A; dims = dims)
        a .= integermultiple(a)
    end
    return Int.(A)
end

integermultiple(A; dims = nothing) =
    isnothing(dims) ? Int.(abs(lcm(denominator.(A))) .* A) : integermultiple!(copy(A); dims)

integermultiple(M::AbstractMatrix) = integermultiple(M, dims = 2)

# linearcombination(M) = sum(eachcol(M))
# linearcombination(A, range) = A * rand(range, size(A, 2), 1)

## Filtering real-positive elements from an array.
isrealof(rtol) = z -> abs(imag(z)) < rtol #
filter_real(A, rtol) = real.(filter(isrealof(rtol), A))
filter_positive(A, ispos_fun = _ispositive) = filter(ispos_fun, A)

function filter_realpositive(A, rtol::Real = 1e-7)
    if isempty(A)
        print("No filtering for empty arrays.\n")
        return A
    end
    # minimag, imin = findmin(abs.(imag.(out)))
    # print("Roots found!! Min image part: $(out[imin])\n")
    out = filter(_ispositive, filter_real(A, rtol))
    if isempty(out)
        print("No one is real and positive... Returning empty array.\n")
        return out
    end
    # print("Real positive found!!\n")
    return out
end

## TODO From Polymake matrices to Juila arrays
## Type T in @convert_to macro!!
# function convert_to_matrix(M::Polymake.Matrix, ::Type{T}=Int) where {T}
#     mat = Polymake.@convert_to Matrix{T} M
#     return T.(permutedims(Matrix(mat)))
# end

# function nonzero_minors(A, k)
#     row_indices = combinations(1:nrows(A), k)
#     col_indices = combinations(1:ncols(A), k)
#     mins = Vector{typeof((first(row_indices), first(col_indices)))}(undef, 0)
#     for ri in row_indices
#         for ci in col_indices
#             if !(iszero(det(A[ri, ci])))
#                 push!(mins, (ri, ci))
#             end
#         end
#     end
#     return mins
# end
