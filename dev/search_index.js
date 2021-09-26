var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SymbolicCRN","category":"page"},{"location":"#SymbolicCRN","page":"Home","title":"SymbolicCRN","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SymbolicCRN.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SymbolicCRN]","category":"page"},{"location":"#SymbolicCRN.tpoly","page":"Home","title":"SymbolicCRN.tpoly","text":"Structure of a tpolynomial T. Vector p of coefficients of T, with v[1] != 0 (nonzero independent term) and integer mindeg idicating the smallest degree of T. mindeg can be negative, we have, with n = length(v)-1,\n\n T = x^mindeg * sum( v .* [1, x,..., x^n])\n\n\n\n\n\n","category":"type"},{"location":"#AbstractAlgebra.matrix-Tuple{AbstractAlgebra.Ring, AbstractVector{T} where T}","page":"Home","title":"AbstractAlgebra.matrix","text":"matrix(R, V::AbstractVector{T}) where {T}\n\nGiven a Julia vector V of entries, construct the corresponding AbstractAlgebra.jl one-column matrix over the given ring R, assuming all the entries can be coerced into R.\n\n\n\n\n\n","category":"method"},{"location":"#SymbolicCRN.Diagonal-Union{Tuple{T}, Tuple{AbstractAlgebra.Ring, AbstractVector{T}}} where T","page":"Home","title":"SymbolicCRN.Diagonal","text":"Diagonal(R, V::AbstractVector{T}) where {T}\n\nGiven a Julia vector V of entries, construct the corresponding AbstractAlgebra.jl diagonal matrix over the given ring R, assuming all the entries can be coerced into R.\n\n\n\n\n\n","category":"method"},{"location":"#SymbolicCRN.Jacobian","page":"Home","title":"SymbolicCRN.Jacobian","text":"Jacobian(R, M)\nJacobian(R, M, vars=1:length(gens(R)))\n\nReturns the Jacobian matrix of a one-column matrix of polynomials M with respect to the generators of R indexed by vars. When vars is omitted all the generators of R are used.\n\n\n\n\n\n","category":"function"},{"location":"#SymbolicCRN.cone_positivenullspace-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T<:Integer","page":"Home","title":"SymbolicCRN.cone_positivenullspace","text":"cone_positivenullspace(N::AbstractMatrix{T}) where {T<:Integer}\n\nReturn the cone (Polymake big object) intersection of the nonnegative orthant and the nullspace of N.\n\n\n\n\n\n","category":"method"},{"location":"#SymbolicCRN.cone_positiveorthant-Tuple{Any}","page":"Home","title":"SymbolicCRN.cone_positiveorthant","text":"cone_positiveorthant(n)\n\nReturn the cone (Polymake big object) corresponding to the nonnegative orthant of R^n.\n\n\n\n\n\n","category":"method"},{"location":"#SymbolicCRN.cone_vectorspace-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T<:Integer","page":"Home","title":"SymbolicCRN.cone_vectorspace","text":"cone_vectorspace(M::AbstractMatrix{T}) where {T <: Integer}\n\nReturn the cone (Polymake big object) corresponding to vector space generated by the columns of M.\n\n\n\n\n\n","category":"method"},{"location":"#SymbolicCRN.dissect-Tuple{AbstractAlgebra.MPolyElem}","page":"Home","title":"SymbolicCRN.dissect","text":"dissect(p::MPolyElem)\n\nIterator for the coefficients and exponent vectors of the given polynomial.\n\n\n\n\n\n","category":"method"},{"location":"#SymbolicCRN.integermultiple-Tuple{Any}","page":"Home","title":"SymbolicCRN.integermultiple","text":"integermultiple(A)\n\nGiven an array of rational numbers A, returns a multiple λA with integer entries, where λ=abs(lcm(denominator.(A))) is the minimal integer with this property.\n\n\n\n\n\n","category":"method"}]
}