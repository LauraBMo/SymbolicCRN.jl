
"""

    cone_positiveorthant(n)

Return the cone (Polymake big object) corresponding to the nonnegative orthant of R^n.

"""
function cone_positiveorthant(n)
    rays = cat(ones(Int, n)..., dims=(1, 2))
    return Polymake.polytope.Cone(INPUT_RAYS=rays)
end

"""

    cone_vectorspace(M::AbstractMatrix{T}) where {T <: Integer}

Return the cone (Polymake big object) corresponding to vector space generated by the columns of M.

"""
function cone_vectorspace(M::AbstractMatrix{T}) where {T <: Integer}
    rays = transpose(hcat(M, -M))
    return Polymake.polytope.Cone(INPUT_RAYS=rays)
end

"""

    cone_positivenullspace(N::AbstractMatrix{T}) where {T<:Integer}

Return the cone (Polymake big object) intersection of the nonnegative orthant and the nullspace of `N`.

"""
function cone_positivenullspace(N::AbstractMatrix{T}) where {T <: Integer}
    nullspace = Nemo.nullspace_right_rational(N)
    return Polymake.polytope.intersection(cone_positiveorthant(size(N, 2)), cone_vectorspace(nullspace))
end

function findnegativepoint(p)
    Newtonp = Newtonpolytope(p)
    ver = verticesof(Newtonp, "Newton polytope")
    negver = Findallcols(isnegative, p, ver)
    points = []
    for v in negver
        rays = rays_outernormalcone(Newtonp, v)
        exponents = integermultiple(rays * ones(Int, size(rays, 2)))
        push!(points, collect_realpositiveroots(p, exponents))
    end
    return points
end

#################################
#################################
#################################
#################################

raysof(cone) = Rational.(transpose(Array(cone.RAYS)))

function rays_outernormalcone(polytope, vertex)
    cone = Polymake.polytope.normal_cone(polytope, vertex - 1, outer=1)
    return raysof(cone)
end

rays_outernormalcone(polytope) =  vertex -> rays_outernormalcone(polytope, vertex)

"""

    integermultiple(A)

Given an array of rational numbers `A`, returns a multiple `λA` with integer entries, where `λ=abs(lcm(denominator.(A)))` is the minimal integer with this property.
"""
function integermultiple(A)
    return Int.(abs(lcm(denominator.(A))) .* A)
end

linearcombination(A) = A * ones(Int, size(A, 2))
linearcombination(A, bound) = A * rand(1:bound, 1, size(A, 2))

outernormalcone(pp::PolyPolyt, vertex) =
    Polymake.polytope.normal_cone(Newtonpolytope(pp), vertex - 1, outer=1)
