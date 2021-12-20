
## Using HomotopyContinuation.jl for bifurcation diagrams

# function starting_solution(system, ksvalues; xs=HC.variables(system), start_system=:total_degree)
#     @debug "Solving system with :total_degree"
#     result = HC.solve(system(xs, ksvalues);
#                       start_system=start_system,
#                       show_progress=false,
#                       )
#     # @debug "Solving system with :polyhedral"
#     # result = HomotopyContinuation.solve(system(xs, ksvalues), show_progress=false)
#     return realpositve_solutions(result)
# end

# function track_solution(system, solution, start, target)
#     result = HC.solve(system, solution; start_parameters=start, target_parameters=target, show_progress=false)
#     return realpositve_solutions(result)
# end

function realpositve_solutions(result;
                               positive_tol = -1e-6,
                               only_nonsingular = false,
                               kwargs...)
    solutions = HC.real_solutions(result; only_nonsingular=only_nonsingular, kwargs...)
    return filter(s -> minimum(s) > positive_tol, solutions)
end

function solve_all!(out, hc_system, K_span;
                    variables = HC.variables(hc_system),
                    start_system = :total_degree,
                    show_progress = false,
                    kwargs...)
    for j in 1:length(K_span)
        @debug "Iteration $j, point and Ks $(Ks[j])"
        result = HC.solve(hc_system(variables, [K_span[j]]);
                          start_system = start_system,
                          show_progress = show_progress,
                          kwargs...)
        # @debug "Solving system with :polyhedral"
        # result = HomotopyContinuation.solve(system(xs, ksvalues), show_progress=false)
        out[j] = realpositve_solutions(result)
        # out[j] = starting_solution(hc_sys, [Ks[j]]; start_system=start_system)
        # solution = track_solution(hc_sys, out[j], [Ks[j]], [Ks[j + 1]])
        @debug "Done."
    end
end

function solve_all(hc_system, K_span; kwargs...)
    out = Vector{Any}(undef, length(K_span))
    solve_all!(out, hc_system, K_span; kwargs...)
    return out
end


# function solve_all(hc_system, Ks;
#                    variables = HC.variables(hc_system),
#                    start_system = :total_degree,
#                    show_progress = false,
#                    kwargs...)
#     return solve_all!(Vector{Any}(undef, length(Ks)), hc_system, Ks;
#                       variables = variables,
#                       start_system = start_system,
#                       show_progress = show_progress,
#                       kwargs...)
# end

## Defining hc_system and hc_homo
# hc_eqs_smooth = F_smooth(hc_vars[new_ind], hc_params)
# hc_system = HC.System(hc_eqs_smooth; variables=hc_vars[new_ind], parameters=[hc_param]);
# hc_homo = HC.Homotopy(hc_eqs_smooth, hc_vars[new_ind], hc_param)

function solve_branch!(out, hc_homo, Xs_init, K_init;
                       K_min = 0.0, K_max = 1.0,
                      dsmin = 1e-8, dsmax = 1e-2, ds = 1e-3,
                      variables = HC.variables(hc_system),
                      show_progress = false,
                      kwargs...)
    tracker = HC.Tracker(hc_homo)
    # track(tracker::Tracker, x::AbstractVector, t₁ = 1.0, t₀ = 0.0; debug::Bool = false)
    k = K_init
    x = X_init
    push!(out, (k, x))
    dds = ds
    k_new = k + dds
    while k <= K_max
        result = HC.track(tracker, x, t₁ = k, t₀ = k_new)
        if HC.is_success(result)
            k = k_new
            x = HC.solution(result)
            push!(out, (k, x))
            dds = min(dsmax, 1.1*dds)
            k_new += dds
        else
            dds = max(dsmin, dds/10)
            k_new = k + dds
        end
    end
    return out
end

# function generatedata_tracking(hc_sys, Ks, start_system=:total_degree)
#     @debug "Preallocate output... "
#     out = Vector{Any}(undef, length(Ks))
#     @debug "Done."
#     ## Compute initial solution to track along all the values of K
#     @debug "Compute initial solutions"
#     out[1] = starting_solution(hc_sys, [first(Ks)]; start_system=start_system)[1]
#     @debug "Done."
#     for j in 1:(length(Ks) - 1)
#         @debug "Iteration $j, point $(out[j]) and Ks $(Ks[j])"
#         solution = track_solution(hc_sys, out[j], [Ks[j]], [Ks[j + 1]])
#         @debug "Done."
#         @debug "Solution computed $(solution)"
#         if length(solution) < 1
#             @debug "Recomputing Initial point"
#             solution = starting_solution(hc_sys, [Ks[j + 1]])[1]
#         end
#         out[j + 1] = solution[1]
#     end
#     return out
# end
