# The default point kernel stays specialized; Optimization.jl is the open path.

"""Allocation-free active-face Levenberg-Marquardt point solver."""
struct ActiveFaceLM end

"""Use any compatible Optimization.jl algorithm for each active simplex face."""
struct OptimizationInverse{A, K}
    algorithm::A
    solve_kwargs::K
    identity::Union{Nothing, String}
end
OptimizationInverse(algorithm; identity = nothing, kwargs...) =
    OptimizationInverse(algorithm, (; kwargs...), identity === nothing ? nothing : String(identity))
solver_description(s::OptimizationInverse) =
    solver_description(SolverSpec(s.algorithm; solve_kwargs = s.solve_kwargs, identity = s.identity))

_algorithm_component(c::OptimizationInverse) = c.algorithm

mutable struct OptimizationScratch{A, K}
    algorithm::A
    solve_kwargs::K
    function_evaluations::Int
    gradient_evaluations::Int
end
solver_scratch(::ActiveFaceLM) = SolverScratch()
solver_scratch(s::OptimizationInverse) = OptimizationScratch(s.algorithm, s.solve_kwargs, 0, 0)

function _polish!(
        sc::OptimizationScratch, model, rgb::SVector{3, T},
        active::NTuple{K, Int}, c::SVector{4, T}, settings
    ) where {T, K}
    K == 1 && return c, _sse(rgb, _mix(model, c)), 0, true
    x = collect(_theta_from_c(active, c))
    concentrations = u -> _c_from_active(active, _softmax_k(SVector{K, eltype(u)}(u..., zero(eltype(u)))))
    f = function (u, p)
        sc.function_evaluations += 1
        return _sse(rgb, _mix(model, concentrations(u)))
    end
    J = zeros(T, 3, 4)
    grad! = function (G, u, p)
        sc.gradient_evaluations += 1
        cc = concentrations(u)
        _jacobian4!(J, model, cc)
        dc = ForwardDiff.jacobian(concentrations, u)
        residual = _mix(model, cc) - rgb
        sc.function_evaluations += 1
        G .= 2 .* (transpose(dc) * transpose(J) * residual)
        return nothing
    end
    optf = Optimization.OptimizationFunction(f; grad = grad!)
    prob = Optimization.OptimizationProblem(optf, x)
    sol = Optimization.solve(
        prob, sc.algorithm;
        merge((maxiters = settings.max_iterations,), sc.solve_kwargs)...
    )
    return concentrations(sol.u), sol.objective, sol.stats.iterations,
        Optimization.SciMLBase.successful_retcode(sol)
end

"""Reusable per-slab solver diagnostics, accumulated without per-point storage."""
Base.@kwdef mutable struct SlabSummary
    points::Int = 0
    iterations::Int = 0
    function_evaluations::Int = 0
    gradient_evaluations::Int = 0
    unconverged::Int = 0
    max_sse::Float64 = 0.0
end
function _record!(s::SlabSummary, r::UnmixResult)
    s.points += 1
    s.iterations += r.iterations
    s.function_evaluations += r.function_evaluations
    s.gradient_evaluations += r.gradient_evaluations
    s.unconverged += !r.converged
    s.max_sse = max(s.max_sse, r.sse)
    return nothing
end

"""Generate an inverse independently of point solver, evaluator, or seed policy."""
function generate_tables_inverse(
        model, st::TableSettings; settings = UnmixSettings{Float64}(50, 1.0e-10, 1.0e-6, 4),
        point_solver = ActiveFaceLM(), fine_settings = nothing, resume = nothing, on_slab = nothing, on_summary = nothing
    )
    T = _scalar_type(model)
    coarse = st.strategy == :coarse ? generate_tables_inverse(
            model,
            TableSettings(st.coarse_n, :reference, st.evaluator, 2, st.fine_max_iterations, st.threads);
            settings, point_solver
        ) : T[]
    fine = fine_settings === nothing ?
        (st.strategy == :coarse ? UnmixSettings{T}(st.fine_max_iterations, settings.tolerance, settings.face_threshold, 1) : settings) : fine_settings
    out = Vector{T}(undef, 3 * st.n^3)
    setup = () -> (solver_scratch(point_solver), Vector{SVector{4, T}}(undef, st.n), Vector{SVector{4, T}}(undef, st.n))
    body = function (ctx, k)
        summary = SlabSummary()
        started = time()
        loaded = resume === nothing ? nothing : resume(k)
        if loaded !== nothing
            copyto!(out, 3 * k * st.n^2 + 1, loaded, 1, 3 * st.n^2)
        else
            _inverse_strategy_slab!(Val(st.strategy), out, model, st, k, ctx, fine, coarse, summary)
            on_slab === nothing || on_slab(k, copy(out[(3 * k * st.n^2 + 1):(3 * (k + 1) * st.n^2)]))
        end
        return if on_summary !== nothing
            report = settings_dict(summary)
            report["resumed"] = loaded !== nothing
            report["seconds"] = time() - started
            on_summary(k, report)
        end
    end
    _run_slabs!(setup, body, st.n, st.threads)
    return out
end

function _inverse_strategy_slab!(strategy, dest, model, st, k, ctx, settings, coarse, summary)
    sc, cur, prev = ctx
    T = eltype(dest)
    uniform = SVector{4, T}(0.25, 0.25, 0.25, 0.25)
    n = st.n
    for j in 0:(n - 1)
        for i in 0:(n - 1)
            rgb = SVector{3, T}(i / (n - 1), j / (n - 1), k / (n - 1))
            r = _strategy_solve(strategy, sc, model, rgb, settings, coarse, st.coarse_n, cur, prev, i, j, uniform)
            cur[i + 1] = r.c
            o = _table_offset(n, i, j, k)
            dest[o + 1], dest[o + 2], dest[o + 3] = r.c[1], r.c[2], r.c[3]
            _record!(summary, r)
        end
        cur, prev = prev, cur
    end
    return nothing
end
_strategy_solve(::Val{:reference}, sc, model, rgb, settings, coarse, cn, cur, prev, i, j, uniform) =
    unmix_reference(model, rgb; settings, scratch = sc)
_strategy_solve(::Val{:coarse}, sc, model, rgb, settings, coarse, cn, cur, prev, i, j, uniform) =
    unmix_bulk!(sc, model, rgb, (coarse_seed(coarse, cn, rgb...), uniform), settings)
function _strategy_solve(::Val{:bulk}, sc, model, rgb, settings, coarse, cn, cur, prev, i, j, uniform)
    if j > 0 && i > 0
        return unmix_bulk!(sc, model, rgb, (cur[i], prev[i + 1]), settings)
    elseif j > 0
        return unmix_bulk!(sc, model, rgb, (prev[i + 1],), settings)
    elseif i > 0
        return unmix_bulk!(sc, model, rgb, (cur[i],), settings)
    end
    return unmix_bulk!(sc, model, rgb, (uniform,), settings)
end
