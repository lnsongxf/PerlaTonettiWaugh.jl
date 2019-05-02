# Dynamic kernel
function f!(residual,du,u,p,t)
    @unpack parameters, settings, L_1, L_2, ω, Ω, Ξ₁, r̃ = p
    @unpack z_ex, T = settings
    z = z_ex[2:end-1]
    @unpack ρ, δ, σ, μ, υ, κ, d, ζ = parameters
    residual .= 0
    P = length(z) # note u[1:P]  = v(t) in the solution iterators
    v_1 = u[1]
    g = u[P+1]
    z_hat = max(u[P+2], 1.)
    E = u[P+3]

    # Get static equilibrium values and calculate the L_tilde growth rate
    @unpack S, L_tilde, z_bar, π_min, π = static_equilibrium(Ξ₁, v_1, g, z_hat, E, Ω(t), z, parameters)

    # Combine (40) with (52) to get A_t = (ρ + δ + L_tilde_log_derivative - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*I - (μ - g + (σ-1)*υ^2)*L_1 - (υ^2/2)*L_2
    # Then building the A(t) * v(t) - v'(t) residual directly for the non-algebraic equations
    residual[1:P] = (r̃(t) - (σ - 1) * (μ - g + (σ - 1) * υ^2 / 2))*u[1:P] # (40) into (C.62)
    residual[1:P] .-= (μ - g + (σ-1)*υ^2)*L_1*u[1:P] # (52) 2nd term in (53)
    residual[1:P] .-= (υ^2/2)*L_2*u[1:P] # (52) third term in (53)
    residual[1:P] .-= π # (53) final term
    residual[1:P] .-= du[1:P] # (53) subtracting the v'(t) to form the residual
    residual[P+1] = Ξ₁*v_1 + ζ - dot(ω, u[1:P])  # (54)
    residual[P+2] = z_hat^(σ-1) - κ * d^(σ-1) / π_min  # (55)
    residual[P+3] = entry_residual(v_1, Ξ₁, parameters)  # (56)
end

# Calculate the transition dynamics given a fixed Ω(t) function
function solve_dynamics(parameters, stationary_sol_T, settings, Ω, r̃)
    @unpack ρ, σ, N, θ, γ, d, κ, ζ, η, Theta, χ, υ, μ, δ = parameters
    @unpack z_ex, T, tstops = settings
    z = z_ex[2:end-1]
    P = length(z)
    @assert γ ≈ 1 # These are the only supported parameters for the transition dynamics at this point
    @assert η == 0
    T >= 50 ? nothing : @warn "T < 50 can give erratic behavior."

    # unpack the stationary solution
    v_T = stationary_sol_T.v_tilde
    g_T = stationary_sol_T.g
    z_hat_T = stationary_sol_T.z_hat
    Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1]))  # (24), with ξ = (σ-1)

    # differential objects 
    ω = ω_weights(z_ex, θ, σ-1) # Quadrature weights.
    bc = (Mixed(σ-1), Mixed(σ-1)) # boundary conditions for differential operators
    L_1 = L₁₋bc(z_ex, bc) # use backward difference as the drift is negative
    L_2 = L₂bc(z_ex, bc)

    # (backwards time) initial conditions and parameter bundle
    u0 = [v_T; g_T; z_hat_T; δ] # i.e. E(T) = δ
    du0 = zeros(P+3)
    p = (parameters = parameters, settings = settings, L_1 = L_1, L_2 = L_2, ω = ω, Ξ₁ = Ξ₁, Ω = Ω, 
            r̃ = r̃)

    dae_prob = DAEProblem(f!, du0, u0, (T, 0.0), p, differential_vars = [trues(P); false; false; false])
    return DifferentialEquations.solve(dae_prob, tstops = tstops)
end

# takes an E_hat guess and spits out a new E_hat guess
function dynamics_fixedpoint(x, params_T, stationary_sol_0, stationary_sol_T, settings)
    @show "one iteration passed."
    @unpack interp, tstops, z_ex = settings
    Nts = length(tstops)
    E_nodes_interior = x[1:Nts-2]
    r_tilde = x[Nts-1:end]

    Ω = Ω_from_E_hat(E_nodes_interior, stationary_sol_T.Ω, stationary_sol_0.Ω, params_T, settings)
    r̃ = interp(tstops, r_tilde) # callable function

    @show r_tilde

    sol = solve_dynamics(params_T, stationary_sol_T, settings, Ω, r̃)
    
    # turn the DAE solution E into an E_hat (i.e., scale it to [-1, 0])
    @unpack ρ, δ = params_T
    E = [sol(t)[end] for t in tstops] # includes boundary points 
    E_hat = (E .- δ)/(δ - E[1])
    new_E_hat_interior = E_hat[2:end-1]
    new_E_hat_interior = E_nodes_interior

    # note that `ts` and all solutions are saved backwards, i.e., solution at `T` comes first
    ts = sort(sol.t)
    gs = (x->x[end-2]).(sol.(tstops))
    z_hats = (x->x[end-1]).(sol.(tstops))
    Es = (x->x[end]).(sol.(tstops))
    Ωs = Ω.(tstops)
    Ss = S.(gs, Ref(params_T))
    L_tildes = L_tilde.(gs, z_hats, Ωs, Es, Ss, Ref(params_T))

    derivatives_forward = diff(log.(1 .- L_tildes))./(diff(tstops))
    derivatives_forward_full = [derivatives_forward; derivatives_forward[end]]
    new_r_tilde = ρ .+ δ .+ derivatives_forward_full # (C.63)

    @show new_r_tilde
    @show L_tildes

    @show norm(new_r_tilde - (ρ .+ δ*ones(length(new_r_tilde)) ))

    return [new_E_hat_interior...; new_r_tilde...]
end

function solve_transition(parameters, settings)
    @unpack T, tstops, fixedpoint_beta, fixedpoint_ftol, fixedpoint_iterations, fixedpoint_show_trace, fixedpoint_m, fixedpoint_x0, z_ex  = settings
    x0 = fixedpoint_x0(parameters, settings)
    @unpack d_0, d_T, σ = parameters
    @assert d_0 !== d_T "Will lead to a divide by 0 error."

    # calculate the steady states at 0 and T
    params_T = merge(parameters, (d = d_T,))
    params_0 = merge(parameters, (d = d_0,))
    stationary_T = stationary_numerical(params_T, settings)
    stationary_0 = stationary_numerical(params_0, settings)

    @assert tstops[1] ≈ 0.0
    @assert tstops[end] ≈ T

    fp_sol = fixedpoint(x -> dynamics_fixedpoint(x, params_T, stationary_0, stationary_T, settings),
                        x0,
                        inplace = false,
                        iterations = fixedpoint_iterations,
                        ftol = fixedpoint_ftol,
                        m = fixedpoint_m, 
                        beta = fixedpoint_beta,
                        show_trace = fixedpoint_show_trace)
    
    converged(fp_sol) || @warn "fixed point didn't converge."
    
    return prepare_results(fp_sol, stationary_T, stationary_0, params_T, settings)
end

function Ω_from_E_hat(E_nodes_interior, Ω_T, Ω_0, parameters, settings)
    @unpack δ = parameters
    @unpack T, tstops, interp = settings
    E_nodes = [-1.0; E_nodes_interior; 0.0] # See footnote 19.  Without pinning a corner there is an extra degree of freedom in the scale
    ts = range(0.0, T, length=length(E_nodes))

    E_hat_interpolation = interp(ts, E_nodes) # might worth trying cubic spline
    E_hat(t) = E_hat_interpolation(t)
    E_hat_integral = quadgk(E_hat, 0, T)[1] # (42)
    Q = log(Ω_T/Ω_0) / E_hat_integral # (42) when Ω_T = Ω_0, then Q = 0 so that E(t) is constant with δ as expected
    E(t) = Q*E_hat(t) + δ # (43)

    Ω_derivative(Ω,p,t) = Q*E_hat(t)*Ω # (44)
    Ω_solution = DifferentialEquations.solve(ODEProblem(Ω_derivative, Ω_0, (0.0, T)), reltol = 1e-15, tstops = tstops) # if this fails, error will be thrown

    function Ω(t)
        if t > T 
            return Ω_T
        elseif t > 0 && t <= T 
            return Ω_solution(t)
        else # useful for fixed point approach 
            return Ω_0
        end 
    end 

    return Ω
end

function prepare_results(fp_sol, stationary_T, stationary_0, params_T, settings)
        @unpack σ = params_T 
        @unpack z_ex, tstops, interp = settings 

        # recreate last DAE solution 
        Nts = length(tstops)
        Ω = Ω_from_E_hat(fp_sol.zero[1:Nts-2], stationary_T.Ω, stationary_0.Ω, params_T, settings)
        r̃ = interp(tstops, fp_sol.zero[Nts-1:end])
        sol = solve_dynamics(params_T, stationary_T, settings, Ω, r̃)
    
        # generate full data 
        z = z_ex[2:end-1]
        Ξ₁ = 1/(1 - (σ-1)*(z[1] - z_ex[1]))
        ts = tstops
        gs = [sol(t)[end - 2] for t in ts]
        z_hats = [sol(t)[end-1] for t in ts]
        Es = [sol(t)[end] for t in ts]
        v_1s = [sol(t)[1] for t in ts]
        Ωs = Ω.(ts)
        Ss = S.(gs, Ref(params_T))
        L_tildes = L_tilde.(gs, z_hats, Ωs, Es, Ss, Ref(params_T))
        L_tilde_xs = L_tilde_x.(z_hats, Ωs, Ref(params_T))
        L_tilde_as = L_tilde_a.(Ωs, Ss, Ref(params_T))
        L_tilde_Es = L_tilde_E.(Ωs, Es, Ref(params_T))
        z_bars = z_bar.(z_hats, Ωs, Ref(params_T))
        ws = w.(z_bars, Ref(params_T))
        xs = x.(ws, Ref(params_T))
        λ_iis = λ_ii.(z_hats, Ref(params_T))
        cs = c.(L_tildes, Ωs, z_bars)
        π_mins = π_min.(L_tildes, z_bars, Ref(params_T))
        π_rats = π_rat.(z_hats, π_mins, Ref(params_T))
        entry_residuals = entry_residual.(v_1s, Ref(Ξ₁), Ref(params_T))
    
        return (results = 
                    (
                        t = ts,
                        r̃ = r̃.(ts),
                        g = gs, 
                        z_hat = z_hats, 
                        E = Es, 
                        v_1 = v_1s, 
                        Ω = Ωs, 
                        S = Ss, 
                        L_tilde = L_tildes,
                        L_tilde_x = L_tilde_xs, 
                        L_tilde_a = L_tilde_as, 
                        L_tilde_E = L_tilde_Es, 
                        λ_ii = λ_iis,
                        z_bar = z_bars, 
                        w = ws, 
                        x = xs, 
                        c = cs, 
                        π_min = π_mins, 
                        π_rat = π_rats, 
                        entry_residual = entry_residuals
                    ),
                sol = sol)
end 
