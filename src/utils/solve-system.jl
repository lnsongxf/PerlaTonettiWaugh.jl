# returns the m-vector zero of a multivariate function h: R^m -> R^n, using nlopt
# constraints_fg! takes (h, x, jacobian_t) and assigns jacobian_t and h[:] given current x.
function solve_system(h, x0; lb = fill(0.0, length(x0)), ub = fill(10e8, length(x0)), 
                    constraints_fg! = nothing, constraints_tol = fill(1e-8, length(x0)),
                    autodiff = (constraints_fg! == nothing) ? :forward : :finite,
                    algorithm = (constraints_fg! == nothing) ? :LD_LBFGS : :LD_SLSQP,
                    iterations = 1000, increasing_x_penalty_coefficient = 0.)
    if (iterations < 1) # if no optimization is required, return the initial solution
        return x0
    end

    function f(x)
        resids = h(x)
        return sum(resids .* resids) +
            ((increasing_x_penalty_coefficient > 0.) ? # add a penalty function for constraints on increasing E
            (increasing_x_penalty_coefficient * sum((max.(0.0, diff(x))).^2)) : 0.) 
    end

    fg!(x::Vector, grad::Vector) = f(x) # fg! for derivative free methods

    # define the optimization problem
    opt = (constraints_fg! == nothing && autodiff != :forward) ? Opt(:LN_NEWUOA_BOUND, length(x0)) : Opt(algorithm, length(x0))

    # specifies that optimization problem is on minimization
    min_objective!(opt, (constraints_fg! == nothing && autodiff != :forward) ? fg! : NLoptAdapter(f, x0, autodiff)) 

    if (lb != nothing)    
        lower_bounds!(opt, lb) # find `x` above lb
    end
    if (ub != nothing)    
        upper_bounds!(opt, ub) # find `x` below ub
    end

    if (constraints_fg! != nothing) # add constraints_fg! if needed
        inequality_constraint!(opt, constraints_fg!, constraints_tol)
    end

    xtol_rel!(opt, 1e-8)
    xtol_abs!(opt, 1e-8)
    maxeval!(opt, iterations)

    # solve the optimization problem
    (minf,minx,ret) = NLopt.optimize(opt, x0)
    return minx
end

struct NLoptAdapter{T} <: Function where T <: AbstractObjective
    nlsolver_base::T
end

# implement fg!; note that the order is reversed
(adapter::NLoptAdapter)(x, df) = adapter.nlsolver_base.fdf(df, x)
(adapter::NLoptAdapter)(result, x, jacobian_transpose) = adapter.nlsolver_base.fdf(result, jacobian_transpose', x)

# constructors
NLoptAdapter(f, x, autodiff = :forward) = NLoptAdapter(OnceDifferentiable(f, x, autodiff = autodiff))
NLoptAdapter(f!, x::Vector, F::Vector, autodiff = :forward) = NLoptAdapter(OnceDifferentiable(f!, x, F, autodiff = autodiff))