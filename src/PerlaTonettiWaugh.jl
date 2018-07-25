module PerlaTonettiWaugh

# package code goes here
using DifferentialEquations, NamedTuples, BandedMatrices, Sundials, Distributions
import Parameters: @with_kw, @unpack

include("diffusionoperators.jl")
include("simple/ODEalgorithm.jl")
include("simple/simple_dynamic_ODE.jl")
include("simple/DAEalgorithm.jl")
include("simple/algebraicstationary.jl")
include("simple/numericalstationary.jl")
include("simple/calculate_residuals.jl")
include("full/algebraicstationary.jl")
include("quadrature.jl")

export f!, diffusionoperators, createsimpleODEproblem, createsimpleDAEproblem, @with_kw, @unpack, stationary_algebraic_simple, stationary_numerical_simple, irregulartrapezoidweights, create_dynamic_ODE, calculate_residuals

end # module
