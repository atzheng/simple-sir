using LinearAlgebra
using Random
using Statistics
using StatsBase
using Distributions
using DataFrames
using ProgressBars
using Functional
using Lazy
using JuMP
using Ipopt
using Setfield
using Optim
using Zygote


# Utils
# -------------------------------------------------------------------------
typedict(x) = Dict(fn => getfield(x, fn) for fn ∈ fieldnames(typeof(x)))

logrange(T, samples) =
    (samples >= T ? collect(1:T) :
     @> 0.:log(T)/samples:log(T) collect exp.() floor.() unique Array{Int64, 1})

function cummean(x; dims=1)
    shape = ones(Int32, length(size(x)))
    shape[dims] = size(x, dims)
    idx = reshape(collect(1:size(x, dims)), tuple(shape...))
    cumsum(x; dims=dims) ./ idx
end

# Simulation
# -------------------------------------------------------------------------
struct SIR
    alpha :: Float64
    beta  :: Float64
    gamma :: Float64
    N     :: Float64
end

struct SIRState
    I :: Float64
    R :: Float64
    T :: Float64 # Last interarrival time
end

is_SIR_terminated(sir, s) = (s.I <= 0) | (s.I + s.R >= sir.N)

function poisson_SIR_step(rng, sir, s::SIRState)
    if is_SIR_terminated(sir, s) return nothing end
    S = sir.N - s.I - s.R
    lam = (s.I .* sir.beta .+ sir.alpha) .* S ./ sir.N
    arrs = rand(rng, Poisson(lam))
    SIRState((1 - sir.gamma) * s.I + arrs,
             sir.gamma * s.I + s.R,
             1)
end

function SIR_interarrival_distr(sir::SIR, s::SIRState)
    S = sir.N .- s.I .- s.R
    lam = (s.I .* sir.beta .+ sir.alpha) .* S ./ sir.N .+ sir.gamma * s.I
    theta = ifelse.(1 ./ lam .> 0, 1 ./ lam, Inf)
    Exponential(theta)
end

function SIR_step(rng, sir::SIR, s::SIRState)
    if is_SIR_terminated(sir, s) return nothing end
    S = sir.N - s.I - s.R
    infection_rate = S * (sir.beta * s.I + sir.alpha)
    p = infection_rate / (infection_rate + sir.N * sir.gamma * s.I)
    is_infection = rand(rng, Bernoulli(p))
    T = rand(rng, SIR_interarrival_distr(sir, s))
    SIRState(s.I + ifelse(is_infection, 1., -1.),
             s.R + ifelse(is_infection, 0., 1.), T)
end

steppers = Dict(:exponential => SIR_step,
                :poisson => poisson_SIR_step)

function sample_SIR(rng::AbstractRNG, sir::SIR;
                    I0=1., R0=0., model=:exponential)
    drop(takeuntil(
        partial(is_SIR_terminated, sir),
        iterated(partial(steppers[model], rng, sir),
                 SIRState(I0, R0, Inf))), 1)
end

# Estimation
# -------------------------------------------------------------------------
SIR_llh(sir::SIR, states) =
    logpdf.(map(partial(SIR_interarrival_distr, sir),
                states[1:end-1]),
            map(rpartial(getfield, :T), states[2:end]))

function grid_search_mle(llh, grid)
    "Finds the maximum of a univariate function llh over a grid."
    llhs = hcat(map(llh, grid)...)
    best_idx = argmax(cummean(llhs; dims=1); dims=2)
    estimate = vec(grid[map(rpartial(getindex, 2), best_idx)])
    DataFrame(k=1:size(llhs, 1),
              estimate=estimate,
              is_ub_binding=estimate .>= maximum(grid) - 1e-6,
              is_lb_binding=estimate .<= minimum(grid) - 1e-6,
              llh=vec(llhs[best_idx]))
end

function SIR_grid_search_mle(sir, param, states;
                             theta_max=getfield(sir, param) * 10,
                             theta_min=0., breaks=1000)
    lens = Setfield.PropertyLens{param}()
    llh = theta -> SIR_llh(set(sir, lens, theta), states)
    grid = LinRange(theta_min, theta_max, breaks)
    grid_search_mle(llh, grid)
end

result_to_dict(r) = Dict(
    :estimate => r.minimizer,
    :converged => r.converged,
    :llh => - r.minimum,
    :is_ub_binding => r.minimizer >= r.initial_upper - (r.initial_upper - r.initial_lower) / 1e5,
    :is_lb_binding => r.minimizer <= r.initial_lower + (r.initial_upper - r.initial_lower) / 1e5)

function SIR_golden_section_mle(sir, param, states;
                                samples=length(states),
                                theta_min=0.,
                                theta_max=getfield(sir, param) * 10)
    Ts = filter(x -> (x > 10) & (x < length(states)),
                logrange(sir.N, samples))
    lens = Setfield.PropertyLens{param}()
    llh(t, theta) = - sum(SIR_llh(set(sir, lens, theta), states[1:t]))
    df = DataFrame([result_to_dict(
        optimize(partial(llh, t), theta_min, theta_max, GoldenSection()))
               for t in Ts])
    df[:k] = Ts
    df
end

function SIR_ipopt_mle(data; a=nothing, b=nothing, g=nothing, N=nothing,
                       Nmax=1e5)
    I = map(rpartial(getfield, :I), data)
    R = map(rpartial(getfield, :R), data)
    T = map(rpartial(getfield, :T), data)
    m = Model(Ipopt.Optimizer)
    if a == nothing @variable(m, a >= 0) end
    if b == nothing @variable(m, b >= 0) end
    if g == nothing @variable(m, g >= 0) end
    if N == nothing
        @variable(m, N >= maximum(I .+ R) + 1)
        @constraint(m, N <= Nmax)
    end
    @NLobjective(
        m, Max,
        sum(log((a + b * I[t]) * (1 - (I[t] + R[t]) / N) + g * I[t])
            - T[t] * (a + b * I[t]) * (1 - (I[t] + R[t]) / N) + g * I[t]
            for t in 1:length(data)))
    optimize!(m)
    Dict(:a => getvalue(a),
         :beta => getvalue(b),
         :gamma => getvalue(g),
         :N => getvalue(N),
         :llh => objective_value(m),
         :status => Int(termination_status(m)))
end


fakeSIRStep(gamma, IR, dS) =
    ((1 - gamma) * IR[1] + dS, IR[2] + gamma * IR[1])
function dStoIR(gamma, dS)
    IRs = accumulate(partial(fakeSIRStep, gamma), dS; init=(0., 0.))
    [map(rpartial(getindex, i), IRs) for i in 1:2]
end

function SIR_ipopt_mle(I, R, X; model_type=:exponential,
                       a=nothing, b=nothing, g=nothing,
                       N=nothing, Nmax=1e5)
    m = Model(Ipopt.Optimizer)
    if a == nothing @variable(m, a >= 0) end
    if b == nothing @variable(m, b >= 0) end
    if g == nothing @variable(m, g >= 0) end
    if N == nothing
        @variable(m, N >= maximum(I .+ R) + 1)
        @constraint(m, N <= Nmax)
    end
    if model_type == :exponential
        W1, W2 = (ones(length(X)), X)
    elseif model_type == :poisson
        W1, W2 = (X, ones(length(X)))
    end
    @NLobjective(
        m, Max,
        sum(W1[t] * log((a + b * I[t]) * (1 - (I[t] + R[t]) / N) + g * I[t])
            - W2[t] * (a + b * I[t]) * (1 - (I[t] + R[t]) / N) + g * I[t]
            for t in 1:length(X)))
    optimize!(m)
    Dict(:a => getvalue(a),
         :beta => getvalue(b),
         :gamma => getvalue(g),
         :N => getvalue(N),
         :llh => objective_value(m),
         :status => Int(termination_status(m)))
end

function interior_point_mle(f, x0)
    g!(G, x) = (G .= gradient(f, x)[1])
    H!(H, x) = (H .= hessian(f, x))
    obj = TwiceDifferentiable(f, g!, H!, x0)
    con = TwiceDifferentiableConstraints(
        zeros(length(x0)), Inf .* ones(length(x0)))
    optimize(obj, con, x0, IPNewton())
end


function run_experiment(sir, param, seed; kwargs...)
    rng = MersenneTwister(seed)
    data = collect(sample_SIR(rng, sir))
    # df = SIR_grid_search_mle(sir, param, data; kwargs...)
    df = SIR_golden_section_mle(sir, param, data; kwargs...)
    df[:seed] = seed
    df
end

run_experiments(sir, param, seeds; kwargs...) =
    vcat(filter(x -> size(x, 1) > 0,
                [run_experiment(sir, param, seed; kwargs...)
                 for seed in tqdm(1:seeds)])...)
