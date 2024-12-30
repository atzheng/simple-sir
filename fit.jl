using ArgParse
using CSV
using Printf
include("sir.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input"
            arg_type = String
            required = true
        "--output"
            arg_type = String
            required = true
        "--model_type"
            arg_type = Symbol
            required = true
        "--gamma"
            arg_type = Float64
            default = 0.1
        "--samples"
            arg_type = Int
            default = 100
        "--N"
            arg_type = Float64
            default = 0.
        "--tol"
            arg_type = Float64
            default = 1e-12
    end
    return parse_args(s)
end

function sample_idx(dS, samples)
    if samples >= length(dS) return collect(1:length(dS)) end
    progress = log.(cumsum(dS)) ./ log(sum(dS))
    breaks = reshape(collect(LinRange(0.4, 1, samples)), 1, samples)
    map(rpartial(getindex, 1),
        vec(argmax(progress .>= breaks; dims=1))) |> unique
end

function fit(I, R, X, Nmax; model_type=:sir, gamma=0., N=0.)
    "Fit an SIR or Bass model to the data."
    kwargs = if model_type == :bass Dict(:g => 0.)
    elseif model_type == :sir Dict(:a => 0., :g => gamma) end
    fit_given_N(N) =
        SIR_ipopt_mle(
            I, R, X;
            model_type=:poisson,
            N=N,
            kwargs...
                )
    if N == 0.  # Fit with N unknown
        opt_N = optimize(
            N -> - fit_given_N(N)[:llh],
            maximum(I .+ R), Nmax, GoldenSection()
        ).minimizer
    else
        opt_N = N
    end
    fit_given_N(opt_N)
end

function main()
    args = parse_commandline()
    df = DataFrame(CSV.File(args["input"]))
    dS = df[:, 2]
    Nmax = maximum(df[:, 3])
    I, R = dStoIR(args["gamma"], df[:, 2])

    # llh(a, b, g) = sum(
    # dS[t] * log((a + b * I[t]) * (1 - (I[t] + R[t]) / N) + g * I[t])
    # - (a + b * I[t]) * (1 - (I[t] + R[t]) / N) + g * I[t]
    # )

    print(args["gamma"])
    ts = if args["samples"] > 1
        filter(relates(>, 1), sample_idx(dS, args["samples"]))
    else
        [length(dS)] # Fit the entire trajectory
    end
    results = [
        fit(I[1:t-1], R[1:t-1], dS[2:t], Nmax;
            model_type=args["model_type"],
            gamma=args["gamma"],
            N=args["N"])
        for t in ts
            ]
    df = DataFrame(results)
    df[!, :t] = ts
    df[!, :m] = cumsum(dS)[ts]
    CSV.write(args["output"], df)
end

main()
