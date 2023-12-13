using ArgParse
using Printf
using CSV
include("sir.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--data"
            arg_type = String
            required = true
        "--estimates"
            arg_type = String
            required = true
        "--output"
            arg_type = String
            required = true
        "--seeds"
            arg_type = Int
            default = 30
    end
    return parse_args(s)
end


function get_sample_peak(seed, sir, I, R)
    generator = sample_SIR(
        MersenneTwister(seed), sir;
        I0=I[end], R0=R[end], model=:poisson
    )
    traj = collect(takewhile(x -> x.I > 1, generator))
    Isim = cat(I, map(rpartial(getfield, :I), traj); dims=1)
    Tsim = cat(1:length(I),
               length(I) .+ cumsum(map(rpartial(getfield, :T), traj));
               dims=1)
    peak_k = argmax(Isim)
    peak_I = Isim[peak_k]
    peak_t = Tsim[peak_k]
    return Dict(:peak_t => peak_t, :peak_I => peak_I)
end

function get_peaks_for_estimate(I, R, x; seeds=10)
    sir = SIR(x[:a], x[:beta], GAMMA, x[:N])
    results = [get_sample_peak(seed, sir, I, R) for seed in 1:seeds]
    df = DataFrame(results)
    Dict(zip(map(Symbol, names(df)), map(mean, eachcol(df))))
end

function get_peaks(I, R, estimates)
    peaks = DataFrame([get_peaks_for_estimate(I[1:x[:t]], R[1:x[:t]], x)
                       for x in estimates])
    df = hcat(DataFrame(estimates), peaks)
    df[!, :actual_peak_I] = maximum(I)
    df[!, :actual_peak_t] = argmax(I)
    df
end

function get_trajectory(sir, seed, I0, R0)
    traj = collect(
        takewhile(
            x -> x.I > 1,
            sample_SIR(
                MersenneTwister(seed), sir;
                model=:poisson, I0=I0, R0=R0
            )
        )
    )
    df = DataFrame(map(typedict, traj))
    df[!, :seed] .= seed
    df[!, :t] .= 1:length(traj)
    df
end

function main()
    args = parse_commandline()
    data = DataFrame(CSV.File(args["data"]))
    df = DataFrame(CSV.File(args["estimates"]))
    sort!(df, [:t])
    N, beta, gamma = df[end, [:N, :beta, :gamma]]
    I, R = dStoIR(gamma, data[!, :dS])
    I0, R0 = I[1], R[1]
    print("N = $N, beta = $beta, gamma = $gamma, I0 = $I0, R0 = $R0 \n")
    sir = SIR(0., beta, gamma, N)
    results = vcat([get_trajectory(sir, seed, I0, R0)
                    for seed in 1:args["seeds"]]...)
    CSV.write(args["output"], results)
end

main()
