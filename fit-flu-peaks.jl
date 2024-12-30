using ArgParse
using Printf
using CSV
# using RCall

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
        "--gamma"
            arg_type = Float64
            default = 0.24
        "--seeds"
            arg_type = Int
            default = 10
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
    return Dict(:peak_t => peak_t, :peak_I => peak_I, :seed => seed)
end

function get_peaks_for_estimate(I, R, x; seeds=10, gamma=0.24)
    sir = SIR(x[:a], x[:beta], gamma, x[:N])
    results = [get_sample_peak(seed, sir, I, R) for seed in 1:seeds]
    df = DataFrame(results)
    Dict(zip(map(Symbol, names(df)), map(median, eachcol(df))))
end

function get_peaks(I, R, estimates; kwargs...)
    peaks = DataFrame([
        get_peaks_for_estimate(I[1:x[:t]], R[1:x[:t]], x; kwargs...)
        for x in estimates
            ])
    df = hcat(DataFrame(estimates), peaks)
    df[!, :actual_peak_I] .= maximum(I)
    df[!, :actual_peak_t] .= argmax(I)
    df
end

function process_id(id; root="flu/results", gamma=0.24)
    data = DataFrame(CSV.File(@sprintf("flu/data/%s.csv", id)))
    estimates = CSV.File(@sprintf("%s/%s.csv", root, id))
    dS = data[:, 2]
    I, R = dStoIR(gamma, dS)
    df = get_peaks(I, R, estimates)
    df[:id] = id
    df
end

function get_trajectory(sir, seed, I0, R0)
    traj = collect(takewhile(x -> x.I > 1,
                             sample_SIR(MersenneTwister(seed), sir;
                                        model=:poisson, I0=I0, R0=R0)))
    df = DataFrame(map(typedict, traj))
    df[:seed] = seed
    df[:t] = 1:length(traj)
    df
end

function main()
    args = parse_commandline()
    data = DataFrame(CSV.File(args["data"]))
    estimates = CSV.File(args["estimates"])
    dS = data[:, 2]
    I, R = dStoIR(args["gamma"], dS)
    df = get_peaks(
        I, R, estimates;
        seeds=args["seeds"], gamma=args["gamma"]
    )
    CSV.write(args["output"], df)
end

main()


# manifest = DataFrame(CSV.File("flu-manifest.csv", types=Dict(:id=>String)))

# results = map(process_id, tqdm(manifest[:id]));
# df = vcat(results...)
# CSV.write("flu/results.csv", df)

# results = map(partial(process_id; root="flu/true-N"), tqdm(manifest[:id]));
# df = vcat(results...)
# CSV.write("flu/true-N-results.csv", df)


# # Spot Check fits
# # ---------------------------------------------------------------------
# id = "2010.01"
# data = DataFrame(CSV.File(@sprintf("flu/data/%s.csv", id)))
# dS = data[:, 2] .- minimum(dS)
# I, R = dStoIR(GAMMA, dS)
# estimates = CSV.File(@sprintf("flu/results-v3/%s.csv", id))
# x = estimates[end]
# # sir = SIR(x[:a], x[:beta], GAMMA, x[:N])
# sir = SIR(x[:a], x[:beta], GAMMA, x[:N])
# t0 = 4
# trajs = vcat([get_trajectory(sir, seed, I[t0], R[t0]) for seed in 1:50]...)
# trajs[:t] .= trajs[:t] .+ t0
# actual = DataFrame(I=I, R=R, t=1:length(I), dataset="actual")
# x  = fit(I[t0:end], R[t0:end], dS[t0:end], 1e6; gamma=0.1)

# R"""
# ($trajs
#     %>% filter(t <= 60)
#     %>% group_by(t) %>% summarise(I=mean(I), R=mean(R))
#     %>% mutate(dataset="predicted" )
#     %>% rbind($actual)
#     %>% group_by(dataset)
#     %>% mutate(C=I+R, dC=C-dplyr::lag(C))
#     %>% gather(var, val, -t, -dataset)
#     %>% ggplot(aes(t, val, color=dataset))
#     + facet_wrap(~ var, scales='free')
#     + geom_line())
# """
