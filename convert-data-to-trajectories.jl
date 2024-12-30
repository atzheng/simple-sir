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
        "--output"
            arg_type = String
            required = true
        "--gamma"
            arg_type = Float64
            default = 0.24
    end
    return parse_args(s)
end


function main()
    args = parse_commandline()
    data = DataFrame(CSV.File(args["data"]))
    dS = data[:, 2]
    I, R = dStoIR(args["gamma"], dS)
    traj_df = DataFrame(t=data[:, :week], I=I, R=R, dS=dS)
    CSV.write(args["output"], traj_df)
end

main()
