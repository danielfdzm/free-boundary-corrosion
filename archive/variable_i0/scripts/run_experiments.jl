#!/usr/bin/env julia
# Run the experiments reported in the paper; write fresh data under outputs/.
using FreeBoundaryNumerics, Printf
include(joinpath(@__DIR__, "cli.jl"))

function main(args)
    options = parse_cli(args; output="data",
        usage="julia --project=. scripts/run_experiments.jl [E1 E2 E3 E4 E5] [--quick] [--out DIR]")
    runners = Dict("E1" => run_E1, "E2" => run_E2, "E3" => run_E3,
        "E4" => run_E4, "E5" => run_E5)
    todo = isempty(options.names) ? ["E1", "E2", "E3", "E4", "E5"] : options.names
    for name in todo
        haskey(runners, name) || error("Unknown experiment $name")
    end
    mkpath(options.outdir)
    for name in todo
        println("$name (quick=$(options.quick)) -> $(options.outdir)")
        started = time()
        runners[name](options.outdir; quick=options.quick)
        @printf("%s finished in %.1f s\n", name, time() - started)
    end
end

main(ARGS)
