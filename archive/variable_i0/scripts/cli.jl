# Shared command-line paths for the reproduction scripts.
function parse_cli(args; output, data=false, usage)
    root = normpath(joinpath(@__DIR__, ".."))
    names = String[]
    quick = false
    datadir = nothing
    outdir = nothing
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("--help", "-h")
            println(usage)
            exit(0)
        elseif arg == "--quick"
            quick = true
        elseif arg == "--out" || (data && arg == "--data")
            i < length(args) || error("$arg requires a directory")
            value = abspath(args[i + 1])
            arg == "--out" ? (outdir = value) : (datadir = value)
            i += 1
        elseif startswith(arg, "-")
            error("Unknown option $arg\n$usage")
        else
            push!(names, arg)
        end
        i += 1
    end
    if isnothing(datadir)
        datadir = quick ? joinpath(root, "outputs", "data", "quick") : joinpath(root, "data")
    end
    if isnothing(outdir)
        outdir = joinpath(root, "outputs", output, quick ? "quick" : "")
    end
    return (; names, quick, datadir=normpath(datadir), outdir=normpath(outdir))
end
