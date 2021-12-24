module DiffEqBayesTuring

    using DocStringExtensions
    using DiffEqBase, Distributions, MacroTools
    using RecursiveArrayTools, ModelingToolkit
    using Parameters, Distributions, Optim, Requires
    using Distances, DocStringExtensions, Random
    using Turing, MCMCChains

    STANDARD_PROB_GENERATOR(prob,p) = remake(prob;u0=eltype(p).(prob.u0),p=p)
    STANDARD_PROB_GENERATOR(prob::EnsembleProblem,p) = EnsembleProblem(remake(prob.prob;u0=eltype(p).(prob.prob.u0),p=p))

    include("turing_inference.jl")

    const src_path = @__DIR__

    """

    # debt_path

    Relative path using the DiffEqBayesTuring src/ directory.

    ### Example to get access to the data subdirectory
    ```julia
    debt_path("..", "data")
    ```

    Note that in the projects, e.g. StatisticalRethinkingStan.jl and StatisticalRethinkingTuring.jl, the
    DrWatson approach is a better choics, i.e: `sr_datadir(filename)`

    """
    debt_path(parts...) = normpath(joinpath(src_path, parts...))

    # DrWatson extension
    """

    # debs_datadir

    Relative path using the StatisticalRethinking src/ directory.

    ### Example to access `Howell1.csv` in StatisticalRethinking:
    ```julia
    df = CSV.read(sr_datadir("Howell1.csv"), DataFrame)
    ```
    """
    debt_datadir(parts...) = debt_path("..", "data", parts...)


    export
        turing_inference,
        debt_path,
        debt_datadir
    

end
