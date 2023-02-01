module MetXGrids

    using MetXBase
    using MetXBase: _find_nearest

    using ProgressMeter
    using Random
    using Distributions

    # TODO: create a EPGridModel (a FluxEPModelT0-like struct but based on a grid-computed distribution not a EP approximated)
    
    #! .
    
    #! include Types
    include("Types/BoxGrid.jl")
    
    #! include Api
    include("Api/entropy.jl")
    include("Api/metxbase.jl")

    #! include Utils
    include("Utils/walk_polytope.jl")

end