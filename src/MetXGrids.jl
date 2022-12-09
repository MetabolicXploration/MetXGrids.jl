module MetXGrids

    using MetXBase

    using ProgressMeter
    using Base.Threads
    using ExtractMacro
    using Statistics
    using Random
    
    #! .
    
    #! include Types
    include("Types/BoxGrid.jl")
    
    #! include Api
    include("Api/metxbase.jl")

    #! include Utils
    include("Utils/walk_polytope.jl")
    # include("Utils/utils.jl")

end