using MetXGrids
using MetXBase
using MetXNetHub
using MetXOptim
using MetXOptim: GLPK
using Test

@testset "MetXGrids.jl" begin
    
    ## --------------------------------------------------------
    # walk_polytope
    let
        net0 = pull_net("toy_net")
        
        # echelon
        net = box(net0, GLPK.Optimizer)
        enet = EchelonMetNet(net)
        tol = 1e-8

        G = enet.net.S[:, enet.idxf]
        be = enet.net.b
        lbd = enet.net.lb[enet.idxd]
        ubd = enet.net.ub[enet.idxd]
        
        bins_pool = 10:1:50
        for nbins in bins_pool
            
            @show nbins
            fbox = BoxGrid(enet, nbins)
            
            Vfea0 = 0
            Vbox0 = prod(size(fbox))
            @time walk_polytope(enet, fbox; tol) do vf
                Vfea0 += 1
            end
            @show Vbox0
            @show Vfea0
            @test Vfea0 <= Vbox0
            
            Vfea1 = 0
            Vbox1 = 0
            @time for vf in fbox
                # Test if feasible
                vd = be - G * vf
                if all(lbd .- tol .<= vd .<= ubd .+ tol)
                    Vfea1 += 1
                end
                Vbox1 += 1
            end
            @show Vbox1
            @show Vfea1
            @test Vbox0 == Vbox1
            @test Vfea0 == Vfea1
            println()
        end

    end
end
