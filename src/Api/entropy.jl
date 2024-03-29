import Distributions.entropy
function entropy(f::Function, enet::EchelonMetNet, fbox::BoxGrid)
    
    Σfi = 0.0
    Σfi_logfi = 0.0
    
    walk_polytope(enet, fbox) do vf
        fi = f(vf)
        Σfi += fi
        # @show fi
        Σfi_logfi += fi * log(fi)
    end
    # @show Σfi
    # @show Σfi_logfi
    ΔV = prod(fbox.dims_steps)
    Σfi_ΔV = Σfi * ΔV
    Σfi_logfi_ΔV = Σfi_logfi * ΔV
    
    return log(Σfi_ΔV) - inv(Σfi_ΔV) * Σfi_logfi_ΔV
end

_one(x) = one(eltype(x))
entropy(enet::EchelonMetNet, fbox::BoxGrid) = 
    entropy(_one, enet, fbox)