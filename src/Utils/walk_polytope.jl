export walk_polytope
function walk_polytope(f::Function, enet::EchelonMetNet, fbox::BoxGrid;
        tol = 1e-6,
        safe_int::Int = 1
    )

    d1 = 1
    
    G = enet.net.S[:, enet.idxf]
    be = enet.net.b
    lbd = enet.net.lb[enet.idxd]
    ubd = enet.net.ub[enet.idxd]

    # convex iter
    _ds = eachindex.(fbox.dims_values)
    _d1 = _ds[d1]
    fixcoors = Iterators.product(_ds[eachindex(_ds) .!= d1]...)
    _d1i0 = firstindex(_d1)
    _d1i1 = lastindex(_d1)
    
    x0k = first(fbox.dims_values[d1])
    dxk = first(fbox.dims_steps[d1])
    
    for dis in fixcoors
        
        # println("-"^10)
        vf = fbox[_d1i0, dis...]
        lbk, ubk = _vf1_extrema!(G, be, lbd, ubd, vf; k = d1)
        # @show lbk, ubk
        # lbk > ubk && continue # Not feasible
        lbki, ubki = _find_nearest.([lbk, ubk], x0k, dxk)
        # @show lbki, ubki

        # correct feasible
        _i0, _i1 = max(_d1i0, lbki - safe_int), min(_d1i1, lbki + safe_int)
        for i in _i0:1:_i1
            vf = fbox[i, dis...]
            _isfea = _isfeasible(G, be, lbd, ubd, vf; tol)
            _isfea && (lbki = i)
            _isfea && break
        end

        _i0, _i1 = max(_d1i0, ubki - safe_int), min(_d1i1, ubki + safe_int)
        for i in _i1:-1:_i0
            vf = fbox[i, dis...]
            _isfea = _isfeasible(G, be, lbd, ubd, vf; tol)
            _isfea && (ubki = i)
            _isfea && break
        end
        # @show lbki, ubki

        # iter feasibles
        for fi1 in lbki:1:ubki
            f(fbox[fi1, dis...]) === true && return nothing
        end

    end # for dis in fixcoors
    
    return nothing

end

function walk_polytope(
        f::Function, enet::EchelonMetNet, bins_or_steps;
        kwargs...
    ) 
    fbox = BoxGrid(enet, bins_or_steps)
    walk_polytope(f, enet, fbox; kwargs...)
end


## --------------------------------------------------------
# Utils

function _vf1_extrema!(G::AbstractMatrix, be::AbstractVector, 
        lbd::AbstractVector, ubd::AbstractVector,
        vf::AbstractVector; k::Int = firstindex(vf)
    )

    Nd = size(G, 1)
    vf[k] = zero(eltype(vf))
    C = be - G * vf
    lbk, ubk = -Inf, Inf
    @inbounds for i in 1:Nd
        Gik = G[i, k]
        iszero(Gik) && continue
        C1_ = inv(Gik)*(C[i] - ubd[i])
        C2_ = inv(Gik)*(C[i] - lbd[i])
        lbk = max(lbk, min(C1_, C2_))
        ubk = min(ubk, max(C1_, C2_))
    end
    
    return lbk, ubk
end

function _vf1_extrema(G::AbstractMatrix, be::AbstractVector, 
        lbd::AbstractVector, ubd::AbstractVector,
        vf::AbstractVector; k::Int = firstindex(vf)
    )
    vfk = vf[k]
    ret = _vf1_extrema!(G, be, lbd, ubd, vf; k)
    vf[k] = vfk
    return ret
end

function _find_nearest(x::Float64, x0::Float64, dx::Float64)
    i, d = divrem(x - x0, dx)
    return d < (dx / 2) ? Int(i)+1 : Int(i)+2
end

# Assumes vf already inbound
function _isfeasible(
        G::AbstractMatrix, be::AbstractVector, 
        lbd::AbstractVector, ubd::AbstractVector,
        vf::AbstractVector; tol = 1e-3
    ) 
    vd = be - G * vf
    all(lbd .- tol .<= vd .<= ubd .+ tol)
end