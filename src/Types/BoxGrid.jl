
_range(v0, v1, N::Int) = range(v0, v1; length = N)
_range(v0, v1, s::Float64) = range(v0, v1; step = s)
_range(b::Tuple, δ) = _range(b..., δ)

## ------------------------------------------------------
export BoxGrid
struct BoxGrid{D}

    dims_ids::Vector{String}
    dims_bounds::Vector{Tuple{Float64, Float64}}
    dims_δs::Vector{Union{Int, Float64}}
    dims_values::Vector{Vector{Float64}}
    dims_steps::Vector{Float64}
    idxs_iter::Iterators.ProductIterator{NTuple{D, Base.OneTo{Int64}}}
    curr_point::Vector{Float64}

    function BoxGrid(dims_ids::Vector, dim_lbs::Vector, dim_ubs::Vector, dims_δs::Vector)

        dims_bounds, dims_values, dims_steps = Tuple[], Vector{Float64}[], Float64[]
        for (id, lb, ub, δ) in zip(dims_ids, dim_lbs, dim_ubs, dims_δs)
            (lb > ub) && 
                error("Bad bounds, lb >= ub. At $(id) got ", (lb, ub))
                
            push!(dims_bounds, (lb, ub))
            r = _range(float(lb), float(ub), δ)
            push!(dims_values, collect(r))
            push!(dims_steps, step(r))
        end
        
        idxs_iter = Iterators.product(eachindex.(dims_values)...)
        D = length(dims_ids)
        curr_point = first.(dims_values)
        new{D}(
            string.(dims_ids), dims_bounds, dims_δs, 
            dims_values, dims_steps, idxs_iter, curr_point
        )
    end
end

function BoxGrid(; dims...)

    dims_ids = String[]
    dim_lbs = Float64[]
    dim_ubs = Float64[]
    dims_δs = Union{Int, Float64}[]

    for (id, (lb, ub, δ)) in dims
        push!(dims_ids, string(id))
        push!(dim_lbs, lb)
        push!(dim_ubs, ub)
        push!(dims_δs, δ)
    end

    return BoxGrid(dims_ids, dim_lbs, dim_ubs, dims_δs)
end

function Base.show(io::IO, b::BoxGrid) 
    println(io, "BoxGrid(;")
    tab = "    "
    println(io, tab, "# id = (lb, ub, δ)")
    params = String[]
    for (id, (lb, ub), δ) in zip(b.dims_ids, b.dims_bounds, b.dims_δs)
        push!(params, string(tab, id, " = ", (lb, ub, δ)))
    end
    println(io, join(params, ",\n"))
    println(io, ")")
end
Base.show(b::BoxGrid) = show(stdout, b)

## ------------------------------------------------------
# Iterator Interface
Base.IteratorSize(::BoxGrid) = Base.IteratorSize(Base.Iterators.ProductIterator)
Base.IteratorEltype(::BoxGrid) = Base.IteratorEltype(Vector{Float64})
Base.eltype(b::BoxGrid) = eltype(b.dims_values)
Base.length(b::BoxGrid) = length(b.idxs_iter)

function _iterate(b::BoxGrid, next)
    (dim_idxs, state) = next
    @inbounds for (dim, idx) in enumerate(dim_idxs)
        b.curr_point[dim] = b.dims_values[dim][idx]
    end
    return (b.curr_point, state)
end
_iterate(::BoxGrid, ::Nothing) = nothing

Base.iterate(b::BoxGrid) = _iterate(b, iterate(b.idxs_iter))
Base.iterate(b::BoxGrid, state0) = _iterate(b, iterate(b.idxs_iter, state0))

function Base.getindex(b::BoxGrid{D}, dim_idxs::Vararg{Int,D}) where {D}
    @inbounds for (dim, idx) in enumerate(dim_idxs)
        b.curr_point[dim] = b.dims_values[dim][idx]
    end
    return b.curr_point
end

function Base.collect(b::BoxGrid)
    vec = Vector{Vector{Float64}}(undef, length(b))
    @inbounds for (i, p) in enumerate(b)
        vec[i] = copy(p)
    end
    return vec
end

Base.size(b::BoxGrid) = tuple(length.(b.dims_values)...)
Base.size(b::BoxGrid, dim) = length(b.dims_values[dim])
Base.step(b::BoxGrid) = b.dims_steps
Base.step(b::BoxGrid, dim) = b.dims_steps[dimindex(b, dim)]

## ------------------------------------------------------
# rand api
function Random.rand(rng::AbstractRNG, b::BoxGrid)
    @inbounds for (dim, values) in enumerate(b.dims_values)
        b.curr_point[dim] = rand(rng, values)
    end
    return b.curr_point
end
Random.rand(b::BoxGrid) = rand(Random.GLOBAL_RNG, b)

function Random.rand(rng::AbstractRNG, b::BoxGrid, n::Int)
    vec = Vector{Vector{Float64}}(undef, n)
    @inbounds for i in 1:n
        vec[i] = copy(rand(rng, b))
    end
    return vec
end
Random.rand(b::BoxGrid, n::Int) = rand(Random.GLOBAL_RNG, b, n)

## ------------------------------------------------------
# dim api

dim_ids(b::BoxGrid) = b.dims_ids
dim_ids(b::BoxGrid, dim) = b.dims_ids[dimindex(b, dim)]

export dimindex
dimindex(b::BoxGrid, dim) = MetXBase._getindex(b, dim_ids, dim)

dim_values(b::BoxGrid, dim) = b.dims_values[dimindex(b, dim)]
