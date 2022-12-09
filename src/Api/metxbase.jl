function BoxGrid(enet::EchelonMetNet, bins_or_steps)

    net = metnet(enet)
    dims_ids = reactions(net, enet.idxf)
    dim_lbs, dim_ubs = bounds(net, enet.idxf)
    bins_or_steps = (bins_or_steps isa Number) ? fill(bins_or_steps, length(enet.idxf)) : bins_or_steps

    return BoxGrid(dims_ids, dim_lbs, dim_ubs, bins_or_steps)

end