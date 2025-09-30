function powerset_inplace(s::Vector{T}, n::Int, storage::Vector{T}) where {T}
    offset = 0
    itr = 1 
    while n != 0
        tz = trailing_zeros(n)
        pos = offset + tz + 1 
        if pos <= length(s)
            storage[itr] = s[pos]
        else
            break
        end
        n >>= (tz+1)
        offset += tz + 1 
        itr += 1 
    end
end
function find_subgroups(G::AbstractGroup)
    sgs = Set([Set([G.id])])
    powers = sort(unique([sum(collect(values(Dict(i for i in factor(num))))) for num in divisors(G.order)]))[2:end-1]
    elem_orders = Dict(a => element_order(a, G.binary_op) for a in G.elems)
    elem_storage = Vector{typeof(G.id)}(undef, 2)
    subset_storage = Vector{typeof(G.id)}(undef, length(G.elems))
    for pow in powers
        # if !isempty(filter(i -> length(i)==G.order, sgs))
        #     break
        # end
        gens = powerset(G.elems, pow, pow)
        for s in gens
            skip = false
            for i in eachindex(s)
                for j in i+1:length(s)-1
                    if elem_orders[s[j]] == elem_orders[s[i]]
                        continue
                    end
                    if max(elem_orders[s[j]],elem_orders[s[i]]) % min(elem_orders[s[j]],elem_orders[s[i]]) == 0
                        elem_storage[1] = s[j]
                        elem_storage[2] = s[i]
                        sort!(elem_storage, by=(i -> elem_orders[i]))
                        for n in 1:elem_orders[elem_storage[1]]
                            if element_power(G.binary_op, elem_storage[1], n) == elem_storage[2]
                                skip = true
                                break
                            end
                        end
                    end
                end
            end
            if skip
                continue
            end
            sg = Set(create_from_generators(ntuple(i -> s[i], length(s)), G.binary_op))
            push!(sgs, sg)
        end
    end
    G.cache.subgroups = [Group(collect(sg), G.binary_op) for sg in sgs]
    return G.cache.subgroups
end

function find_inverse(elem, G::AbstractGroup)
    for e in G.elems
        if G.binary_op(elem, e) == G.id
            return e
        end
    end
end

function dict_to_cycle(d::Dict)
    n = maximum(keys(d))
    visited = falses(n)
    cycles = Vector{Vector{Int64}}()
    for start in 1:n
        if !visited[start]
            cycle = Int64[]
            current = start
            while !visited[current]
                push!(cycle, current)
                visited[current] = true
                current = d[current]
            end
            if length(cycle) > 1
                push!(cycles, cycle)
            end
        end
    end
    return isempty(cycles) ? "()" : join(["(" * join(cycle, " ") * ")" for cycle in cycles], "")
end
function dict_to_cycle(dv)
    new_dv = []
    for d in dv
        n = maximum(keys(d))
        visited = falses(n)
        cycles = Vector{Vector{Int64}}()
        for start in 1:n
            if !visited[start]
                cycle = Int64[]
                current = start
                while !visited[current]
                    push!(cycle, current)
                    visited[current] = true
                    current = d[current]
                end
                if length(cycle) > 1
                    push!(cycles, cycle)
                end
            end
        end
        push!(new_dv, isempty(cycles) ? "()" : join(["(" * join(cycle, " ") * ")" for cycle in cycles], ""))
    end
    return new_dv
end

function is_normal(G::AbstractGroup, H::AbstractGroup)
    for elem in G.elems
        if !issetequal(G.binary_op.(G.binary_op.(Ref(elem), H.elems), Ref(find_inverse(elem,G))), H.elems)
            return false
        end
    end
    return true
end

function find_normal_subgroups(G::AbstractGroup)
    normal_subgroups = []
    if isnothing(G.cache.subgroups)
        find_subgroups(G)
    end
    for sg in G.cache.subgroups
        if is_normal(G, sg)
            push!(normal_subgroups, sg)
        end
    end
    G.cache.normal_subgroups = normal_subgroups
    return G.cache.normal_subgroups
end

function is_abelian(G::AbstractGroup)
    for a in G.elems, b in G.elems
        if G.binary_op(a, b) != G.binary_op(b, a)
            G.is_abelian = false
            return false
        end
    end
    G.is_abelian = true
    return true
end

function factor_group_binary_op(a, b, parent_binary_op, N::AbstractGroup, cosets::Dict)
    result = parent_binary_op(a, b)
    if result ∈ N.elems
        return N.id
    elseif result ∉ collect(values(cosets))
        for key in keys(cosets)
            if result ∈ key
                return cosets[key]
            end
        end
    else
        return result
    end
end
function find_equivalences(G::AbstractGroup, N::AbstractGroup)
    elems = setdiff(G.elems, N.elems)
    m = length(elems)
    for idx in 1:m
        setdiff!(elems, G.binary_op.(Ref(elems[idx]), setdiff(N.elems, [N.id])))
        if idx ≥ length(elems)
            break
        end
    end
    return elems
end
function find_factor_group(G::AbstractGroup, N::AbstractGroup)
    @assert is_normal(G,N) "Subgroup must be normal in G"
    elems = [G.id] ∪ find_equivalences(G, N)
    cosets = Dict(G.binary_op.(Ref(i), N.elems) => i for i in elems)
    factor_binary_op = (a,b) -> factor_group_binary_op(a, b, G.binary_op, N, cosets)
    return FactorGroup(G, N, elems, factor_binary_op, cosets)
end

function find_center(G::AbstractGroup)
    if isnothing(G.is_abelian)
        is_abelian(G)
    end
    if G.is_abelian
        G.cache.center = Group(G.elems, G.binary_op)
    else
        center = Dict{Any, Any}[]
        for a in G.elems
            for b in G.elems
                if G.binary_op(a,b) != G.binary_op(b,a)
                    break
                end
                if findfirst(==(b), G.elems) == length(G.elems)
                    push!(center, a)
                end
            end
        end
        G.cache.center = Group(center, G.binary_op)
    end
    return G.cache.center
end

function find_commutator(G::AbstractGroup)
    if isnothing(G.is_abelian)
        is_abelian(G)
    end
    if G.is_abelian
        G.cache.commutator = Group([G.id], G.binary_op)
    else
        commutators = Dict{Any, Any}[]
        for a in G.elems, b in G.elems
            ainv = find_inverse(a, G)
            binv = find_inverse(b, G)
            push!(commutators, G.binary_op(G.binary_op(G.binary_op(ainv, binv), a), b))
        end
        unique!(commutators)
        G.cache.commutator = Group(commutators, G.binary_op)
    end
    return G.cache.commutator
end

function find_derived_series(G::AbstractGroup)
    series = AbstractGroup[G]
    for i in 1:G.order
        push!(series, find_commutator(series[i]))
    end
    G.cache.derived_series = unique!(x -> x.order, series)
    return G.cache.derived_series
end
function is_solvable(G::AbstractGroup)
    if isnothing(G.cache.derived_series)
        find_derived_series(G)
    end
    if !isempty(filter(x -> x.order == 1, G.cache.derived_series))
        G.is_solvable = true
        return true
    else
        G.is_solvable = false
        return false
    end
end

function find_inn(G::AbstractGroup)
    if isnothing(G.cache.center)
        find_center(G)
    end
    return find_factor_group(G, G.cache.center)
end

function is_isomorphic(G::AbstractGroup, H::AbstractGroup)
    if isnothing(G.is_abelian)
        is_abelian(G)
    end
    if isnothing(H.is_abelian)
        is_abelian(H)
    end
    if isnothing(G.is_solvable)
        is_solvable(G)
    end
    if isnothing(H.is_solvable)
        is_solvable(H)
    end
    if isnothing(G.cache.derived_series)
        find_derived_series(G)
    end
    if isnothing(H.cache.derived_series)
        find_derived_series(H)
    end
    if G.order != H.order
        return false
    elseif G.is_abelian != H.is_abelian
        return false
    elseif G.is_solvable != H.is_solvable
        return false
    elseif length(G.cache.derived_series) != length(H.cache.derived_series)
        return false
    end

    # All checks pass, now move to creating mappings between same order elements
    Gelemsorders = Dict(i => element_order(i, G.binary_op) for i in G.elems)
    Helemsorders = Dict(i => element_order(i, H.binary_op) for i in H.elems)

        if !issetequal(collect(values(Gelemsorders)), collect(values(Helemsorders)))
        return false
    end

    function try_bijection(G_elems, H_elems, mapping=Dict(), used=Set())
        if length(mapping) == G.order
            for g1 in G_elems, g2 in G_elems
                if H.binary_op(mapping[g1], mapping[g2]) != mapping[G.binary_op(g1, g2)]
                    return false, nothing
                end
            end
            return true, mapping
        end

        g = G_elems[length(mapping)+1]
        for h in H_elems
            if h ∉ used && Gelemsorders[g] == Helemsorders[h]
                push!(used, h)
                mapping[g] = h
                success, result = try_bijection(G_elems, H_elems, mapping, used)
                if success
                    return true, result
                end
                delete!(mapping, g)
                pop!(used)
            end
        end
        return false, nothing
    end

    success, bijection = try_bijection(G.elems, H.elems)
    return success, bijection
end

function conjugacy_classes(G::AbstractGroup)
    unused = copy(G.elems)
    if isnothing(G.cache.conjugacy_classes)
        G.cache.conjugacy_classes = []
    else
        return G.cache.conjugacy_classes
    end
    for a ∈ G.elems
        if a ∉ unused
            continue
        end
        conj_class = [reduce(G.binary_op, [g, a, find_inverse(g, G)]) for g in G.elems]
        push!(G.cache.conjugacy_classes, unique!(conj_class))
        setdiff!(unused, conj_class)
    end
    return G.cache.conjugacy_classes
end
