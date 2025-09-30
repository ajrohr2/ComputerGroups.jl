mutable struct Cache
    subgroups::Union{Vector,Nothing}
    normal_subgroups::Union{Vector,Nothing}
    center::Union{Any,Nothing}
    commutator::Union{Any,Nothing}
    derived_series::Union{Vector,Nothing}
    conjugacy_classes::Union{Vector, Nothing}
    generators::Union{Tuple, Nothing}
end

abstract type AbstractGroup end

mutable struct Group <: AbstractGroup
    elems::Vector
    id
    binary_op::Function
    order::Int
    is_abelian::Union{Bool,Nothing}
    is_solvable::Union{Bool,Nothing}
    cache::Cache
end

mutable struct FactorGroup <: AbstractGroup
    parent::Group
    factor::Group
    elems::Vector
    id
    binary_op::Function
    order::Int
    is_abelian::Union{Bool,Nothing}
    is_solvable::Union{Bool,Nothing}
    cosets::Dict
    cache::Cache
end

function Group(elements::Vector, star::Function)
    id = nothing
    for elem in elements
        found = true
        for i in elements
            if star(elem, i) != i
                found = false
                break
            end
        end
        if found
            id = elem
            break
        end
    end
    if id == nothing
        error("No identity found!")
    end
    for elem in elements
        inv = nothing
        for i in elements
            if star(elem, i) == id
                inv = i
            end
        end
        if inv == nothing
            error("$elem has no inverse!")
        end
    end
    return Group(elements, id, star, length(elements), nothing, nothing, Cache(nothing, nothing, nothing, nothing, nothing, nothing, nothing))
end
function Group(generators::Tuple, star::Function)
    elements = create_from_generators(generators, star)
    id = nothing
    for elem in elements
        found = true
        for i in elements
            if star(elem, i) != i
                found = false
                break
            end
        end
        if found
            id = elem
            break
        end
    end
    if id == nothing
        error("No identity found!")
    end
    for elem in elements
        inv = nothing
        for i in elements
            if star(elem, i) == id
                inv = i
            end
        end
        if inv == nothing
            error("$elem has no inverse!")
        end
    end
    return Group(elements, id, star, length(elements), nothing, nothing, Cache(nothing, nothing, nothing, nothing, nothing, nothing, generators))
end

function FactorGroup(parent, factor, elements, star, cosets)
    return FactorGroup(parent, factor, elements, parent.id, star, length(elements), nothing, nothing, cosets, Cache(nothing, nothing, nothing, nothing, nothing, nothing, nothing))
end

function Base.show(io::IO, G::Group)
    order = G.order
    num_subgroups = isnothing(G.cache.subgroups) ? "unknown" : length(G.cache.subgroups)
    num_norm_subgroups = isnothing(G.cache.normal_subgroups) ? "unknown" : length(G.cache.normal_subgroups)
    max_display = 3 
    print(io, "Group of order $order")
    print(io, "\n---")
    print(io, "\nElements (first $max_display): $(G.elems[1:min(3, length(G.elems))])")
    print(io, "\n$num_subgroups subgroups, $num_norm_subgroups normal")
    print(io, "\nAbelian: $(G.is_abelian)")
    print(io, "\nSolvable: $(G.is_solvable)")
    print(io, "\n---")
end
function Base.show(io::IO, G::FactorGroup)
    order = G.order
    parent_order = G.parent.order
    factor_order = G.factor.order
    num_subgroups = isnothing(G.cache.subgroups) ? "unknown" : length(G.cache.subgroups)
    num_norm_subgroups = isnothing(G.cache.normal_subgroups) ? "unknown" : length(G.cache.normal_subgroups)
    max_display = 3 
    print(io, "FactorGroup of order $order")
    print(io, "\n---")
    print(io, "\nParent group order: $parent_order, Factor group order: $factor_order")
    print(io, "\nElements (first $max_display): $(G.elems[1:min(3, length(G.elems))])")
    print(io, "\n$num_subgroups subgroups, $num_norm_subgroups normal")
    print(io, "\nAbelian: $(G.is_abelian)")
    print(io, "\nSolvable: $(G.is_solvable)")
    print(io, "\n---")
end

function element_order(a, star::Function)
    limit = 1000
    order = 1
    tmp = a
    for _ in 1:limit
        tmp = star(tmp, a)
        order += 1 
        if tmp == a
            return order - 1 
        end
    end
end
function element_power(star::Function, elem::T, n::Int) where {T}
    result = elem
    for _ in 2:n
        result = star(result, elem)
    end
    return result
end
function create_from_generators(generators::Tuple, star::Function)
    elements = Set([element_power(star, generators[1], element_order(generators[1], star))])
    queue = [element_power(star, generators[1], element_order(generators[1], star))]
    while !isempty(queue)
        current = pop!(queue)
        for g in generators
            next = star(current, g)
            if !(next ∈ elements)
                push!(elements, next)
                push!(queue, next)
            end
        end
    end
    return collect(elements)
end
