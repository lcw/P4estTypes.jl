@inline pxest_t(::Val{4}) = p4est_t
@inline pxest_t(::Val{8}) = p8est_t
@inline pxest_new_ext(::Val{4}) = p4est_new_ext
@inline pxest_new_ext(::Val{8}) = p8est_new_ext
@inline pxest_balance_ext(::Val{4}) = p4est_balance_ext
@inline pxest_balance_ext(::Val{8}) = p8est_balance_ext
@inline pxest_coarsen_ext(::Val{4}) = p4est_coarsen_ext
@inline pxest_coarsen_ext(::Val{8}) = p8est_coarsen_ext
@inline pxest_refine_ext(::Val{4}) = p4est_refine_ext
@inline pxest_refine_ext(::Val{8}) = p8est_refine_ext
@inline pxest_partition_ext(::Val{4}) = p4est_partition_ext
@inline pxest_partition_ext(::Val{8}) = p8est_partition_ext
@inline pxest_partition_lnodes(::Val{4}) = p4est_partition_lnodes
@inline pxest_partition_lnodes(::Val{8}) = p8est_partition_lnodes
@inline pxest_ghost_t(::Val{4}) = p4est_ghost_t
@inline pxest_ghost_t(::Val{8}) = p8est_ghost_t
@inline pxest_ghost_new(::Val{4}) = p4est_ghost_new
@inline pxest_ghost_new(::Val{8}) = p8est_ghost_new
@inline pxest_lnodes_t(::Val{4}) = p4est_lnodes
@inline pxest_lnodes_t(::Val{8}) = p8est_lnodes
@inline pxest_lnodes_new(::Val{4}) = p4est_lnodes_new
@inline pxest_lnodes_new(::Val{8}) = p8est_lnodes_new
@inline pxest_quadrant_t(::Val{4}) = p4est_quadrant
@inline pxest_quadrant_t(::Val{8}) = p8est_quadrant
@inline pxest_tree_t(::Val{4}) = p4est_tree
@inline pxest_tree_t(::Val{8}) = p8est_tree
@inline pxest_iter_volume_info_t(::Val{4}) = p4est_iter_volume_info_t
@inline pxest_iter_volume_info_t(::Val{8}) = p8est_iter_volume_info_t

@inline CONNECT_FULL(::Val{4}) = P4EST_CONNECT_FULL
@inline CONNECT_FULL(::Val{8}) = P8EST_CONNECT_FULL
@inline CONNECT_FACE(::Val{4}) = P4EST_CONNECT_FACE
@inline CONNECT_FACE(::Val{8}) = P8EST_CONNECT_FACE
@inline CONNECT_CORNER(::Val{4}) = P4EST_CONNECT_CORNER
@inline CONNECT_CORNER(::Val{8}) = P8EST_CONNECT_CORNER
@inline CONNECT_EDGE(::Val{8}) = P8EST_CONNECT_CORNER

const Locidx = P4est.p4est_locidx_t
const Gloidx = P4est.p4est_gloidx_t

# X = 4 or 8
# T = user data type
struct Quadrant{X,T,P}
    pointer::P
end

@inline level(quadrant::Quadrant) =
    GC.@preserve quadrant unsafe_load(quadrant.pointer).level
@inline function coordinates(quadrant::Quadrant{4})
    GC.@preserve quadrant begin
        qs = unsafe_load(quadrant.pointer)
        return (qs.x, qs.y)
    end
end
@inline function unsafe_which_tree(quadrant::Quadrant)
    return GC.@preserve quadrant unsafe_load(quadrant.pointer).p.piggy3.which_tree + 0x1
end
@inline function coordinates(quadrant::Quadrant{8})
    GC.@preserve quadrant begin
        qs = unsafe_load(quadrant.pointer)
        return (qs.x, qs.y, qs.z)
    end
end
function storeuserdata!(quadrant::Quadrant{X,T}, data::T) where {X,T}
    GC.@preserve quadrant begin
        qs = unsafe_load(quadrant.pointer)
        unsafe_store!(Ptr{T}(qs.p.user_data), data)
    end
end
function loaduserdata(quadrant::Quadrant{X,T}) where {X,T}
    GC.@preserve quadrant begin
        unsafe_load(Ptr{T}(unsafe_load(quadrant.pointer).p.user_data))
    end
end

struct Tree{X,T,P,Q} <: AbstractArray{Quadrant,1}
    pointer::P
    forest::Q
end

Base.size(t::Tree) =
    (convert(Int, (GC.@preserve t unsafe_load(t.pointer).quadrants.elem_count)),)
function Base.getindex(t::Tree{X,T}, i::Int) where {X,T}
    @boundscheck checkbounds(t, i)
    GC.@preserve t begin
        Q = pxest_quadrant_t(Val(X))
        quadrant = Ptr{Q}(unsafe_load(t.pointer).quadrants.array + sizeof(Q) * (i - 1))
        return Quadrant{X,T,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Tree) = IndexLinear()

offset(tree::Tree) = GC.@preserve tree unsafe_load(tree.pointer).quadrants_offset

mutable struct Pxest{X,T,P,C} <: AbstractArray{Tree,1}
    pointer::P
    connectivity::C
    comm::MPI.Comm
    function Pxest{4}(
        pointer::Ptr{p4est_t},
        connectivity::Connectivity{4},
        comm::MPI.Comm,
        ::Type{T},
    ) where {T}
        forest = new{4,T,typeof(pointer),typeof(connectivity)}(pointer, connectivity, comm)
        finalizer(forest) do p
            p4est_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
    function Pxest{8}(
        pointer::Ptr{p8est_t},
        connectivity::Connectivity{8},
        comm::MPI.Comm,
        ::Type{T},
    ) where {T}
        forest = new{8,T,typeof(pointer),typeof(connectivity)}(pointer, connectivity, comm)
        finalizer(forest) do p
            p8est_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
end

function unsafe_get_user_pointer_pointer(forest::Pxest{X,T,P}) where {X,T,P}
    return unsafe_get_user_pointer_pointer(forest.pointer)
end

function unsafe_get_user_pointer_pointer(ptr::Ptr{F}) where {F<:Union{p4est_t,p8est_t}}
    i = findfirst(isequal(:user_pointer), fieldnames(F))
    offset = fieldoffset(F, i)
    ptrtype = Ptr{fieldtype(F, i)}
    return reinterpret(ptrtype, ptr + offset)
end

function pxest(
    connectivity::Connectivity{X};
    comm = MPI.COMM_WORLD,
    min_quadrants = 0,
    min_level = 0,
    fill_uniform = true,
    data_type = Nothing,
    init_function = nothing,
) where {X}
    MPI.Initialized() || MPI.Init()

    pointer = (pxest_new_ext(Val(X)))(
        comm,
        connectivity,
        min_quadrants,
        min_level,
        fill_uniform,
        sizeof(data_type),
        C_NULL,
        C_NULL,
    )

    forest = Pxest{X}(pointer, connectivity, comm, data_type)

    if !isnothing(init_function)
        init(forest, _, quadrant, _, treeid, _) = init_function(forest, treeid, quadrant)
        iterateforest(forest; volume = init)
    end

    return forest
end

quadrantstyle(::Pxest{X}) where {X} = X
quadrantndims(::Pxest{4}) = 2
quadrantndims(::Pxest{8}) = 3
typeofquadrantuserdata(::Pxest{X,T}) where {X,T} = T
lengthoflocalquadrants(p::Pxest) =
    GC.@preserve p PointerWrapper(p.pointer).local_num_quadrants[]
lengthofglobalquadrants(p::Pxest) =
    GC.@preserve p PointerWrapper(p.pointer).global_num_quadrants[]
comm(p::Pxest) = p.comm
connectivity(p::Pxest) = p.connectivity

function Base.unsafe_convert(::Type{Ptr{p4est_t}}, p::Pxest{4,T,Ptr{p4est_t}}) where {T}
    return p.pointer
end
function Base.unsafe_convert(::Type{Ptr{p8est_t}}, p::Pxest{8,T,Ptr{p8est_t}}) where {T}
    return p.pointer
end

Base.size(p::Pxest) =
    (convert(Int, (GC.@preserve p PointerWrapper(p.pointer).trees.elem_count[])),)
function Base.getindex(p::Pxest{X,T}, i::Int) where {X,T}
    @boundscheck checkbounds(p, i)
    GC.@preserve p begin
        TR = pxest_tree_t(Val(X))
        tree =
            Ptr{TR}(pointer(PointerWrapper(p.pointer).trees.array) + sizeof(TR) * (i - 1))
        return Tree{X,T,Ptr{TR},typeof(p)}(tree, p)
    end
end
Base.IndexStyle(::Pxest) = IndexLinear()

function iterate_volume_callback(info, _)
    info = unsafe_load(info)
    data =
        unsafe_pointer_to_objref(unsafe_load(unsafe_get_user_pointer_pointer(info.p4est)))[]
    X = quadrantstyle(data.forest)
    T = typeofquadrantuserdata(data.forest)
    quadrant = Quadrant{X,T,Ptr{pxest_quadrant_t(Val(X))}}(info.quad)
    data.volume(
        data.forest,
        data.ghost,
        quadrant,
        info.quadid + 1,
        info.treeid + 1,
        data.userdata,
    )
    return
end

@generated function generate_volume_callback(::Val{X}) where {X}
    I = pxest_iter_volume_info_t(Val(X))
    quote
        @cfunction(iterate_volume_callback, Cvoid, (Ptr{$I}, Ptr{Cvoid}))
    end
end

function iterateforest(
    forest::Pxest{X};
    ghost = nothing,
    volume = nothing,
    face = nothing,
    edge = nothing,
    corner = nothing,
    userdata = nothing,
) where {X}
    data = Ref((; forest, ghost, volume, face, edge, corner, userdata))

    ghost = isnothing(ghost) ? C_NULL : ghost
    volume::Ptr{Cvoid} = isnothing(volume) ? C_NULL : generate_volume_callback(Val(X))
    @assert face === nothing
    @assert edge === nothing
    @assert corner === nothing

    GC.@preserve data begin
        user_pointer_pointer = unsafe_get_user_pointer_pointer(forest)
        unsafe_store!(user_pointer_pointer, pointer_from_objref(data))

        if X == 4
            p4est_iterate(forest, ghost, C_NULL, volume, C_NULL, C_NULL)
        elseif X == 8
            p8est_iterate(forest, ghost, C_NULL, volume, C_NULL, C_NULL, C_NULL)
        else
            error("Not implemented")
        end

        unsafe_store!(user_pointer_pointer, C_NULL)
    end

    return
end

function init_callback(forest, treeid, quadrant)
    data = unsafe_pointer_to_objref(unsafe_load(unsafe_get_user_pointer_pointer(forest)))[]

    X = quadrantstyle(data.forest)
    T = typeofquadrantuserdata(data.forest)
    Q = pxest_quadrant_t(Val(X))

    quadrant = Quadrant{X,T,Ptr{Q}}(quadrant)

    data.init(data.forest, treeid + 1, quadrant)

    return
end

@generated function generate_init_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(init_callback, Cvoid, (Ptr{$P}, p4est_topidx_t, Ptr{$Q}))
    end
end

function replace_callback(forest, treeid, num_outgoing, outgoing, num_incoming, incoming)
    data = unsafe_pointer_to_objref(unsafe_load(unsafe_get_user_pointer_pointer(forest)))[]

    X = quadrantstyle(data.forest)
    T = typeofquadrantuserdata(data.forest)
    Q = pxest_quadrant_t(Val(X))

    outgoing = unsafe_wrap(Array, outgoing, num_outgoing)
    outgoing = ntuple(i -> Quadrant{X,T,Ptr{Q}}(outgoing[i]), num_outgoing)

    incoming = unsafe_wrap(Array, incoming, num_incoming)
    incoming = ntuple(i -> Quadrant{X,T,Ptr{Q}}(incoming[i]), num_incoming)

    data.replace(data.forest, treeid + 1, outgoing, incoming)

    return
end

@generated function generate_replace_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(
            replace_callback,
            Cvoid,
            (Ptr{$P}, p4est_topidx_t, Cint, Ptr{Ptr{$Q}}, Cint, Ptr{Ptr{$Q}})
        )
    end
end

function coarsen_callback(forest, treeid, children)
    data = unsafe_pointer_to_objref(unsafe_load(unsafe_get_user_pointer_pointer(forest)))[]

    X = quadrantstyle(data.forest)
    T = typeofquadrantuserdata(data.forest)
    Q = pxest_quadrant_t(Val(X))

    children = unsafe_wrap(Array, children, X)
    children = ntuple(i -> Quadrant{X,T,Ptr{Q}}(children[i]), Val(X))
    return data.coarsen(data.forest, treeid + 1, children) ? one(Cint) : zero(Cint)
end

@generated function generate_coarsen_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(coarsen_callback, Cint, (Ptr{$P}, p4est_topidx_t, Ptr{Ptr{$Q}}))
    end
end

function coarsen!(
    forest::Pxest{X};
    recursive = false,
    coarsen = (_...) -> false,
    init = nothing,
    replace = nothing,
) where {X}
    data = Ref((; forest, coarsen, init, replace))

    # Right now we do not support loading a orphans
    callback_orphans = false

    coarsen::Ptr{Cvoid} = isnothing(coarsen) ? C_NULL : generate_coarsen_callback(Val(X))
    init::Ptr{Cvoid} = isnothing(init) ? C_NULL : generate_init_callback(Val(X))
    replace::Ptr{Cvoid} = isnothing(replace) ? C_NULL : generate_replace_callback(Val(X))

    GC.@preserve data begin
        user_pointer_pointer = unsafe_get_user_pointer_pointer(forest)
        unsafe_store!(user_pointer_pointer, pointer_from_objref(data))

        (pxest_coarsen_ext(Val(X)))(
            forest,
            recursive,
            callback_orphans,
            coarsen,
            init,
            replace,
        )

        unsafe_store!(user_pointer_pointer, C_NULL)
    end

    return
end

function refine_callback(forest, treeid, quadrant)
    data = unsafe_pointer_to_objref(unsafe_load(unsafe_get_user_pointer_pointer(forest)))[]

    X = quadrantstyle(data.forest)
    T = typeofquadrantuserdata(data.forest)
    Q = pxest_quadrant_t(Val(X))

    quadrant = Quadrant{X,T,Ptr{Q}}(quadrant)
    return data.refine(data.forest, treeid + 1, quadrant) ? one(Cint) : zero(Cint)
end

@generated function generate_refine_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(refine_callback, Cint, (Ptr{$P}, p4est_topidx_t, Ptr{$Q}))
    end
end

function refine!(
    forest::Pxest{X};
    recursive = false,
    maxlevel = -1,
    refine = (_...) -> false,
    init = nothing,
    replace = nothing,
) where {X}
    data = Ref((; forest, refine, init, replace))

    refine::Ptr{Cvoid} = isnothing(refine) ? C_NULL : generate_refine_callback(Val(X))
    init::Ptr{Cvoid} = isnothing(init) ? C_NULL : generate_init_callback(Val(X))
    replace::Ptr{Cvoid} = isnothing(replace) ? C_NULL : generate_replace_callback(Val(X))

    GC.@preserve data begin
        user_pointer_pointer = unsafe_get_user_pointer_pointer(forest)
        unsafe_store!(user_pointer_pointer, pointer_from_objref(data))

        (pxest_refine_ext(Val(X)))(forest, recursive, maxlevel, refine, init, replace)

        unsafe_store!(user_pointer_pointer, C_NULL)
    end

    return
end

function balance!(
    forest::Pxest{X};
    connect = CONNECT_FULL(Val(X)),
    init = nothing,
    replace = nothing,
) where {X}
    data = Ref((; forest, init, replace))

    init::Ptr{Cvoid} = isnothing(init) ? C_NULL : generate_init_callback(Val(X))
    replace::Ptr{Cvoid} = isnothing(replace) ? C_NULL : generate_replace_callback(Val(X))

    GC.@preserve data begin
        user_pointer_pointer = unsafe_get_user_pointer_pointer(forest)
        unsafe_store!(user_pointer_pointer, pointer_from_objref(data))

        (pxest_balance_ext(Val(X)))(forest, connect, init, replace)

        unsafe_store!(user_pointer_pointer, C_NULL)
    end

    return
end

function weight_callback(forest, treeid, quadrant)
    data = unsafe_pointer_to_objref(unsafe_load(unsafe_get_user_pointer_pointer(forest)))[]

    X = quadrantstyle(data.forest)
    T = typeofquadrantuserdata(data.forest)
    Q = pxest_quadrant_t(Val(X))

    quadrant = Quadrant{X,T,Ptr{Q}}(quadrant)
    return data.weight(data.forest, treeid + 1, quadrant)
end

@generated function generate_weight_callback(::Val{X}) where {X}
    P = pxest_t(Val(X))
    Q = pxest_quadrant_t(Val(X))

    quote
        @cfunction(weight_callback, Cint, (Ptr{$P}, p4est_topidx_t, Ptr{$Q}))
    end
end

function partition!(
    forest::Pxest{X};
    ghost = nothing,
    lnodes_degree = nothing,
    allow_for_coarsening = false,
    weight = nothing,
) where {X}
    if !isnothing(lnodes_degree)
        if isnothing(ghost)
            ghost = ghostlayer(forest)
        end
        (pxest_partition_lnodes(Val(X)))(forest, ghost, lnodes_degree, allow_for_coarsening)
    else
        data = Ref((; forest, weight))

        weight::Ptr{Cvoid} = isnothing(weight) ? C_NULL : generate_weight_callback(Val(X))

        GC.@preserve data begin
            user_pointer_pointer = unsafe_get_user_pointer_pointer(forest)
            unsafe_store!(user_pointer_pointer, pointer_from_objref(data))

            (pxest_partition_ext(Val(X)))(forest, allow_for_coarsening, weight)

            unsafe_store!(user_pointer_pointer, C_NULL)
        end
    end

    return
end

import Base: show
function Base.show(io::IO, forest::P4estTypes.Pxest{X}) where {X}
    print(io, "Forest{$X} with $(length(forest)) trees.")
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, forest::P4estTypes.Pxest)
    print_tree(io, forest)
end

function Base.show(io::IO, tree::Tree{X}) where {X}
    print(io, "Tree{$X} with $(length(tree)) quadrants.")
end

function Base.show(io::IO, q::Quadrant{X}) where {X}
    print(io, "Quadrant{$X}: level $(level(q)), coordinates $(coordinates(q)).")
end
