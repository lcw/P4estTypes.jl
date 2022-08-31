@inline pxest_t(::Val{4}) = p4est
@inline pxest_t(::Val{8}) = p8est
@inline pxest_new_ext(::Val{4}) = p4est_new_ext
@inline pxest_new_ext(::Val{8}) = p8est_new_ext
@inline pxest_quadrant_t(::Val{4}) = p4est_quadrant
@inline pxest_quadrant_t(::Val{8}) = p8est_quadrant
@inline pxest_tree_t(::Val{4}) = p4est_tree
@inline pxest_tree_t(::Val{8}) = p8est_tree
@inline pxest_iter_volume_info_t(::Val{4}) = p4est_iter_volume_info_t
@inline pxest_iter_volume_info_t(::Val{8}) = p8est_iter_volume_info_t

const Locidx = P4est.p4est_locidx_t
const Gloidx = P4est.p4est_gloidx_t

struct Quadrant{X,T,P}
    pointer::P
end

level(quadrant::Quadrant) = quadrant.pointer.level
storeuserdata!(quadrant::Quadrant{X,T}, data::T) where {X,T} =
    unsafe_store!(Ptr{T}(quadrant.pointer.p.user_data), data)
loaduserdata(quadrant::Quadrant{X,T}) where {X,T} =
    unsafe_load(Ptr{T}(quadrant.pointer.p.user_data))

struct Tree{X,T,P} <: AbstractArray{Quadrant,1}
    pointer::P
end

Base.size(t::Tree) = (t.pointer.quadrants.elem_count,)
function Base.getindex(t::Tree{X,T}, i::Int) where {X,T}
    @boundscheck checkbounds(t, i)
    GC.@preserve t begin
        Q = pxest_quadrant_t(Val(X))
        quadrant = Ptr{Q}(t.pointer.quadrants.array + sizeof(Q) * (i - 1))
        return Quadrant{X,T,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Tree) = IndexLinear()

offset(tree::Tree) = tree.pointer.quadrants_offset

mutable struct Pxest{X,T,P,C} <: AbstractArray{Tree,1}
    pointer::P
    connectivity::C
    comm::MPI.Comm
    function Pxest{4}(
        pointer::Ptr{p4est},
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
        pointer::Ptr{p8est},
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
        init(forest, _, quadrant, _, treeid) = init_function(forest, treeid, quadrant)
        iterateforest(forest; volume = init)
    end

    return forest
end

quadrantstyle(::Pxest{X}) where {X} = X
quadrantndims(::Pxest{4}) = 2
quadrantndims(::Pxest{8}) = 3
typeofquadrantuserdata(::Pxest{X,T}) where {X,T} = T
lengthoflocalquadrants(p::Pxest) = p.pointer.local_num_quadrants
lengthofglobalquadrants(p::Pxest) = p.pointer.global_num_quadrants

function Base.unsafe_convert(::Type{Ptr{p4est}}, p::Pxest{4,T,Ptr{p4est}}) where {T}
    return p.pointer
end
function Base.unsafe_convert(::Type{Ptr{p8est}}, p::Pxest{8,T,Ptr{p8est}}) where {T}
    return p.pointer
end

Base.size(p::Pxest) = (p.pointer.trees.elem_count,)
function Base.getindex(p::Pxest{X,T}, i::Int) where {X,T}
    @boundscheck checkbounds(p, i)
    GC.@preserve p begin
        tree = unsafe_load(Ptr{pxest_tree_t(Val(X))}(p.pointer.trees.array), i)
        return Tree{X,T,pxest_tree_t(Val(X))}(tree)
    end
end
Base.IndexStyle(::Pxest) = IndexLinear()

function iterate_volume_callback(info, _)
    data = unsafe_pointer_to_objref(info.p4est.user_pointer)[]
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
    forest.pointer.user_pointer = pointer_from_objref(data)

    volume::Ptr{Cvoid} = isnothing(volume) ? C_NULL : generate_volume_callback(Val(X))
    @assert face === nothing
    @assert edge === nothing
    @assert corner === nothing

    GC.@preserve data begin
        if X == 4
            p4est_iterate(forest, ghost, C_NULL, volume, C_NULL, C_NULL)
        elseif X == 8
            p8est_iterate(forest, ghost, C_NULL, volume, C_NULL, C_NULL, C_NULL)
        else
            error("Not implemented")
        end
    end

    forest.pointer.user_pointer = C_NULL

    return
end
