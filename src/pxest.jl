
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
        Q = X == 4 ? p4est_quadrant : p8est_quadrant
        quadurant = unsafe_load(Ptr{Q}(t.pointer.quadrants.array), i)
        return Quadrant{X,T,Q}(quadurant)
    end
end
Base.IndexStyle(::Tree) = IndexLinear()

offset(tree::Tree) = tree.pointer.quadrants_offset

mutable struct Pxest{X,T,P,C} <: AbstractArray{Tree,1}
    pointer::P
    connectivity::C
    comm::MPI.Comm
    function Pxest{4}(
        pointer::Ptr{P4est.LibP4est.p4est},
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
        pointer::Ptr{P4est.LibP4est.p8est},
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
    connectivity::Connectivity{4};
    comm = MPI.COMM_WORLD,
    min_quadrants = 0,
    min_level = 0,
    fill_uniform = true,
    data_type = Nothing,
    init_function = nothing,
)
    MPI.Initialized() || MPI.Init()

    pointer = p4est_new_ext(
        comm,
        connectivity,
        min_quadrants,
        min_level,
        fill_uniform,
        sizeof(data_type),
        C_NULL,
        C_NULL,
    )

    forest = Pxest{4}(pointer, connectivity, comm, data_type)

    if !isnothing(init_function)
        init(forest, _, quadrant, _, treeid) = init_function(forest, treeid, quadrant)
        iterateforest(forest; volume = init)
    end

    return forest
end

function pxest(
    connectivity::Connectivity{8};
    comm = MPI.COMM_WORLD,
    min_quadrants = 0,
    min_level = 0,
    fill_uniform = true,
    data_type = Nothing,
    init_function = nothing,
)
    MPI.Initialized() || MPI.Init()

    @assert isnothing(init_function)

    pointer = p8est_new_ext(
        comm,
        connectivity,
        min_quadrants,
        min_level,
        fill_uniform,
        sizeof(data_type),
        C_NULL,
        C_NULL,
    )
    return Pxest{8}(pointer, connectivity, comm, data_type)
end

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
        tree = unsafe_load(Ptr{p4est_tree}(p.pointer.trees.array), i)
        return Tree{X,T,p4est_tree}(tree)
    end
end
Base.IndexStyle(::Pxest) = IndexLinear()

function _p4est_volume_callback(info, _)
    data = unsafe_pointer_to_objref(info.p4est.user_pointer)[]
    T = typeofquadrantuserdata(data.forest)
    quadrant = Quadrant{4,T,Ptr{p4est_quadrant}}(info.quad)
    data.volume(data.forest, data.ghost, quadrant, info.quadid + 1, info.treeid + 1)
    return
end

function iterateforest(
    forest::Pxest{4};
    ghost = nothing,
    volume = nothing,
    face = nothing,
    corner = nothing,
)
    data = Ref((; forest, ghost, volume, face, corner))
    forest.pointer.user_pointer = pointer_from_objref(data)

    volume_callback::Ptr{Cvoid} =
        isnothing(volume) ? C_NULL :
        @cfunction(
            _p4est_volume_callback,
            Cvoid,
            (Ptr{p4est_iter_volume_info_t}, Ptr{Cvoid})
        )
    @assert face === nothing
    @assert corner === nothing

    GC.@preserve data begin
        p4est_iterate(forest, C_NULL, C_NULL, volume_callback, C_NULL, C_NULL)
    end

    forest.pointer.user_pointer = C_NULL

    return
end
