"""
    GhostLayer{X,T,P}

Stores a ghost layer of quadrants that neighbor the domain local to the rank
for a `Pxest{X,T}`.  Also stores the corresponding local domain quadrants,
mirrors, that are in other rank's ghost layers.

# Fields
$(DocStringExtensions.FIELDS)

# See also
- [`ghostlayer`](@ref): a function used to construct a `GhostLayer`
"""
mutable struct GhostLayer{X,T,P}
    """The pointer (of type `P`) can be a pointer to either a
    `P4estTypes.P4est.p4est_ghost_t` or a
    `P4estTypes.P4est.p8est_ghost_t`.  See the help
    documentation for these types for more information about the
    underlying p4est structures. """
    pointer::P
    function GhostLayer{4}(pointer::Ptr{p4est_ghost_t}, ::Type{T}) where {T}
        ghost = new{4,T,typeof(pointer)}(pointer)
        finalizer(ghost) do p
            p4est_ghost_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
    function GhostLayer{8}(pointer::Ptr{p8est_ghost_t}, ::Type{T}) where {T}
        ghost = new{8,T,typeof(pointer)}(pointer)
        finalizer(ghost) do p
            p8est_ghost_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
end

"""
    ghostlayer(forest::Pxest{X}; connection=P4estTypes.CONNECT_FULL(Val(X)))

Construct a ghost layer of quadrants that neighbor the local to the rank
for the given `forest`.  Here `connection` determines what neighboring
quadrants to include (across face, edge, corner, or full) and can take the
values:
- `P4estTypes.CONNECT_FULL(Val(4))`: get face and corner neighbors.
- `P4estTypes.CONNECT_FULL(Val(8))`: get face, edge, and corner neighbors.
- `P4estTypes.CONNECT_FACE(Val(4))`: get face neighbors.
- `P4estTypes.CONNECT_FACE(Val(8))`: get face neighbors.
- `P4estTypes.CONNECT_EDGE(Val(8))`: get face and edge neighbors.
- `P4estTypes.CONNECT_CORNER(Val(4))`: get face and corner neighbors.
- `P4estTypes.CONNECT_CORNER(Val(8)): `get face, edge, and corner neighbors.
"""
function ghostlayer(forest::Pxest{X,T}; connection = CONNECT_FULL(Val(X))) where {X,T}
    return GhostLayer{X}((pxest_ghost_new(Val(X)))(forest, connection), T)
end

struct Ghosts{X,T,P,Q} <: AbstractArray{Quadrant,1}
    pointer::P
    ghostlayer::GhostLayer{X,T,Q}
end

Base.size(g::Ghosts) = (PointerWrapper(g.pointer).elem_count[],)
Base.@propagate_inbounds function Base.getindex(g::Ghosts{X,T}, i::Int) where {X,T}
    @boundscheck checkbounds(g, i)
    GC.@preserve g begin
        Q = pxest_quadrant_t(Val(X))
        quadrant = Ptr{Q}(pointer(PointerWrapper(g.pointer).array) + sizeof(Q) * (i - 1))
        return Quadrant{X,T,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Ghosts) = IndexLinear()

struct Mirrors{X,T,P,Q} <: AbstractArray{Quadrant,1}
    pointer::P
    ghostlayer::GhostLayer{X,T,Q}
end

Base.size(m::Mirrors) = (PointerWrapper(m.pointer).elem_count[],)
Base.@propagate_inbounds function Base.getindex(m::Mirrors{X,T}, i::Int) where {X,T}
    @boundscheck checkbounds(m, i)
    GC.@preserve m begin
        Q = pxest_quadrant_t(Val(X))
        quadrant = Ptr{Q}(pointer(PointerWrapper(m.pointer).array) + sizeof(Q) * (i - 1))
        return Quadrant{X,T,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Mirrors) = IndexLinear()

"""
    ghosts(gl::GhostLayer)

Returns an array-like structure with the [`Quadrant`](@ref)s that neighbor the
domain of the local rank.
"""
function ghosts(gl::GhostLayer{X,T}) where {X,T}
    gs = pointer(PointerWrapper(gl.pointer).ghosts)
    return Ghosts{X,T,typeof(gs),typeof(gl.pointer)}(gs, gl)
end

"""
    mirrors(gl::GhostLayer)

Returns an array-like structure with the [`Quadrant`](@ref)s in the local
domain that are in neighboring rank's ghost layers.
"""
function mirrors(gl::GhostLayer{X,T}) where {X,T}
    ms = pointer(PointerWrapper(gl.pointer).mirrors)
    return Mirrors{X,T,typeof(ms),typeof(gl.pointer)}(ms, gl)
end

function Base.unsafe_convert(
    ::Type{Ptr{p4est_ghost_t}},
    p::GhostLayer{4,T,Ptr{p4est_ghost_t}},
) where {T}
    return p.pointer
end

function Base.unsafe_convert(
    ::Type{Ptr{p8est_ghost_t}},
    p::GhostLayer{8,T,Ptr{p8est_ghost_t}},
) where {T}
    return p.pointer
end
