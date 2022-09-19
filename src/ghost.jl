mutable struct GhostLayer{X,T,P}
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

function ghostlayer(forest::Pxest{X,T}; connection = CONNECT_FULL(Val(X))) where {X,T}
    return GhostLayer{X}((pxest_ghost_new(Val(X)))(forest, connection), T)
end

struct Ghosts{X,T,P} <: AbstractArray{Quadrant,1}
    pointer::P
end

Base.size(g::Ghosts) = (g.pointer.elem_count,)
Base.@propagate_inbounds function Base.getindex(g::Ghosts{X,T}, i::Int) where {X,T}
    @boundscheck checkbounds(g, i)
    GC.@preserve g begin
        Q = pxest_quadrant_t(Val(X))
        quadrant = Ptr{Q}(g.pointer.array + sizeof(Q) * (i - 1))
        return Quadrant{X,T,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Ghosts) = IndexLinear()

struct Mirrors{X,T,P} <: AbstractArray{Quadrant,1}
    pointer::P
end

Base.size(m::Mirrors) = (m.pointer.elem_count,)
Base.@propagate_inbounds function Base.getindex(m::Mirrors{X,T}, i::Int) where {X,T}
    @boundscheck checkbounds(m, i)
    GC.@preserve m begin
        Q = pxest_quadrant_t(Val(X))
        quadrant = Ptr{Q}(m.pointer.array + sizeof(Q) * (i - 1))
        return Quadrant{X,T,Ptr{Q}}(quadrant)
    end
end
Base.IndexStyle(::Mirrors) = IndexLinear()

ghosts(gl::GhostLayer{X,T}) where {X,T} =
    Ghosts{X,T,typeof(gl.pointer.ghosts)}(gl.pointer.ghosts)
mirrors(gl::GhostLayer{X,T}) where {X,T} =
    Mirrors{X,T,typeof(gl.pointer.mirrors)}(gl.pointer.mirrors)

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
