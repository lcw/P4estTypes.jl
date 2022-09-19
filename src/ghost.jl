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
