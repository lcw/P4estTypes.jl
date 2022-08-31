mutable struct GhostLayer{X,P}
    pointer::P
    function GhostLayer{4}(pointer::Ptr{p4est_ghost_t})
        ghost = new{4,typeof(pointer)}(pointer)
        finalizer(ghost) do p
            p4est_ghost_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
    function GhostLayer{8}(pointer::Ptr{p8est_ghost_t})
        ghost = new{8,typeof(pointer)}(pointer)
        finalizer(ghost) do p
            p8est_ghost_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
end

function ghostlayer(forest::Pxest{X}; connection = CONNECT_FULL(Val(X))) where {X}
    return GhostLayer{X}((pxest_ghost_new(Val(X)))(forest, connection))
end

function Base.unsafe_convert(
    ::Type{Ptr{p4est_ghost_t}},
    p::GhostLayer{4,Ptr{p4est_ghost_t}},
)
    return p.pointer
end

function Base.unsafe_convert(
    ::Type{Ptr{p8est_ghost_t}},
    p::GhostLayer{8,Ptr{p8est_ghost_t}},
)
    return p.pointer
end
