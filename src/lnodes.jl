mutable struct LNodes{X,P}
    pointer::P
    function LNodes{4}(pointer::Ptr{P4est.LibP4est.p4est_lnodes})
        nodes = new{4,typeof(pointer)}(pointer)
        finalizer(nodes) do p
            p4est_lnodes_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
    function LNodes{8}(pointer::Ptr{P4est.LibP4est.p8est_lnodes})
        nodes = new{8,typeof(pointer)}(pointer)
        finalizer(nodes) do p
            p8est_lnodes_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
end

function lnodes(forest::Pxest{X}; ghost = nothing, degree = 1) where {X}
    if isnothing(ghost)
        ghost = ghostlayer(forest)
    end
    return LNodes{X}((pxest_lnodes_new(Val(X)))(forest, ghost, degree))
end

function Base.unsafe_convert(::Type{Ptr{p4est_lnodes}}, p::LNodes{4,Ptr{p4est_lnodes}})
    return p.pointer
end

function Base.unsafe_convert(::Type{Ptr{p8est_lnodes}}, p::LNodes{8,Ptr{p8est_lnodes}})
    return p.pointer
end
