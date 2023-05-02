"""
    LNodes{X,P}

Stores a parallel numbering of Lobatto nodes for a `Pxest{X}`.

# Fields
$(DocStringExtensions.FIELDS)

# See also
- [`lnodes`](@ref): a function used to construct `LNodes`
"""
mutable struct LNodes{X,P}
    """The pointer (of type `P`) can be a pointer to either a
    `P4estTypes.P4est.p4est_lnodes_t` or a
    `P4estTypes.P4est.p8est_lnodes_t`.  See the help
    documentation for these types for more information about the
    underlying p4est structures. """
    pointer::P
    """The MPI Communicator that includes the ranks participating in the
    lnods."""
    comm::MPI.Comm
    function LNodes{4}(pointer::Ptr{P4est.LibP4est.p4est_lnodes}, comm::MPI.Comm)
        nodes = new{4,typeof(pointer)}(pointer, comm)
        finalizer(nodes) do p
            p4est_lnodes_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
    function LNodes{8}(pointer::Ptr{P4est.LibP4est.p8est_lnodes}, comm::MPI.Comm)
        nodes = new{8,typeof(pointer)}(pointer, comm)
        finalizer(nodes) do p
            p8est_lnodes_destroy(p.pointer)
            p.pointer = C_NULL
            return
        end
    end
end

"""
    lnodes(forest::Pxest; ghost = nothing, degree = 1)

Construct a parallel node numbering for the `forest`. If the
ghost layer `ghost` is not provided it will be constructed.

A `degree > 0` indicates that degree `degree` Lobotto nodes will be
constructed.  A `degree < 0` indicates that the boundary objects
(faces, edges, and corners) will be numbered.

See `@doc P4estTypes.P4est.p4est_lnodes_t` and
`@doc P4estTypes.P4est.p8est_lnodes_t` for a more detailed discussion
of the numbering based on `degree`.
"""
function lnodes(forest::Pxest{X}; ghost = nothing, degree = 1) where {X}
    if isnothing(ghost)
        ghost = ghostlayer(forest)
    end
    return LNodes{X}((pxest_lnodes_new(Val(X)))(forest, ghost, degree), comm(forest))
end

function Base.unsafe_convert(::Type{Ptr{p4est_lnodes}}, p::LNodes{4,Ptr{p4est_lnodes}})
    return p.pointer
end

function Base.unsafe_convert(::Type{Ptr{p8est_lnodes}}, p::LNodes{8,Ptr{p8est_lnodes}})
    return p.pointer
end
