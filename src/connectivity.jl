mutable struct Connectivity{X,P}
    pointer::P
    function Connectivity{4}(pointer::Ptr{P4est.LibP4est.p4est_connectivity})
        connectivity = new{4,typeof(pointer)}(pointer)
        finalizer(connectivity) do c
            p4est_connectivity_destroy(c.pointer)
            c.pointer = C_NULL
            return
        end
    end
    function Connectivity{8}(pointer::Ptr{P4est.LibP4est.p8est_connectivity})
        connectivity = new{8,typeof(pointer)}(pointer)
        finalizer(connectivity) do c
            p8est_connectivity_destroy(c.pointer)
            c.pointer = C_NULL
            return
        end
    end
end

function Connectivity{4}(name::Symbol, args...)
    if name == :unitsquare
        pointer = p4est_connectivity_new_unitsquare()
    elseif name == :periodic
        pointer = p4est_connectivity_new_periodic()
    elseif name == :rotwrap
        pointer = p4est_connectivity_new_rotwrap()
    elseif name == :twotrees
        pointer = p4est_connectivity_new_twotrees(args...)
    elseif name == :corner
        pointer = p4est_connectivity_new_corner()
    elseif name == :pillow
        pointer = p4est_connectivity_new_pillow()
    elseif name == :moebius
        pointer = p4est_connectivity_new_moebius()
    elseif name == :star
        pointer = p4est_connectivity_new_star()
    elseif name == :cubed
        pointer = p4est_connectivity_new_cubed()
    elseif name == :disk_nonperiodic
        pointer = p4est_connectivity_new_disk_nonperiodic()
    elseif name == :disk
        pointer = p4est_connectivity_new_disk(args...)
    elseif name == :icosahedron
        pointer = p4est_connectivity_new_icosahedron()
    elseif name == :shell2d
        pointer = p4est_connectivity_new_shell2d()
    elseif name == :disk2d
        pointer = p4est_connectivity_new_disk2d()
    else
        throw(ArgumentError("Unknown connectivity"))
    end
    return Connectivity{4}(pointer)
end

function Connectivity{8}(name::Symbol, args...)
    if name == :unitcube
        pointer = p8est_connectivity_new_unitcube()
    elseif name == :periodic
        pointer = p8est_connectivity_new_periodic()
    elseif name == :rotcubes
        pointer = p8est_connectivity_new_rotcubes()
    elseif name == :rotwrap
        pointer = p8est_connectivity_new_rotwrap()
    elseif name == :shell
        pointer = p8est_connectivity_new_shell()
    elseif name == :sphere
        pointer = p8est_connectivity_new_sphere()
    elseif name == :torus
        pointer = p8est_connectivity_new_torus(args...)
    elseif name == :twocubes
        pointer = p8est_connectivity_new_twocubes()
    elseif name == :twotrees
        pointer = p8est_connectivity_new_twotrees(args...)
    elseif name == :twowrap
        pointer = p8est_connectivity_new_twowrap()
    else
        throw(ArgumentError("Unknown connectivity"))
    end
    return Connectivity{8}(pointer)
end

function Base.unsafe_convert(
    ::Type{Ptr{p4est_connectivity}},
    c::Connectivity{4,Ptr{p4est_connectivity}},
)
    return c.pointer
end
function Base.unsafe_convert(
    ::Type{Ptr{p8est_connectivity}},
    c::Connectivity{8,Ptr{p8est_connectivity}},
)
    return c.pointer
end

function Base.:(==)(x::Connectivity{4}, y::Connectivity{4})
    p4est_connectivity_is_equivalent(x.pointer, y.pointer) == 1
end
function Base.:(==)(x::Connectivity{8}, y::Connectivity{8})
    p8est_connectivity_is_equivalent(x.pointer, y.pointer) == 1
end

Connectivity{4}(name::String) = Connectivity{4}(p4est_connectivity_read_inp(name))
Connectivity{8}(name::String) = Connectivity{8}(p8est_connectivity_read_inp(name))

Base.isvalid(c::Connectivity{4}) = p4est_connectivity_is_valid(c.pointer) == 1
Base.sizeof(c::Connectivity{4}) = Int(p4est_connectivity_memory_used(c.pointer))

Base.isvalid(c::Connectivity{8}) = p8est_connectivity_is_valid(c.pointer) == 1
Base.sizeof(c::Connectivity{8}) = Int(p8est_connectivity_memory_used(c.pointer))

function unsafe_vertices(c::Connectivity{X}) where {X}
    v = unsafe_wrap(
        Matrix{Cdouble},
        c.pointer.vertices,
        (3, Int(c.pointer.num_vertices)),
        own = false,
    )
    return reinterpret(reshape, NTuple{3,Cdouble}, v)
end

function unsafe_trees(c::Connectivity{X}) where {X}
    ttv = unsafe_wrap(
        Matrix{p4est_topidx_t},
        c.pointer.tree_to_vertex,
        (X, Int(c.pointer.num_trees)),
        own = false,
    )
    return reinterpret(reshape, NTuple{X,p4est_topidx_t}, ttv)
end

function unsafe_tree_to_tree(c::Connectivity{X}) where {X}
    sides = X == 4 ? 4 : 6
    ttt = unsafe_wrap(
        Matrix{p4est_topidx_t},
        c.pointer.tree_to_tree,
        (sides, Int(c.pointer.num_trees)),
        own = false,
    )
    return reinterpret(reshape, NTuple{sides,p4est_topidx_t}, ttt)
end

function unsafe_tree_to_face(c::Connectivity{X}) where {X}
    sides = X == 4 ? 4 : 6
    ttf = unsafe_wrap(
        Matrix{Int8},
        c.pointer.tree_to_face,
        (sides, Int(c.pointer.num_trees)),
        own = false,
    )
    return reinterpret(reshape, NTuple{sides,Int8}, ttf)
end

function unsafe_tree_to_corner(c::Connectivity{X}) where {X}
    s = c.pointer.tree_to_corner == C_NULL ? 0 : Int(c.pointer.num_trees)
    ttc = unsafe_wrap(Matrix{p4est_topidx_t}, c.pointer.tree_to_corner, (X, s), own = false)
    return reinterpret(reshape, NTuple{X,p4est_topidx_t}, ttc)
end

function unsafe_ctt_offset(c::Connectivity{X}) where {X}
    ctt = unsafe_wrap(
        Vector{p4est_topidx_t},
        c.pointer.ctt_offset,
        (Int(c.pointer.num_corners + 1),),
        own = false,
    )
    return ctt
end

function unsafe_corner_to_tree_array(c::Connectivity{X}) where {X}
    off = unsafe_ctt_offset(c)
    s = length(off) > 0 ? last(off) : 0
    return unsafe_wrap(Vector{p4est_topidx_t}, c.pointer.corner_to_tree, (s,), own = false)
end

function unsafe_corner_to_corner_array(c::Connectivity{X}) where {X}
    off = unsafe_ctt_offset(c)
    s = length(off) > 0 ? last(off) : 0
    return unsafe_wrap(Vector{Int8}, c.pointer.corner_to_corner, (s,), own = false)
end

function unsafe_corner_to_tree(c::Connectivity{X}) where {X}
    o = unsafe_ctt_offset(c)
    ctt = P4estTypes.unsafe_corner_to_tree_array(c)
    ctc = P4estTypes.unsafe_corner_to_corner_array(c)
    s = SparseMatrixCSC(c.pointer.num_trees, c.pointer.num_corners, o .+ 1, ctt .+ 1, ctc)
    return s
end

function refine(c::Connectivity{4}, nedge)
    return Connectivity{4}(p4est_connectivity_refine(c.pointer, nedge))
end
function refine(c::Connectivity{8}, nedge)
    return Connectivity{8}(p8est_connectivity_refine(c.pointer, nedge))
end

reduce!(c::Connectivity{4}) = p4est_connectivity_reduce(c.pointer)
reduce!(c::Connectivity{8}) = p8est_connectivity_reduce(c.pointer)

complete!(c::Connectivity{4}) = p4est_connectivity_complete(c.pointer)
complete!(c::Connectivity{8}) = p8est_connectivity_complete(c.pointer)

function brick(n::Tuple{Integer,Integer}, p::Tuple{Bool,Bool} = (false, false))
    return Connectivity{4}(p4est_connectivity_new_brick(n..., p...))
end
function brick(
    n::Tuple{Integer,Integer,Integer},
    p::Tuple{Bool,Bool,Bool} = (false, false, false),
)
    return Connectivity{8}(p8est_connectivity_new_brick(n..., p...))
end
brick(l::Integer, m::Integer, p::Bool = false, q::Bool = false) = brick((l, m), (p, q))
function brick(
    l::Integer,
    m::Integer,
    n::Integer,
    p::Bool = false,
    q::Bool = false,
    r::Bool = false,
)
    return brick((l, m, n), (p, q, r))
end

function Base.show(io::IO, c::Connectivity{X}) where {X}
    print(io, "Connectivity{", string(X), "}")
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, c::Connectivity{X}) where {X}
    show(io, c)
    if !get(io, :compact, false)
        GC.@preserve c begin
            # Should we print everthing 1 based?
            println(io, "\nnote: the following entries are zero-based")
            println(io, "\ntrees:")
            show(io, mime, unsafe_trees(c))
            println(io, "\nvertices:")
            show(io, mime, unsafe_vertices(c))
            println(io, "\ntree to tree:")
            show(io, mime, unsafe_tree_to_tree(c))
            println(io, "\ntree to face:")
            show(io, mime, unsafe_tree_to_face(c))
            println(io, "\ntree to corner:")
            show(io, mime, unsafe_tree_to_corner(c))
            println(io, "\ncorners:")
            show(io, mime, unsafe_corner_to_tree(c))
        end
    end
end

@recipe function f(c::Connectivity{4})
    verts = unsafe_vertices(c)

    xlims = extrema(getindex.(verts, 1))
    ylims = extrema(getindex.(verts, 2))
    zlims = extrema(getindex.(verts, 3))

    isconstz = zlims[1] == zlims[2]

    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"

    aspect_ratio --> :equal
    legend --> false
    grid --> false

    @series begin
        seriestype --> :path
        linecolor --> :gray
        linewidth --> 1

        x = []
        y = []
        z = []
        for tree in unsafe_trees(c)
            for i in (1, 2, 4, 3, 1)
                xi, yi, zi = verts[tree[i]+1]

                push!(x, xi)
                push!(y, yi)
                push!(z, zi)
            end

            push!(x, NaN)
            push!(y, NaN)
            push!(z, NaN)
        end
        if isconstz
            x, y
        else
            x, y, z
        end
    end
end

@recipe function f(c::Connectivity{8})
    verts = unsafe_vertices(c)

    xlims = extrema(getindex.(verts, 1))
    ylims = extrema(getindex.(verts, 2))
    zlims = extrema(getindex.(verts, 3))

    xlabel --> "x"
    ylabel --> "y"
    zlabel --> "z"

    aspect_ratio --> :equal
    legend --> false
    grid --> false

    @series begin
        seriestype --> :path
        linecolor --> :gray
        linewidth --> 1

        x = []
        y = []
        z = []
        for tree in unsafe_trees(c)
            for j in (0, 4)
                for i in (1 + j, 2 + j, 4 + j, 3 + j, 1 + j)
                    xi, yi, zi = verts[tree[i]+1]

                    push!(x, xi)
                    push!(y, yi)
                    push!(z, zi)
                end

                push!(x, NaN)
                push!(y, NaN)
                push!(z, NaN)
            end

            for j = 0:3
                for i in (1 + j, 5 + j)
                    xi, yi, zi = verts[tree[i]+1]

                    push!(x, xi)
                    push!(y, yi)
                    push!(z, zi)
                end

                push!(x, NaN)
                push!(y, NaN)
                push!(z, NaN)
            end
        end

        x, y, z
    end
end
