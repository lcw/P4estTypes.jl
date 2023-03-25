module P4estTypes

using AbstractTrees: print_tree
using MPI
using P4est
using RecipesBase
using SparseArrays

export Connectivity, Locidx, Gloidx

export refine, reduce!, complete!, brick, pxest
export lengthoflocalquadrants, lengthofglobalquadrants # Find better names
export level, storeuserdata!, loaduserdata
export offset
export iterateforest
export refine!, coarsen!, balance!, partition!
export lnodes, ghostlayer

include("sc.jl")

initialized() = P4est.package_id() >= 0
function initialize()
    if !SC.initialized()
        sc_init(MPI.COMM_NULL, 0, 0, C_NULL, SC_LP_ERROR)
    end
    if !initialized()
        p4est_init(C_NULL, SC_LP_ERROR)
    end
end

function setverbosity(logpriority::SC.LP.LogPriority)
    sc_package_set_verbosity(P4est.package_id(), logpriority)
end

uses_mpi() = P4est.uses_mpi()

include("connectivity.jl")
include("pxest.jl")
include("lnodes.jl")
include("ghost.jl")

function __init__()
    initialize()
end

using SnoopPrecompile
@precompile_setup begin
    @precompile_all_calls begin
        MPI.Initialized() || MPI.Init()
        initialize()
        p = pxest(brick(3, 4))
    end
end

end # module
