module P4estTypes

using CBinding
using MPI
using P4est
using RecipesBase
using SparseArrays

export Connectivity, Locidx, Gloidx

export refine, reduce!, complete!, brick, pxest, local_size, global_size, global_axes
export level, storeuserdata!, loaduserdata
export offset
export iterateforest

include("sc.jl")

initialized() = p4est_package_id()[] >= 0

function setverbosity(logpriority::SC.LP.LogPriority)
    sc_package_set_verbosity(p4est_package_id()[], logpriority)
end

uses_mpi() = P4est.uses_mpi()

include("connectivity.jl")
include("pxest.jl")

function __init__()
    if !SC.initialized()
        sc_init(MPI.COMM_NULL, 0, 0, C_NULL, SC_LP_ERROR)
    end
    if !initialized()
        p4est_init(C_NULL, SC_LP_ERROR)
    end
end

using SnoopPrecompile
@precompile_setup begin
    @precompile_all_calls begin
        MPI.Initialized() || MPI.Init()
        p = pxest(brick(3, 4))
    end
end

end # module
