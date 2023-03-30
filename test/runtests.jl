using MPI
using MPIPreferences
using P4estTypes
using Test

nprocs_str = get(ENV, "P4ESTTYPES_TEST_NPROCS", "")
nprocs = nprocs_str == "" ? clamp(Sys.CPU_THREADS, 2, 4) : parse(Int, nprocs_str)

@info "Running MPI-based tests" nprocs MPIPreferences.abi MPIPreferences.binary

testdir = @__DIR__
istest(f) = endswith(f, ".jl") && startswith(f, "test_")
testfiles = sort(filter(istest, readdir(testdir)))

@testset "$f" for f in testfiles
    mpiexec() do mpirun
        run(
            `$mpirun -n $nprocs $(Base.julia_cmd()) --startup-file=no $(joinpath(testdir, f))`,
        )
        @test true
    end
end
