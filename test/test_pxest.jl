using MPI
using MPIPreferences
using P4estTypes
using Test

MPI.Init()

const comm = MPI.COMM_WORLD

let
    forest = pxest(brick(3, 4); comm)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 12
end

let
    forest = pxest(brick(2, 4, 3); comm)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 24
end

let
    VXYZ = [
        -1.0 -1.0 0.0
        -1.0 0.0 0.0
        -1.0 1.0 0.0
        0.0 -1.0 0.0
        0.0 0.0 0.0
        0.0 1.0 0.0
        1.0 -1.0 0.0
        1.0 0.0 0.0
        1.0 1.0 0.0
    ]

    EToV = [
        0 3 1 4
        3 6 4 7
        1 4 2 5
        4 7 5 8
    ]

    conn = Connectivity{4}(VXYZ, EToV)
    @test isvalid(conn)

    forest = pxest(conn, min_level = 2)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 64
end

MPI.Finalize()
