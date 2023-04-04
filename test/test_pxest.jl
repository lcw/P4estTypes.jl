using MPI
using MPIPreferences
using P4estTypes
using Test

if haskey(ENV, "P4ESTTYPES_TEST_BINARY")
    @test ENV["P4ESTTYPES_TEST_BINARY"] == MPIPreferences.binary
end

MPI.Init()

const comm = MPI.COMM_WORLD

let
    forest = pxest(brick(3, 4); comm)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 12
    @test MPI.Allreduce(lengthoflocalquadrants(forest), +, comm) == 12
    @test lengthofglobalquadrants(forest) == 12
    @test connectivity(forest) == brick(3, 4)
    @test size(forest) == (12,)
    @test P4estTypes.quadrantstyle(forest) == 4
    @test P4estTypes.quadrantndims(forest) == 2
end

let
    forest = pxest(brick(2, 4, 3); comm)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 24
    @test MPI.Allreduce(lengthoflocalquadrants(forest), +, comm) == 24
    @test lengthofglobalquadrants(forest) == 24
    @test connectivity(forest) == brick(2, 4, 3)
    @test size(forest) == (24,)
    @test P4estTypes.quadrantstyle(forest) == 8
    @test P4estTypes.quadrantndims(forest) == 3

    coarsen!(forest)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 24
    refine!(forest)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 24
    balance!(forest)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 24
    partition!(forest)
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

    struct Foo
        a::Int32
        b::Float64
        c::Float32
    end

    function foo_init(_, treeid, quadrant)
        data = Foo(treeid, rand(Float64), rand(Float32))
        storeuserdata!(quadrant, data)
        @test loaduserdata(quadrant) == data
    end

    function foo_check(_, _, quadrant, _, treeid, _)
        @test level(quadrant) == 2
        @test coordinates(quadrant) isa NTuple{2}
        @test loaduserdata(quadrant).a == treeid
    end

    forest = pxest(conn, min_level = 2, data_type = Foo, init_function = foo_init)
    iterateforest(forest; volume = foo_check)

    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 64
    @test MPI.Allreduce(lengthoflocalquadrants(forest), +, comm) == 64
    @test lengthofglobalquadrants(forest) == 64
    @test connectivity(forest) == conn
    @test size(forest) == (4,)
    @test P4estTypes.quadrantstyle(forest) == 4
    @test P4estTypes.quadrantndims(forest) == 2
    @test P4estTypes.typeofquadrantuserdata(forest) == Foo

    coarsen!(forest)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 64
    refine!(forest)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 64
    balance!(forest)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 64
    partition!(forest)
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 64

    function replace(_, _, outgoing, incoming)
        for q in incoming
            d = loaduserdata(q)
            storeuserdata!(q, Foo(d.a, length(outgoing), 5.0))
        end
    end
    partition!(forest; allow_for_coarsening = true)
    coarsen!(forest; init = foo_init, replace, coarsen = (_, t, _) -> iseven(t))
    iterateforest(forest; volume = (_, _, q, _, t, _) -> (@test loaduserdata(q).a == t))
    iterateforest(
        forest;
        volume = (_, _, q, _, t, _) -> (@test isodd(t) || loaduserdata(q).b == 4.0),
    )
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 40
    refine!(forest; init = foo_init, replace, refine = (_, t, _) -> iseven(t))
    iterateforest(forest; volume = foo_check)
    iterateforest(
        forest;
        volume = (_, _, q, _, t, _) -> (@test isodd(t) || loaduserdata(q).b == 1.0),
    )
    @test MPI.Allreduce(sum(length.(forest)), +, comm) == 64

    refine!(
        forest;
        init = foo_init,
        replace,
        refine = (_, t, q) -> t == 1 && coordinates(q) == (0, 0),
    )
    balance!(forest; init = foo_init, replace)
    iterateforest(forest; volume = (_, _, q, _, t, _) -> (@test loaduserdata(q).a == t))
    partition!(forest)
    iterateforest(forest; volume = (_, _, q, _, t, _) -> (@test loaduserdata(q).a == t))
end

MPI.Finalize()
