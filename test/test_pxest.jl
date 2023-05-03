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

    @test_nowarn P4estTypes.savevtk("basicbrick", forest)
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

    ghost = ghostlayer(forest)
    @test mirrors(ghost) isa AbstractArray{Quadrant,1}
    @test ghosts(ghost) isa AbstractArray{Quadrant,1}
    @test size(mirrors(ghost), 1) >= 0
    @test size(ghosts(ghost), 1) >= 0

    # Make sure we can index into the ghosts and mirrors array
    @test collect(mirrors(ghost)) isa AbstractArray{Quadrant,1}
    @test collect(ghosts(ghost)) isa AbstractArray{Quadrant,1}

    @test lnodes(forest; ghost, degree = 3) isa LNodes

    @test_nowarn P4estTypes.savevtk("basicconn", forest)
end

let
    # Create a comm with only 2 ranks
    worldgroup = MPI.Comm_group(MPI.COMM_WORLD)
    group = MPI.Group_incl(worldgroup, Cint[0, 1])
    twocomm = MPI.Comm_create(MPI.COMM_WORLD, group)

    if twocomm != MPI.COMM_NULL
        tworank = MPI.Comm_rank(twocomm)
        twosize = MPI.Comm_size(twocomm)
        @test twosize == 2

        forest = pxest(brick(2, 2); comm = twocomm)
        refine!(forest; refine = (_, tid, _) -> tid == 4)

        tworank == 0 && @test lengthoflocalquadrants(forest) == 2
        tworank == 1 && @test lengthoflocalquadrants(forest) == 5

        ghost = ghostlayer(forest)
        gs = ghosts(ghost)
        ms = mirrors(ghost)
        nodes = lnodes(forest; ghost, degree = 2)

        if tworank == 0
            @test lengthoflocalquadrants(forest) == 2
            @test P4estTypes.coordinates.(gs) ==
                  Tuple{Int32,Int32}[(0, 0), (0, 0), (536870912, 0)]
            @test P4estTypes.coordinates.(ms) == Tuple{Int32,Int32}[(0, 0), (0, 0)]
            @test P4estTypes.unsafe_which_tree.(gs) == Int32[3, 4, 4]
            @test P4estTypes.unsafe_which_tree.(ms) == Int32[1, 2]
            @test P4estTypes.unsafe_local_num.(gs) == Int32[1, 2, 3]
            @test P4estTypes.unsafe_local_num.(ms) == Int32[1, 2]

            GC.@preserve forest ghost nodes begin
                @test P4estTypes.unsafe_mirror_proc_offsets(ghost) == Int32[0, 0, 2]
                @test P4estTypes.unsafe_proc_offsets(ghost) == Int32[0, 0, 3]
                @test P4estTypes.unsafe_mirror_proc_mirrors(ghost) == Int32[0, 1]
                @test P4estTypes.unsafe_global_owned_count(nodes) == Int32[15, 22]
                @test P4estTypes.unsafe_face_code(nodes) == Int8[0, 0]
                @test P4estTypes.unsafe_global_first_quadrant(forest) == [0, 2, 7]
            end
        end

        if tworank == 1
            @test lengthoflocalquadrants(forest) == 5
            @test P4estTypes.coordinates.(gs) == Tuple{Int32,Int32}[(0, 0), (0, 0)]
            @test P4estTypes.coordinates.(ms) ==
                  Tuple{Int32,Int32}[(0, 0), (0, 0), (536870912, 0)]
            @test P4estTypes.unsafe_which_tree.(gs) == Int32[1, 2]
            @test P4estTypes.unsafe_which_tree.(ms) == Int32[3, 4, 4]
            @test P4estTypes.unsafe_local_num.(gs) == Int32[1, 2]
            @test P4estTypes.unsafe_local_num.(ms) == Int32[1, 2, 3]

            GC.@preserve forest ghost nodes begin
                @test P4estTypes.unsafe_mirror_proc_offsets(ghost) == Int32[0, 3, 3]
                @test P4estTypes.unsafe_proc_offsets(ghost) == Int32[0, 2, 2]
                @test P4estTypes.unsafe_mirror_proc_mirrors(ghost) == Int32[0, 1, 2]
                @test P4estTypes.unsafe_global_owned_count(nodes) == Int32[15, 22]
                @test P4estTypes.unsafe_face_code(nodes) == Int8[0, 12, 9, 6, 0]
                @test P4estTypes.unsafe_global_first_quadrant(forest) == [0, 2, 7]
            end
        end
    end
end

MPI.Finalize()
