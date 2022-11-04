@testset "Connectivity" begin
    @test isvalid(Connectivity{4}(:star))
    @test isvalid(Connectivity{8}(:twotrees, 1, 2, 3))
    @test isvalid(Connectivity{4}("data/hole_2d_gmsh.inp"))
    @test isvalid(brick((2, 3)))
    @test isvalid(brick((3, 2, 3)))
    @test brick((3, 4)) == brick((3, 4))
    @test brick((3, 4, 5), (false, false, true)) != brick((3, 4, 5))
    @test isvalid(refine(Connectivity{4}(:star), 2))
    @test isvalid(refine(Connectivity{8}(:twotrees, 1, 2, 3), 2))

    for b in (brick((4, 3)), brick((2, 3, 1)))
        reduce!(b)
        @test isvalid(b)
        complete!(b)
        @test isvalid(b)
    end

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

    # test default `refine` argument
    p4est = pxest(conn, min_level=2)
    @test length(p4est)==4 # 4 trees
    @test all(length.(p4est) .== 16) # each tree has 16 quadrants due to min_level = 2
end
