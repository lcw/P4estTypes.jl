# Usage

P4estTypes.jl is a high-level Julia interface to [p4est](https://p4est.org) a
distributed (via [MPI](http://www.mpi-forum.org/)) forest-of-octrees C library.
P4estTypes.jl uses the lower-level wrappers from
[P4est.jl](https://github.com/trixi-framework/P4est.jl) to interface with
p4est.

## Installation

It is required to install [MPI.jl](http://github.com/JuliaParallel/MPI.jl)
alongside P4estTypes.jl.

```julia
julia> using Pkg; Pkg.add(["MPI", "P4estTypes"])
```

This uses the default MPI library distributed with MPI.jl.  If you would like
to use another MPI installation you need to make sure that MPI.jl points to the
correct MPI installation.  P4est.jl must also point to a p4est C library
compiled for the desired MPI installation.  Detailed instructions for doing
this may be found
[here](https://github.com/trixi-framework/P4est.jl#installation).

## Basic example

Here we present a basic example of a forest-of-quadtrees distributed among 2 MPI
ranks.  We will use [tmpi](https://github.com/Azrael3000/tmpi) to interactively
investigate the forest.  In practice MPI-based codes are launched in a
non-interactive mode via `mpiexec`.

We create a project with MPI.jl and P4estTypes.jl and ensure everything is
precompiled (we can run into trouble having multiple ranks trying to precompile
at the same time).

```sh
julia --project=. -e 'using Pkg; Pkg.add(["MPI", "P4estTypes"]); Pkg.API.precompile()'
```

We begin by launch an interactive session with 2 MPI ranks from our shell
prompt.
```sh
$ tmpi 2 julia --project
```

Next, we on both ranks we declare the packages we are using and initialize MPI.
Importantly, MPI must be initialized before P4estTypes functions are called.
```
julia> using P4estTypes, MPI; MPI.Init();
```

On each rank, a three element (in the x-direction) nonperiodic connectivity is
constructed with
```
julia> conn = brick(3, 1)
Connectivity{4}
note: the following entries are zero-based

trees:
3-element reinterpret(reshape, NTuple{4, Int32}, ::Matrix{Int32}) with eltype NTuple{4, Int32}:
 (0, 1, 2, 3)
 (1, 4, 3, 5)
 (4, 6, 5, 7)
vertices:
8-element reinterpret(reshape, Tuple{Float64, Float64, Float64}, ::Matrix{Float64}) with eltype Tuple{Float64, Float64, Float64}:
 (0.0, 0.0, 0.0)
 (1.0, 0.0, 0.0)
 (0.0, 1.0, 0.0)
 (1.0, 1.0, 0.0)
 (2.0, 0.0, 0.0)
 (2.0, 1.0, 0.0)
 (3.0, 0.0, 0.0)
 (3.0, 1.0, 0.0)
tree to tree:
3-element reinterpret(reshape, NTuple{4, Int32}, ::Matrix{Int32}) with eltype NTuple{4, Int32}:
 (0, 1, 0, 0)
 (0, 2, 1, 1)
 (1, 2, 2, 2)
tree to face:
3-element reinterpret(reshape, NTuple{4, Int8}, ::Matrix{Int8}) with eltype NTuple{4, Int8}:
 (0, 0, 2, 3)
 (1, 0, 2, 3)
 (1, 1, 2, 3)
tree to corner:
0-element reinterpret(reshape, NTuple{4, Int32}, ::Matrix{Int32}) with eltype NTuple{4, Int32}
corners:
3×0 SparseArrays.SparseMatrixCSC{Int8, Int64} with 0 stored entries
```
The REPL shows the underlying data stored for this connectivity (note that it
is zero-based).

Next we can construct a uniform forest with 1 level of refinement from this
connectivity.  We now show the REPL output of each rank individually as the
output is different.

!!! output "Rank 0 output building the forest"
    ```
    julia> forest = pxest(conn; min_level=1)
    Forest{4} with 3 trees.
    ├─ Tree{4} with 4 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ├─ Tree{4} with 2 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 0).
    └─ Tree{4} with 0 quadrants.
    ```

!!! output "Rank 1 output building the forest"
    ```
    julia> forest = pxest(conn; min_level=1)
    Forest{4} with 3 trees.
    ├─ Tree{4} with 0 quadrants.
    ├─ Tree{4} with 2 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    └─ Tree{4} with 4 quadrants.
       ├─ Quadrant{4}: level 1, coordinates (0, 0).
       ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
       ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
       └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ```

The forest is an array-of-arrays data structure containing all the quadrants.
Note that each rank stores only the quadrants in its rank-local domain.

Above, we see that rank 0 has an array for each tree in `conn` but only some of
those arrays have `Quadrants` in them.  Specifically, all of the quadrants from
tree 1 and half of the quadrants from tree 2.  Rank 1 contains the remaining
quadrants.  The coordinates of the quadrants are given in tree-local integer
coordinates.

We can get a visual representation of the forest by saving VTK files and viewing
them in [Paraview](https://paraview.org).  Running
```
julia> savevtk("initialforest", forest)
```
generates one `.vtu` file for each rank, here `initialforest_0000.vtu` and
`initialforest_0001.vtu`, and the collections file `initialforest.pvtu` and
`initialforest.visit`.

We currently have the sibling quadrants (quadrants that have the same parent)
in tree 2 distributed among both ranks.  We are not able to coarsen siblings
that live on different ranks.  So we may want to partition the mesh in a way
that keeps siblings together.  This can be done with
```
julia> partition!(forest; allow_for_coarsening=true)
```

We now show the mesh on each rank.

!!! output "Rank 0 output of the newly partitioned mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 4 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ├─ Tree{4} with 4 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    └─ Tree{4} with 0 quadrants.
    ```

!!! output "Rank 1 output of the newly partitioned mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 0 quadrants.
    ├─ Tree{4} with 0 quadrants.
    └─ Tree{4} with 4 quadrants.
       ├─ Quadrant{4}: level 1, coordinates (0, 0).
       ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
       ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
       └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ```

Rank 0 gets the first 2 sibling groups and rank 1 gets the other.

We can now coarsen the elements in tree 2 with
```
julia> coarsen!(forest; coarsen = (_, treeid, _) -> treeid == 2 )
```
which gives the distribution

!!! output "Rank 0 output of the coarsened mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 4 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ├─ Tree{4} with 1 quadrants.
    │  └─ Quadrant{4}: level 0, coordinates (0, 0).
    └─ Tree{4} with 0 quadrants.
    ```

!!! output "Rank 1 output of the coarsened mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 0 quadrants.
    ├─ Tree{4} with 0 quadrants.
    └─ Tree{4} with 4 quadrants.
       ├─ Quadrant{4}: level 1, coordinates (0, 0).
       ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
       ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
       └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ```

We see that tree 2 only has 1 quadrant now at the root of the tree.

We are now going to refine the quadrant in tree 1 with coordinates
`(536870912, 0)` with
```
julia> refine!(forest;
               refine = (_, treeid, quad) -> treeid == 1 &&
                                             coordinates(quad) == (536870912, 0)
              )
```
which gives

!!! output "Rank 0 output of the refined mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 7 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (805306368, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (536870912, 268435456).
    │  ├─ Quadrant{4}: level 2, coordinates (805306368, 268435456).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ├─ Tree{4} with 1 quadrants.
    │  └─ Quadrant{4}: level 0, coordinates (0, 0).
    └─ Tree{4} with 0 quadrants.
    ```

!!! output "Rank 1 output of the refined mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 0 quadrants.
    ├─ Tree{4} with 0 quadrants.
    └─ Tree{4} with 4 quadrants.
       ├─ Quadrant{4}: level 1, coordinates (0, 0).
       ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
       ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
       └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ```

We see that we have level 2 elements now in tree 1.  These neighbor the level 0
element in tree 2.  The mesh is unbalanced because there are quadrants differing
by more than one level that share an edge.  We bring this back to a 1 level
difference by enforcing the 2-to-1 constraint with
```
julia> balance!(forest)
```
which gives

!!! output "Rank 0 output of the balanced mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 7 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (805306368, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (536870912, 268435456).
    │  ├─ Quadrant{4}: level 2, coordinates (805306368, 268435456).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ├─ Tree{4} with 4 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    └─ Tree{4} with 0 quadrants.
    ```

!!! output "Rank 1 output of the balanced mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 0 quadrants.
    ├─ Tree{4} with 0 quadrants.
    └─ Tree{4} with 4 quadrants.
       ├─ Quadrant{4}: level 1, coordinates (0, 0).
       ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
       ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
       └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ```

There are now 11 quadrants on rank 0 and 4 on rank 1.  We can more evenly
distribute the quadrants with
```
julia> partition!(forest; allow_for_coarsening=true)
```
which gives

!!! output "Rank 0 output of the partitioned mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 7 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (805306368, 0).
    │  ├─ Quadrant{4}: level 2, coordinates (536870912, 268435456).
    │  ├─ Quadrant{4}: level 2, coordinates (805306368, 268435456).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ├─ Tree{4} with 0 quadrants.
    └─ Tree{4} with 0 quadrants.
    ```

!!! output "Rank 1 output of the partitioned mesh"
    ```
    julia> forest
    Forest{4} with 3 trees.
    ├─ Tree{4} with 0 quadrants.
    ├─ Tree{4} with 4 quadrants.
    │  ├─ Quadrant{4}: level 1, coordinates (0, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
    │  ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
    │  └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    └─ Tree{4} with 4 quadrants.
       ├─ Quadrant{4}: level 1, coordinates (0, 0).
       ├─ Quadrant{4}: level 1, coordinates (536870912, 0).
       ├─ Quadrant{4}: level 1, coordinates (0, 536870912).
       └─ Quadrant{4}: level 1, coordinates (536870912, 536870912).
    ```

This is a more even partitioning with 7 quadrants on rank 0 and 8 quadrants on
rank 2.

Now that we have the forest adapted as desired we would like to access the
quadrants of the forest.  There are currently two API for doing this.  The
first is the arrays-of-arrays approach.  Let's get the 2nd quadrant in
tree 1.

!!! output "Rank 0 output accessing 2nd quadrant in tree 1"
    ```
    julia> forest[1][2]
    Quadrant{4}: level 2, coordinates (536870912, 0).
    ```

!!! output "Rank 1 output accessing 2nd quadrant in tree 1"
    ```
    julia> forest[1][2]
    ERROR: BoundsError: attempt to access 0-element P4estTypes.Tree{4, Nothing, Ptr{P4est.LibP4est.p4est_tree}, Pxest{4, Nothing, Ptr{P4est.LibP4est.p4est}, Connectivity{4, Ptr{P4est.LibP4est.p4est_connectivity}}}} at index [2]
    Stacktrace:
     [1] throw_boundserror(A::P4estTypes.Tree{4, Nothing, Ptr{P4est.LibP4est.p4est_tree}, Pxest{4, Nothing, Ptr{P4est.LibP4est.p4est}, Connectivity{4, Ptr{P4est.LibP4est.p4est_connectivity}}}}, I::Tuple{Int64})
       @ Base ./abstractarray.jl:744
     [2] checkbounds
       @ ./abstractarray.jl:709 [inlined]
     [3] getindex(t::P4estTypes.Tree{4, Nothing, Ptr{P4est.LibP4est.p4est_tree}, Pxest{4, Nothing, Ptr{P4est.LibP4est.p4est}, Connectivity{4, Ptr{P4est.LibP4est.p4est_connectivity}}}}, i::Int64)
       @ P4estTypes ~/research/code/julia/P4estTypes/src/pxest.jl:315
     [4] top-level scope
       @ REPL[35]:1
    ```
    We get an error message here because rank 1 doesn't have any elements
    in tree 1.

We can also access the quadrants through the p4est iterator function,
[`iterateforest`](@ref).  This function iterates over the volumes, faces,
edges, and corners of the forest via callback functions.  To show the
quadrants with tree coordinates `(0, 0)` we can run
```
julia> iterateforest(forest;
                     volume = (_, _, quad, _, treeid, _) ->
                         (coordinates(quad) == (0, 0) &&
                          (@show (treeid, quad))))
```
to print

!!! output "Rank 0 output calling `iterateforest`"
    ```
    (treeid, quad) = (1, Quadrant{4}: level 1, coordinates (0, 0).)
    ```

!!! output "Rank 1 output calling `iterateforest`"
    ```
    (treeid, quad) = (2, Quadrant{4}: level 1, coordinates (0, 0).)
    (treeid, quad) = (3, Quadrant{4}: level 1, coordinates (0, 0).)
    ```
