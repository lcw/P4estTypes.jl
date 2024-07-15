using P4estTypes, MPI;
MPI.Init();

forest = pxest(brick(3, 2))
refine = (_, t, q) -> t == 1 && coordinates(q) == (0, 0)
refine!(forest; refine, recursive = true, maxlevel = 5)
balance!(forest)
partition!(forest)

P4estTypes.savevtk("forest_2d", forest)

forest = pxest(brick(3, 2, 4))
refine = (_, t, q) -> t == 1 && coordinates(q) == (0, 0, 0)
refine!(forest; refine, recursive = true, maxlevel = 5)
balance!(forest)
partition!(forest)

P4estTypes.savevtk("forest_3d", forest)
