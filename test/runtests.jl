using MPI
using P4estTypes
using Test

MPI.Initialized() || MPI.Init()

include("connectivity.jl")
