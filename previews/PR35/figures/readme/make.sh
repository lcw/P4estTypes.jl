#!/bin/sh

mpiexec -n 5 julia --project=../../../.. forest.jl
pvbatch forest_2d.py
pvbatch forest_3d.py
