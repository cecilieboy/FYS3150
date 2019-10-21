#!/bin/bash
for n in 1 2 4
    do
    for mc in 100 1000 10000 100000 1000000 10000000
        do 
            mpirun -n $n python P3MPI.py 100 $mc
        done  
    done

