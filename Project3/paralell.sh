#!/bin/bash
for n in 1 2 4
    do
    count=0
    while [ $count -lt 100 ]
        do 
        for mc in 100 1000 10000 100000 1000000
            do 
                mpirun -n $n python P3MPI.py $mc
            done
        count=$(( $count + 1 ))
        done
    done

