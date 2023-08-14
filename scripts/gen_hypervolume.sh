#!/bin/bash

for i in `seq 1 $3`;
do
    hv=`python3 scripts/fitness_convergence.py $1$i`
    echo $i $hv >> $2
done