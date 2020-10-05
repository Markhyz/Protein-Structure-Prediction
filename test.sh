#!/bin/bash

for i in `seq 1 10`;
do
    OMP_NUM_THREADS=2 ./run.sh 4 1hhp 1 results/1hhp.rmsd
done