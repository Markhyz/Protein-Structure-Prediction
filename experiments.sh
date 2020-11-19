#!/bin/bash

#OMP_PLACES=cores OMP_NUM_THREADS=10 numactl -N 1 -m 1 ./run.sh mobrkga_main 1hhp 1 ./ results/decoys_2/
#OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=4 ./run.sh mobrkga_main T0868 46 ./ results/decoys/

proteins=(1crn  1gb1  1hhp  1i6c  1vii  2imu  T0868  T0900)
offset=(1 1 1 1 1 1 46 5)

mv ./build/$1 $ROSETTA_BUILD_DIR/$1 2>/dev/null

mkdir -p results

for protein_idx in ${!proteins[@]};
do
    protein_name=${proteins[$protein_idx]}
    protein_offset=${offset[$protein_idx]}
    result_dir="results/$protein_name"

    mkdir -p $result_dir

    for i in `seq 1 1`;
    do
        cmd="./$ROSETTA_BUILD_DIR/$1"
        args="$protein_name $protein_offset $result_dir/ $result_dir/decoys/"
        args_list=""
        for j in `seq 1 2`;
        do
            mkdir -p $result_dir/decoys/$(($i + $j - 1))/
            args_list+="$args$(($i + $j - 1))/ "
        done
        printf '%s' "$args_list" | xargs -n 4 -P 2 $cmd
    done
done