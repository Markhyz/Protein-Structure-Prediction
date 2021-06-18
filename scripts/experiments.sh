#!/bin/bash

#OMP_PLACES=cores OMP_NUM_THREADS=10 numactl -N 1 -m 1 ./run.sh mobrkga_main 1hhp 1 ./ results/decoys_2/
#OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=4 ./run.sh mobrkga_main T0868 46 ./ results/decoys/

proteins=(1zdd 1enh 2mr9)
offset=(1 1 1)

mv ./build/$1 $ROSETTA_BUILD_DIR/$1 2>/dev/null

for protein_idx in ${!proteins[@]};
do
    protein_name=${proteins[$protein_idx]}
    protein_offset=${offset[$protein_idx]}

    for i in `seq 1 5`;
    do
        cmd="numactl"
        args_list=""
        for j in `seq 1 4`;
        do
            result_dir="results_convergence/$protein_name/$((($i - 1) * 4 + $j))"

            mkdir -p $result_dir/scores $result_dir/decoys $result_dir/algorithm $result_dir/graphs
        
            echo "
                FASTA_PATH = proteins/$protein_name/fasta
                PDB_PATH   = proteins/$protein_name/pdb
                FRAG3_PATH = proteins/$protein_name/frag3
                FRAG9_PATH = proteins/$protein_name/frag9
                SS_PATH    = proteins/$protein_name/ss2
                CM_PATH    = proteins/$protein_name/con      

                SCORE_OUTPUT_DIR      = $result_dir/scores
                DECOY_OUTPUT_DIR      = $result_dir/decoys
                ALGORITHM_OUTPUT_DIR  = $result_dir/algorithm

                ITERATIONS      = 1000
                POPULATION_SIZE = 500
            " > predictor_$j.conf

            args_list+="-N $((j - 1)) -m $((j - 1)) $ROSETTA_BUILD_DIR/$1 predictor_$j.conf "
        done
        printf '%s' "$args_list" | xargs -n 6 -P 4 $cmd
    done
done
