#!/bin/bash

#OMP_PLACES=cores OMP_NUM_THREADS=10 numactl -N 1 -m 1 ./run.sh mobrkga_main 1hhp 1 ./ results/decoys_2/
#OMP_PROC_BIND=close OMP_PLACES=cores OMP_NUM_THREADS=4 ./run.sh mobrkga_main T0868 46 ./ results/decoys/

proteins=(1zdd 1gb1 1ail 1hhp)
offset=(1 1 1 1)
threads=(1 2 4 6 8 10 20 40 80)

mv ./build/$1 $ROSETTA_BUILD_DIR/$1 2>/dev/null

mkdir -p results

for thr in ${threads[@]};
do
    for protein_idx in ${!proteins[@]};
    do
      protein_name=${proteins[$protein_idx]}
      protein_offset=${offset[$protein_idx]}
    
      for i in `seq 1 10`;
      do
        result_dir="results_parallel/$thr/$protein_name/$i"

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
            POPULATION_SIZE = 100
        " > predictor.conf

        OMP_NUM_THREADS=$thr $ROSETTA_BUILD_DIR/$1 predictor.conf
      done
    done
done
