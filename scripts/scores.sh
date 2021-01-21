#!/bin/bash

proteins=(1crn  1gb1  1hhp  1i6c  1vii  T0868  T0900)

program="population_score"
mv ./build/$program $ROSETTA_BUILD_DIR/$program 2>/dev/null

for protein in ${proteins[@]};
do
    for i in `seq 1 32`;
    do
        $ROSETTA_BUILD_DIR/$program results/$protein/decoys/$i proteins/$protein/pdb 500 rmsd >> results/$protein/best_ca_rmsd2
        $ROSETTA_BUILD_DIR/$program results/$protein/decoys/$i proteins/$protein/pdb 500 gdt >> results/$protein/best_ca_gdt2
        #$ROSETTA_BUILD_DIR/$program results/$protein/decoys/$i proteins/$protein/pdb 500 mean_rmsd >> results/$protein/mean_ca_rmsd
        #$ROSETTA_BUILD_DIR/$program results/$protein/decoys/$i proteins/$protein/pdb 500 mean_gdt >> results/$protein/mean_ca_gdt
    done
done

