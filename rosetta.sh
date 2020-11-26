#!/bin/bash

proteins=(1crn  1gb1  1hhp  1i6c  1vii  T0868  T0900)

program="$ROSETTA_BUILD_DIR/score_calc"
mv ./build/score_calc $ROSETTA_BUILD_DIR/score_calc 2>/dev/null

for protein in ${proteins[@]};
do
    for i in `seq 1 32`;
    do
        $program results/$protein/rosetta/decoy_$i.pdb proteins/$protein/pdb rmsd >> results/$protein/rosetta_rmsd
        $program results/$protein/rosetta/decoy_$i.pdb proteins/$protein/pdb gdt >> results/$protein/rosetta_gdt
    done
done

