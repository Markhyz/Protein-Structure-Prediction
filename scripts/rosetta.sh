#!/bin/bash

#proteins=(1crn  1gb1  1hhp  1i6c  1vii  T0868  T0900)
#proteins=(1acw  1ail  1crn  1enh  1rop  1zdd  2mr9  2p81)

proteins=(1ab1 1acw 1ail 1aly 1bdd 1crn 1dfn 1enh 1gb1 1hhp 1i6c 1rop 1zdd 2kdl 2mr9 2p81 T0868 T0900 T0968s1 T1010)

program="$ROSETTA_BUILD_DIR/score_calc"
mv ./build/score_calc $ROSETTA_BUILD_DIR/score_calc 2>/dev/null

for protein in ${proteins[@]};
do
    for i in `seq 1 20`;
    do
	mkdir -p statistics/predictor/$protein
        $program results_predictor/rosetta/$protein/decoy_$i.pdb proteins/$protein/pdb rmsd >> statistics/predictor/$protein/rosetta_rmsd
        $program results_predictor/rosetta/$protein/decoy_$i.pdb proteins/$protein/pdb gdt >> statistics/predictor/$protein/rosetta_gdt
    done
done

