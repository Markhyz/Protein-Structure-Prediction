#!/bin/bash

proteins=(1ab1 1acw 1ail 1aly 1bdd 1crn 1dfn 1enh 1gb1 1hhp 1i6c 1rop 1zdd 2kdl 2mr9 2p81 T0868 T0900 T0968s1 T1010)

program="$ROSETTA_BUILD_DIR/score_calc"
mv ./build/score_calc $ROSETTA_BUILD_DIR/score_calc 2>/dev/null

for protein in ${proteins[@]};
do
    for i in `seq 1 20`;
    do
        external/mufold -P results_predictor/$protein/$i/decoys/res/
        selected_decoy=$(head -n 10 res.log | grep -Po "^INFO.*1 :[^res]*\Kres_decoy_[0-9]+.pdb")
        $program results_predictor/$protein/$i/decoys/res/$selected_decoy proteins/$protein/pdb rmsd >> statistics/predictor/$protein/mufold_rmsd
        $program results_predictor/$protein/$i/decoys/res/$selected_decoy proteins/$protein/pdb gdt >> statistics/predictor/$protein/mufold_gdt
	rm res.log
    done
done

