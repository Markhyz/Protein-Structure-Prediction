#!/bin/bash

proteins=(1ab1 1acw 1ail 1aly 1bdd 1crn 1dfn 1enh 1gb1 1hhp 1i6c 1rop 1zdd 2kdl 2mr9 2p81 T0868 T0900 T0968s1 T1010)

for protein in ${proteins[@]};
do
    best_run=$(python3 scripts/best_score.py $protein)
    external/mufold -P results_predictor/$protein/$best_run/decoys/res/
    selected_decoy=$(head -n 10 res.log | grep -Po "^INFO.*1 :[^res]*\Kres_decoy_[0-9]+.pdb")
    cp results_predictor/$protein/$best_run/decoys/res/$selected_decoy statistics/predictor/$protein/mufold.pdb
    echo "$protein $selected_decoy"
    rm res.log
done

