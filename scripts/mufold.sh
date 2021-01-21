#!/bin/bash

proteins=(1crn  1gb1  1hhp  1i6c  1vii  T0868  T0900)

program="$ROSETTA_BUILD_DIR/score_calc"
mv ./build/score_calc $ROSETTA_BUILD_DIR/score_calc 2>/dev/null

for protein in ${proteins[@]};
do
    for i in `seq 1 32`;
    do
        external/mufold -P results/$protein/decoys/$i/
        selected_decoy=$(head -n 5 $i.log | grep -Po "^INFO[ ]+:[ ]+1[ ]+:[^decoy]*\Kdecoy_[0-9]+.pdb")
        $program results/$protein/decoys/$i/$selected_decoy proteins/$protein/pdb rmsd >> results/$protein/mufold_rmsd
        $program results/$protein/decoys/$i/$selected_decoy proteins/$protein/pdb gdt >> results/$protein/mufold_gdt
    done
done

