#!/bin/bash

proteins=(1acw  1ail  1crn  1enh  1rop  1zdd  2mr9  2p81)

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

