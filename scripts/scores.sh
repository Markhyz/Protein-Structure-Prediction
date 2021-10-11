#!/bin/bash

#proteins=(1crn  1gb1  1hhp  1i6c  1vii  T0868  T0900)
#proteins=(1acw  1ail  1crn  1enh  1rop  1zdd  2mr9  2p81)
proteins=(1ab1 1acw 1ail 1aly 1bdd 1crn 1dfn 1enh 1gb1 1hhp 1i6c 1rop 1zdd 2kdl 2mr9 2p81 T0868 T0900 T0968s1 T1010)

for protein in ${proteins[@]};
do
     	cat results_predictor/$protein/*/scores/res.ca_rmsd > statistics/predictor/$protein/best_rmsd
    	cat results_predictor/$protein/*/scores/res.gdtts > statistics/predictor/$protein/best_gdt
    	cat results_predictor/$protein/*/scores/res.mean_rmsd > statistics/predictor/$protein/mean_rmsd
    	cat results_predictor/$protein/*/scores/res.mean_gdt > statistics/predictor/$protein/mean_gdt
   ##for i in `seq 1 20`;
    #do
   # 	cat results_predictor/$protein/$i/scores/res.ca_rmsd >> statistics/predictor/best_rmsd
   # 	cat results_predictor/$protein/$i/scores/res.ca_gdtts >> statistics/predictor/best_gdt
   # 	cat results_predictor/$protein/$i/scores/res.mean_rmsd >> statistics/predictor/mean_rmsd
   # 	cat results_predictor/$protein/$i/scores/res.mean_gdt >> statistics/predictor/mean_gdt
   # done
done

