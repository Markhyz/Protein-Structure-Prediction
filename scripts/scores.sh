#!/bin/bash

#proteins=(1crn  1gb1  1hhp  1i6c  1vii  T0868  T0900)
#proteins=(1acw  1ail  1crn  1enh  1rop  1zdd  2mr9  2p81)
proteins=(1ab1 1acw 1ail 1aly 1bdd 1crn 1dfn 1enh 1gb1 1hhp 1i6c 1rop 1zdd 2kdl 2mr9 2p81 T0868 T0900 T0968s1 T1010)

for protein in ${proteins[@]};
do
      mkdir -p statistics/$2/$protein

     	cat $2/$protein/*/scores/res.ca_rmsd > statistics/predictor/$protein/best.rmsd
    	cat $2/$protein/*/scores/res.gdtts > statistics/predictor/$protein/best.gdt
    	cat $2/$protein/*/scores/res.mean_rmsd > statistics/predictor/$protein/mean.rmsd
    	cat $2/$protein/*/scores/res.mean_gdt > statistics/predictor/$protein/mean.gdt

      cp $2/$protein/1/algorithm/frag_gen_1.norm_fitness statistics/$2/$protein/frag_first_frontier.norm_fitness
      cp $2/$protein/1/algorithm/frag_gen_$4.norm_fitness statistics/$2/$protein/frag_last_frontier.norm_fitness
      cp $2/$protein/1/algorithm/res_gen_1.norm_fitness statistics/$2/$protein/res_first_frontier.norm_fitness
      cp $2/$protein/1/algorithm/res_gen_$4.norm_fitness statistics/$2/$protein/res_last_frontier.norm_fitness

      python3 scripts/mean_diversity.py $2/$protein $3 algorithm/total.diversity > statistics/$2/$protein/mean.diversity
      
      ./scripts/hypervolume_mean.sh $2/$protein $3 $4 statistics/$2/$protein
done

