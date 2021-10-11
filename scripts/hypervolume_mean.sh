#!/bin/bash

for i in `seq 1 $2`;
do
    for j in `seq 1 $3`;
    do
        frag_hv=`python3 scripts/fitness_convergence.py $1/$i/algorithm/frag_gen_$j`
        res_hv=`python3 scripts/fitness_convergence.py $1/$i/algorithm/res_gen_$j`
        echo $j $frag_hv >> $1/$i/algorithm/frag.hv
        echo $j $res_hv >> $1/$i/algorithm/res.hv
    done
done

python3 scripts/mean_hypervolume.py $1 $2 algorithm/frag.hv > $4/mean_frag.hv
python3 scripts/mean_hypervolume.py $1 $2 algorithm/res.hv > $4/mean_res.hv
