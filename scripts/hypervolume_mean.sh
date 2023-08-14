#!/bin/bash

for i in `seq 1 $2`;
do
    python3 scripts/fitness_convergence.py $1/$i/algorithm/frag_gen_ $3 > $1/$i/algorithm/frag.hv;
    python3 scripts/fitness_convergence.py $1/$i/algorithm/res_gen_ $3 > $1/$i/algorithm/res.hv;
done

python3 scripts/mean_hypervolume.py $1 $2 algorithm/frag.hv > $4/mean_frag.hv
python3 scripts/mean_hypervolume.py $1 $2 algorithm/res.hv > $4/mean_res.hv
