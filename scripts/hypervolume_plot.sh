#!/bin/bash

python3 scripts/hypervolume_plot.py $1/graphs/frag_hypervolume.png  "FRAG hypervolume" $1/algorithm/frag.hv
python3 scripts/hypervolume_plot.py $1/graphs/res_hypervolume.png "RES hypervolume" $1/algorithm/res.hv
python3 scripts/hypervolume_plot.py $1/graphs/total_hypervolume.png "FRAG+RES Hypervolume" $1/algorithm/frag.hv $1/algorithm/res.hv
