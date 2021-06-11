#!/bin/bash

python3 scripts/diversity_plot.py $1/graphs/frag_diversity.png $1/algorithm/frag.diversity "FRAG Diversity"
python3 scripts/diversity_plot.py $1/graphs/res_diversity.png $1/algorithm/res.diversity "RES Diversity"
python3 scripts/diversity_plot.py $1/graphs/total_diversity.png $1/algorithm/total.diversity "FRAG+RES Diversity"
