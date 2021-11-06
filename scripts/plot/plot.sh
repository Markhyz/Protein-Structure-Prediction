#!/bin/bash

proteins=(1ab1 1acw 1ail 1aly 1bdd 1crn 1dfn 1enh 1gb1 1hhp 1i6c 1rop 1zdd 2kdl 2mr9 2p81 T0868 T0900 T0968s1 T1010)

for protein in ${proteins[@]};
do
    python3 scripts/plot/fitness_plot.py statistics/$1/$protein/${protein}_frag_scatter_1_full.$2 "${protein^^} fitness (FRAG decoder gen 1)" 1 statistics/$1/$protein frag_first_frontier
    python3 scripts/plot/fitness_plot.py statistics/$1/$protein/${protein}_frag_scatter_1000_full.$2 "${protein^^} fitness (FRAG decoder gen 1000)" 1 statistics/$1/$protein frag_last_frontier 
    python3 scripts/plot/fitness_plot.py statistics/$1/$protein/${protein}_frag_scatter_1_zoom.$2 "${protein^^} fitness (FRAG decoder gen 1)" 0 statistics/$1/$protein frag_first_frontier
    python3 scripts/plot/fitness_plot.py statistics/$1/$protein/${protein}_frag_scatter_1000_zoom.$2 "${protein^^} fitness (FRAG decoder gen 1000)" 0 statistics/$1/$protein frag_last_frontier 

    python3 scripts/plot/fitness_plot.py statistics/$1/$protein/${protein}_res_scatter_1_full.$2 "${protein^^} fitness (RES decoder gen 1)" 1 statistics/$1/$protein res_first_frontier
    python3 scripts/plot/fitness_plot.py statistics/$1/$protein/${protein}_res_scatter_1000_full.$2 "${protein^^} fitness (RES decoder gen 1000)" 1 statistics/$1/$protein res_last_frontier
    python3 scripts/plot/fitness_plot.py statistics/$1/$protein/${protein}_res_scatter_1_zoom.$2 "${protein^^} fitness (RES decoder gen 1)" 0 statistics/$1/$protein res_first_frontier
    python3 scripts/plot/fitness_plot.py statistics/$1/$protein/${protein}_res_scatter_1000_zoom.$2 "${protein^^} fitness (RES decoder gen 1000)" 0 statistics/$1/$protein res_last_frontier

    python3 scripts/plot/hypervolume_plot.py statistics/$1/$protein/${protein}_hv.$2 "${protein^^} hypervolume" statistics/$1/$protein/mean_frag.hv statistics/$1/$protein/mean_res.hv
    python3 scripts/plot/diversity_plot.py statistics/$1/$protein/${protein}_diversity.$2 statistics/$1/$protein/mean.diversity "${protein^^} diversity"
done

