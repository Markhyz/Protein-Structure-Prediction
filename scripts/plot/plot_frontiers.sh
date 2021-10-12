#!/bin/bash

#python3 scripts/fitness_plot.py $1/graphs/frag_scatter_zoom.png "Fitness (FRAG decoder, gen 1 and 1000)" 0 $1/algorithm frag_gen_1 frag_gen_1000
python3 scripts/plot/fitness_plot.py $1/graphs/frag_scatter_1_full.pdf "Fitness (FRAG decoder gen 1)" 1 $1/algorithm frag_gen_1
python3 scripts/plot/fitness_plot.py $1/graphs/frag_scatter_1000_full.pdf "Fitness (FRAG decoder gen 1000)" 1 $1/algorithm frag_gen_1000 
python3 scripts/plot/fitness_plot.py $1/graphs/frag_scatter_1_zoom.pdf "Fitness (FRAG decoder gen 1)" 0 $1/algorithm frag_gen_1
python3 scripts/plot/fitness_plot.py $1/graphs/frag_scatter_1000_zoom.pdf "Fitness (FRAG decoder gen 1000)" 0 $1/algorithm frag_gen_1000 

#python3 scripts/fitness_plot.py $1/graphs/res_scatter_zoom.png "Fitness (RES decoder, gen 1 and 1000)" 0 $1/algorithm res_gen_1 res_gen_1000
python3 scripts/plot/fitness_plot.py $1/graphs/res_scatter_1_full.pdf "Fitness (RES decoder gen 1)" 1 $1/algorithm res_gen_1
python3 scripts/plot/fitness_plot.py $1/graphs/res_scatter_1000_full.pdf "Fitness (RES decoder gen 1000)" 1 $1/algorithm res_gen_1000
python3 scripts/plot/fitness_plot.py $1/graphs/res_scatter_1_zoom.pdf "Fitness (RES decoder gen 1)" 0 $1/algorithm res_gen_1
python3 scripts/plot/fitness_plot.py $1/graphs/res_scatter_1000_zoom.pdf "Fitness (RES decoder gen 1000)" 0 $1/algorithm res_gen_1000

#python3 scripts/fitness_plot.py $1/graphs/frag_res_scatter_zoom.png "Fitness (FRAG and RES decoders, gen 1000)" 0 $1/algorithm frag_gen_1000 res_gen_1000
#python3 scripts/fitness_plot.py $1/graphs/frag_res_scatter_full.png "Fitness (FRAG and RES decoder, gen 1000)" 1 $1/algorithm frag_gen_1000 res_gen_1000
