#!/bin/bash

python3 scripts/fitness_plot.py $1/graphs/frag_scatter_zoom.png "Fitness (FRAG decoder, gen 1 and 1000)" 0 $1/algorithm frag_gen_1 frag_gen_1000
python3 scripts/fitness_plot.py $1/graphs/frag_scatter_full.png "Fitness (FRAG decoder, gen 1 and 1000)" 1 $1/algorithm frag_gen_1 frag_gen_1000
python3 scripts/fitness_radviz.py $1/algorithm/frag_gen_1 $1/graphs/frag_radviz_1.png "Fitness (FRAG decoder, gen 1)"
python3 scripts/fitness_radviz.py $1/algorithm/frag_gen_1000 $1/graphs/frag_radviz_1000.png "Fitness (FRAG decoder, gen 1000)"

python3 scripts/fitness_plot.py $1/graphs/res_scatter_zoom.png "Fitness (RES decoder, gen 1 and 1000)" 0 $1/algorithm res_gen_1 res_gen_1000
python3 scripts/fitness_plot.py $1/graphs/res_scatter_full.png "Fitness (RES decoder, gen 1 and 1000)" 1 $1/algorithm res_gen_1 res_gen_1000
python3 scripts/fitness_radviz.py $1/algorithm/res_gen_1 $1/graphs/res_radviz_1.png "Fitness (RES decoder, gen 1)"
python3 scripts/fitness_radviz.py $1/algorithm/res_gen_1000 $1/graphs/res_radviz_1000.png "Fitness (RES decoder, gen 1000)"

python3 scripts/fitness_plot.py $1/graphs/frag_res_scatter_zoom.png "Fitness (FRAG and RES decoders, gen 1000)" 0 $1/algorithm frag_gen_1000 res_gen_1000
python3 scripts/fitness_plot.py $1/graphs/frag_res_scatter_full.png "Fitness (FRAG and RES decoder, gen 1000)" 1 $1/algorithm frag_gen_1000 res_gen_1000