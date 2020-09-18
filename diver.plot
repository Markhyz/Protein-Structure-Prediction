# set terminal pngcairo
set xlabel "Generation"
set ylabel "Diversity"
set key top left
set autoscale xfixmin
set autoscale xfixmax
set yrange [0:0.5]
plot "res.diver" using 1:2 title 'Diversity' with lines