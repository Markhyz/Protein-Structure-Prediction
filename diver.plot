# set terminal pngcairo
set xlabel "Generation"
set ylabel "Diversity"
set key top left
plot "res.diver" using 1:2 title 'Diversity' with lines