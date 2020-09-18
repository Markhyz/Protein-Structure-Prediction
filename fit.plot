# set terminal pngcairo
set xlabel "Generation"
set ylabel "Fitness"
set key top left
set autoscale xfixmin
set autoscale xfixmax
plot "res.fit" using 1:2 title 'Best fitness' with lines, "res.fit" using 1:3 title 'Mean Fitness'  with lines