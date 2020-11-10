# set terminal pngcairo
set xlabel "Energy"
set ylabel "SS"
set key outside top right
plot "1gb1_3000.fit" title '3000', \
     "1gb1_4000.fit" title '4000', \
     "1gb1_2000.fit" title '2000', \
     "1gb1_1000.fit" title '1000', \