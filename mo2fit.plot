# set terminal pngcairo
set xlabel "Energy"
set ylabel "SS"
set key outside top right
plot "res_0.fit" title '1', \
     "res_1.fit" title '2', \
     "res_2.fit" title '3', \
     "res_3.fit" title '4', \
     "res_4.fit" title '5', \