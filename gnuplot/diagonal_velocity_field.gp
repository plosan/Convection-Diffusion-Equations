set terminal wxt
plot "output/velocityDiagonal.dat" with image
set xrange [0:1]
set yrange [0:1]
set size ratio 1
unset key
set palette rgb 33,13,10
set cbrange [9:11]
replot


plot "output/velocityDiagonal.dat" with image, "output/velocityDiagonal.dat" using 1:2:($4/75):($5/75) every 20:20 with vectors lc -1 filled
