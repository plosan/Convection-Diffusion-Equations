set terminal epslatex color colortext size 15cm,15cm
set output "figures/case_diagonal_flow/multiplot.tex"

set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xlabel ("$x \\ \\mathrm{m}$")
set ylabel ("$y \\ \\mathrm{m}$")
set cblabel ("$\\phi$")
set xtics axis
set ytics axis
unset key
set palette rgb 33,13,10

set multiplot layout 1,2 rowsfirst

set title ("\\texbt{Diagonal flow} ($\\mathrm{Pe} = 1$)")
plot "input/diagonal_N200_Pe1.0e+00.dat" with image

set title ("\\texbt{Diagonal flow} ($\\mathrm{Pe} = 10^3$)")
plot "input/diagonal_N200_Pe1.0e+03.dat" with image

unset multiplot
unset output
unset terminal
