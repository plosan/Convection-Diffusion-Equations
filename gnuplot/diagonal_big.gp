set terminal wxt

plot "input/diagonal_N200_Pe1.0e+00.dat" with image
set xrange [0:1]
set yrange [0:1]
set size ratio 1

set xtics axis
set xtics format "%.1f"
set xlabel ("$x \\ (\\mathrm{m})$")

set ytics axis
set ytics format "%.1f"
set ylabel ("$y \\ (\\mathrm{m})$")

set cblabel ("$\\phi$")
set cbtics format "%.1f"

unset key
set title ("\\textbf{Diagonal case} $(\\mathrm{Pe} = 1)$")
set palette rgb 33,13,10

replot

set terminal epslatex color colortext size 12cm,12cm
set rmargin 2
set output "figures/case_diagonal_flow/diagonal_N200_Pe1.0e+00.tex"
replot
set output
