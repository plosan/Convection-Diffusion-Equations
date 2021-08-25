set terminal wxt

plot "input/smith_hutton_N201_Pe1.0e-06.dat" with image
set xrange [-1:1]
set yrange [0:1]
set size ratio 0.5

set xtics axis
set xtics format "%.1f"
set xlabel ("$x \\ (\\mathrm{m})$")

set ylabel ("$y \\ (\\mathrm{m})$")

set cblabel ("$\\phi$")
set cbtics format "%.1f"

unset key
set title ("\\textbf{Smith--Hutton case} $(\\mathrm{Pe} = 10^{-6})$")
set palette rgb 33,13,10

replot

set terminal epslatex color colortext size 13cm,7cm
set rmargin 2
set output "figures/case_smith_hutton/smith_hutton_N201_Pe1.0e-06.tex"
replot
unset output
unset terminal
