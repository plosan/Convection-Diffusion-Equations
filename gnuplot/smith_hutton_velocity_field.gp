set terminal wxt

set xrange [-1:1]
set yrange [0:1]
set size ratio 0.5

set xtics axis
set xtics format "%.1f"
set xlabel ("$x \\ (\\mathrm{m})$")

set ytics border out nomirror
set ylabel ("$y \\ (\\mathrm{m})$")

set cblabel ("$\\norm{\\vb{v}} \\ (\\mathrm{m} / \\mathrm{s})$")
set cbtics format "%.1f"

unset key
set title ("\\textbf{Smith--Hutton case -- Velocity field}")
set palette rgb 33,13,10

plot "input/smith_hutton_velocity_mod.dat" with image, "input/smith_hutton_velocity_vec.dat" using 1:2:($4/sqrt(($4)**2 + ($5)**2)/30):($5/sqrt(($4)**2 + ($5)**2)/30) every 10:10 with vectors lc -1 filled

replot

set terminal epslatex color colortext size 13cm,7cm
set rmargin 2
set output "figures/case_smith_hutton/smith_hutton_velocity_field.tex"
replot
unset output
unset terminal
