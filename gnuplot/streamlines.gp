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

set style line 1 \
    linecolor rgb '#000000' \
    linetype -1 linewidth 1 \
    pointtype 7 pointsize 0

scale = 50;

plot    'input/smith_hutton_velocity_mod.dat' with image, \
        'input/smith_hutton_streamline_0.10.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.20.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.30.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.40.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.50.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.60.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.70.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.80.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.90.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_0.99.dat' with linespoints linestyle 1, \
        'input/smith_hutton_streamline_arrow_0.10.dat' using 1:2:($3/scale):($4/scale) every 250:250 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.20.dat' using 1:2:($3/scale):($4/scale) every 250:250 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.30.dat' using 1:2:($3/scale):($4/scale) every 200:200 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.40.dat' using 1:2:($3/scale):($4/scale) every 200:200 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.50.dat' using 1:2:($3/scale):($4/scale) every 200:200 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.60.dat' using 1:2:($3/scale):($4/scale) every 200:200 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.70.dat' using 1:2:($3/scale):($4/scale) every 200:200 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.80.dat' using 1:2:($3/scale):($4/scale) every 200:200 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.90.dat' using 1:2:($3/scale):($4/scale) every 200:200 with vectors lc -1 filled, \
        'input/smith_hutton_streamline_arrow_0.99.dat' using 1:2:($3/scale):($4/scale) every 200:200 with vectors lc -1 filled

set terminal epslatex color colortext size 13cm,7cm
set rmargin 2
set output "figures/case_smith_hutton/smith_hutton_N201_streamlines.tex"
replot
unset output
unset terminal
