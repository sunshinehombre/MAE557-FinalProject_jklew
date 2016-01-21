# Gnuplot script file to plot the temporal and spatial convergence of the vorticity-stream function method and the
# G-equation solver.
# This file is called convergenceplotter.p and needs to be placed in FinalProject/data

# set term pngcairo dashed enhanced
set terminal postscript eps size 8,6 enhanced color \
    font 'Helvetica,32' linewidth 1
    
    set autoscale
    unset log
    unset label
    set xtic auto
    set ytic auto
    set size square

    set style line 1 lc rgb "black" lt 2 lw 1 # Dotted line
    set style line 2 lc rgb "blue" pt 5 lt -1 lw 1 # Solid line w filled blue square
    set style line 3 lc rgb "green" pt 7 lt -1 lw 1 # Solid line w filled green circles

    # *** Temporal Convergence ***
    set output "temporal_convergence.eps"
    set logscale x
    set logscale y
    set xrange [1E-8:1E-3]
    set yrange [1E-8:1.0]
    set xlabel "{/Symbol D}t"
    set ylabel "|{/Symbol e}|"
    set key right bottom

    set arrow from 1E-7,1E-6 to 0.0001,1E-3 nohead front ls 1
    set label "\\~ {/Symbol D}t" at 1E-5,1E-5 right

    plot "timeErrVortSF.txt" u 1:2 title "Vorticity-Stream Function" w linespoints ls 2 , \
    "timeErrGEqn.txt" u 1:2 title "G-Eqn Solver" w linespoints ls 3

    # *** Spatial Convergence ***
    unset arrow
    unset label
    set output "spatial_convergence.eps"
    set logscale x
    set logscale y
    set xrange [0.001:1.0]
    set yrange [1E-4:1]
    set xlabel "{/Symbol D}x"
    # set ylabel "|{/Symbol e}|"
    set ylabel "L_{1} Norm"
    set key right top

    # set arrow from 1E-6,1E-4 to 0.001,0.1 nohead front ls 1
    # set label "\\~ {/Symbol D}x" at 1E-5,0.002 right

    set arrow from 0.01,3e-2 to 0.1,3e-1 nohead front ls 1
    set label "\\~ {/Symbol D}x" at 0.03,0.03 right

    set arrow from 0.01,1e-4 to 0.1,1e-2 nohead front ls 1
    set label "\\~ {/Symbol D}x^2" at 0.04,9e-4 left

    plot "spceErrVortSF.txt" u 1:2 title "Vorticity-Stream Function" w linespoints ls 2 , \
    "spceErrGEqn.txt" u 1:2 title "G-Eqn Solver" w linespoints ls 3
