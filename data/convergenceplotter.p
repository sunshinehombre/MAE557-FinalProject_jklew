# Gnuplot script file to plot the temporal and spatial convergence of the vorticity-stream function method
# This file is called convergenceplotter.p and needs to be placed in FinalProject/data

set term pngcairo dashed enhanced
    
    set autoscale
    unset log
    unset label
    set xtic auto
    set ytic auto
    set key font ",10"
    set size square

    set style line 1 lc rgb "black" pt 7 lt 0 lw 1
    set style line 2 lc rgb "blue" pt 7 lt 0 lw 1
    set style line 3 lc rgb "black" lt -1 lw 1
    set style line 4 lc rgb "blue" lt -1 lw 1
    set style line 5 lc rgb "green" lt -1 lw 1
    set style line 6 lc rgb "orange" lt -1 lw 1
    set style line 7 lc rgb "cyan" lt -1 lw 1
    set style line 8 lc rgb "violet" lt -1 lw 1
    set style line 9 lc rgb "red" lt -1 lw 1
    set style line 10 lc rgb "black" lt 3 lw 1
    set style line 11 lc rgb "blue" pt 5 lt -1 lw 1
    set style line 12 lc rgb "green" pt 7 lt -1 lw 1
    set style line 13 lc rgb "red" pt 9 lt -1 lw 1
    set style line 14 lc rgb "grey" pt 13 lt -1 lw 1
    set style line 15 lc rgb "black" lt 2 lw 1  
    set style line 16 lc rgb "violet" pt 13 lt -1 lw 1
    set style line 17 lc rgb "blue" lt 3 lw 1
    set style line 18 lc rgb "blue" lt 4 lw 1

    # # *** Temporal Convergence ***
    # set output "temporal_convergence.png"
    # set logscale x
    # set logscale y
    # set xrange [1E-6:1E-3]
    # set yrange [1E-6:1E-1]
    # set xlabel "{/Symbol D}t"
    # set ylabel "|{/Symbol e}|"
    # set key right bottom

    # set arrow from 1E-6,1E-4 to 0.001,0.1 nohead front ls 15
    # set label "\\~ {/Symbol D}t" at 1E-5,0.002 right

    # plot "timeErrVortSF.txt" u 1:2 title "Vorticity-Stream Function" w linespoints ls 11 , \
    # "timeErrGEqn.txt" u 1:2 title "G-Eqn Solver" w linespoints ls 12

    # *** Spatial Convergence ***
    unset arrow
    unset label
    set output "spatial_convergence.png"
    set logscale x
    set logscale y
    set xrange [0.001:1.0]
    set yrange [1E-12:1]
    set xlabel "{/Symbol D}x"
    # set ylabel "|{/Symbol e}|"
    set ylabel "L_{1} Norm"
    set key right bottom

    # set arrow from 1E-6,1E-4 to 0.001,0.1 nohead front ls 15
    # set label "\\~ {/Symbol D}x" at 1E-5,0.002 right

    set arrow from 0.01,1e-4 to 0.1,1e-3 nohead front ls 15
    set label "\\~ {/Symbol D}x" at 0.03,0.001 right

    set arrow from 0.01,1e-8 to 0.1,1e-6 nohead front ls 15
    set label "\\~ {/Symbol D}x^2" at 0.02,5e-7 left

    plot "spceErrVortSF.txt" u 1:2 title "Vorticity-Stream Function" w linespoints ls 11 , \
    "spceErrGEqn.txt" u 1:2 title "G-Eqn Solver" w linespoints ls 12
