# Gnuplot script file to plot the primary and secondary conservation of the vorticity-stream function method.
# This file is called conservationplotter.p and needs to be placed in FinalProject/data

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
    set style line 16 lc rgb "blue" lt -1 lw 1
    set style line 17 lc rgb "blue" lt 3 lw 1
    set style line 18 lc rgb "blue" lt 4 lw 1

    tottime = 0.3

    # *** Primary Conservation ***
    unset arrow
    unset label
    set output "primary.png"
    set xrange [0:tottime]
    set yrange [-4e-10:4e-10]
    set xlabel "t [s]"
    set ylabel "{/Symbol S} u_{j,k}^{n+1} - u_{j,k}^n"
    
    plot 'primary.txt' u 1:2 every 100 title "" w lines ls 9
    
    # *** Secondary Conservation ***
    unset arrow
    unset label
    set output "secondary.png"
    set xrange [0:tottime]
    set yrange [-4e-8:4e-8]
    set xlabel "t [s]"
    set ylabel "{/Symbol S} (u_{j,k}^{n+1}+v_{j,k}^{n+1})^2 - (u_{j,k}^n+v_{j,k}^n)^2"

    plot 'secondary.txt' u 1:2 every 100 title "" w lines ls 9
