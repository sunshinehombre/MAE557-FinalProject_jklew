# Gnuplot script file to plot the primary and secondary conservation of the vorticity-stream function method.
# This file is called conservationplotter.p and needs to be placed in FinalProject/data

# set term pngcairo dashed enhanced
set terminal postscript eps size 8,6 enhanced color \
    font 'Helvetica,32' linewidth 1
    
    set autoscale
    unset log
    unset label
    set xtic auto
    set ytic auto
    set key font ",10"
    set size square
    set encoding utf8

    set style line 1 lc rgb "red" lt -1 lw 1 # Solid line

    tottime = 1.0

    # *** Primary Conservation ***
    unset arrow
    unset label
    set output "primary.eps"
    set xrange [0:tottime]
    set yrange [-4e-10:4e-10]
    set xlabel "t [s]"
    set ylabel "{/Symbol S} u_{j,k}^{n+1} - u_{j,k}^n"
    
    plot 'primary.txt' u 1:2 every 100 title "" w lines ls 1
    
    # *** Secondary Conservation ***
    unset arrow
    unset label
    set output "secondary.eps"
    set xrange [0:tottime]
    set yrange [-4e-8:4e-8]
    set xlabel "t [s]"
    set ylabel "{/Symbol S} (u_{j,k}^{n+1}+v_{j,k}^{n+1})^2 - (u_{j,k}^n+v_{j,k}^n)^2"

    plot 'secondary.txt' u 1:2 every 100 title "" w lines ls 1
