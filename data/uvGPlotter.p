# Gnuplot script file to plot the evolution of u(x,y,t), v(x,y,t), and G(x,y,t)
# This file is called uvGPlotter.p and needs to be placed in FinalProject/data

# set term pngcairo dashed enhanced
set terminal postscript eps size 8,6 enhanced color \
    font 'Helvetica,32' linewidth 1
    
    set autoscale
    unset log
    unset label
    set xtic auto
    set ytic auto
    set size square
    set encoding utf8

    # *** Set variables *************************************************************
    h = 0.015625 # Spatial step
    lside = 1.0 # Length & width of domain
    tss = lside/h # Total spatial steps

    udata = system('find . -maxdepth 1 -name "*1_var*" -print')
    vdata = system('find . -maxdepth 1 -name "*2_var*" -print')
    gdata = system('find . -maxdepth 1 -name "*3_var*" -print')

    # *** Set plot properties *******************************************************
    set xrange[0:tss]
    set yrange[0:tss]
    set xlabel "x [m]"
    set ylabel "y [m]"
    set xtics ("0" 0, sprintf("%3.1f",lside/5.0) tss/5, \
	       sprintf("%3.1f",2*lside/5.0) 2*tss/5, \
	       sprintf("%3.1f",3*lside/5.0) 3*tss/5, \
	       sprintf("%3.1f",4*lside/5.0) 4*tss/5, \
	       sprintf("%3.1f",lside) tss)
    set ytics ("0" 0, sprintf("%3.1f",lside/5.0) tss/5, \
	       sprintf("%3.1f",2*lside/5.0) 2*tss/5, \
	       sprintf("%3.1f",3*lside/5.0) 3*tss/5, \
	       sprintf("%3.1f",4*lside/5.0) 4*tss/5, \
	       sprintf("%3.1f",lside) tss)

    # *** Plot u(x,y,t) *************************************************************
    do for [ufile in udata] {
	unset log
	    unset cblabel
	    # set output sprintf('%s.png', ufile)
	    set output sprintf('%s.eps', ufile)
	    set cbrange[0:25]
	    set cblabel "u(x,y,t) [m/s]"
	    plot ufile matrix with image title ""
}

    # *** Plot v(x,y,t) *************************************************************
    do for [vfile in vdata] {
	unset log
	    unset cblabel
	    # set output sprintf('%s.png', vfile)
	    set output sprintf('%s.eps', vfile)
	    set cbrange[-0.1:0.1]
	    set cblabel "v(x,y,t) [m/s]"
	    plot vfile matrix with image title ""
}

    # *** Plot G(x,y,t) *************************************************************
do for [gfile in gdata] {
    unset log
	unset cblabel
	# set output sprintf('%s.png', gfile)
	set output sprintf('%s.eps', gfile)
	set cbrange [-1:1]
	# set cbtics ('G<G_{o}' -1, 'G=G_{o}' 0, 'G>G_{o}' 1)
	set cbtics ('G<G_{o}' -1, 'G>G_{o}' 1)
	set palette defined (-1 "blue", -0.00001 "blue", 0 "black", 0.00001 "red", 1 "red")
	plot gfile matrix with image title ""
}

