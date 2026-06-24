set terminal svg size 900,500 enhanced background "#ffffff"

set style line 1 lc rgb "#4c72b0" lw 2 pt 7 ps 1.4
set style line 2 lc rgb "#dd8452" lw 2 pt 9 ps 1.4

set grid lt 0 lw 0.5 lc rgb "#cccccc"
set border 3
set tics nomirror
set key top left font "Sans,10"


files  = "sqrt(x) 1_sqrt(x) ln(x)_sqrt(x)"
titles = "sqrt(x)=2/3 1/sqrt(x)=2 ln(x)/sqrt(x)=-4"

do for [i=1:3] {
    name  = word(files,  i)
    title = word(titles, i)
    file  = "data/sweep_".name.".data"


    
    set output "plots/" . name . "_sweep.svg"
    set multiplot layout 1,2 title "N sweep — ".title font "Sans,12"

    set title "Function evaluations vs N"
    set xlabel "N (points per level)"
    set ylabel "ncalls (avg over 8 runs)"
    set logscale x 2
    set logscale y
    set xrange [1.5:260]
    set xtics (2, 4, 8, 16, 32, 64, 128, 256)
    plot file using 1:2 with linespoints ls 1 title "random+CC"

    set title "Actual error vs N"
    set ylabel "act\\_err"
    set logscale y
    plot file using 1:3 with linespoints ls 2 title "random+CC"

    unset multiplot
    unset logscale
}

set output
