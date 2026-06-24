# Two-panel SVG per integral: ncalls and act_err per method.

set terminal svg size 1200,550 enhanced background "#ffffff"


set style fill solid 0.85 border -1
set boxwidth 0.7
set grid ytics lt 0 lw 0.5 lc rgb "#cccccc"
set border 3
set tics nomirror

files  = "sqrt(x) 1_sqrt(x) ln(x)_sqrt(x)"
titles = "sqrt(x)=2/3 1/sqrt(x)=2 ln(x)/sqrt(x)=-4"

# Method labels match the first column of the .data files (rows 1-10)
set xtics ("fixed" 0, "fixed CC" 1, \
           "rand N2" 2, "rand+CC N2" 3, \
           "rand N4" 4, "rand+CC N4" 5, \
           "rand N8" 6, "rand+CC N8" 7, \
           "rand N16" 8, "rand+CC N16" 9) \
    rotate by 90 nomirror right \


do for [i=1:3] {
    name  = word(files,  i)
    title = word(titles, i)
    file  = "data/".name.".data"

    set output "plots/" . name . "_comparison.svg"
    set multiplot layout 1,2 title title font "Sans,13"

    # Left panel: ncalls
    set title "Function evaluations (ncalls)"
    set ylabel "ncalls"
    set logscale y
    set yrange [1:*]
    unset key
    plot file using 2:xtic(1) with boxes lc rgb "#4c72b0" notitle

    # Right panel: actual error
    set title "Actual error  |result - exact|"
    set ylabel "act\\_err"
    set yrange [*:*]
    set logscale y
    plot file using 3:xtic(1) with boxes lc rgb "#dd8452" notitle

    unset multiplot
    unset logscale
    unset yrange
}

set output
