set terminal svg background rgb "white"
set output "erf_out.svg"


set title "erf(x) comparison"
set xlabel "x"
set ylabel "erf(x)"
set key left top
plot \
	"erf_out.data" using 1:2 with linespoints title "my erf(x)", \
	"erf_out.data" using 1:3 with lines title "tabulated erf(x)", \
	"erf_out.data" using 1:4 with lines title "abs diff"