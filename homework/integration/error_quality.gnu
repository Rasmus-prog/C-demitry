set terminal svg background rgb "white"
set output "error_quality.svg"

set key left top
set logscale x
set logscale y
set format x "10^{%L}"
set format y "10^{%L}"
set xlabel "estimated error"
set ylabel "actual error"
set grid
set title "Estimated error vs actual error"

plot \
	"error_quality.data" every ::0::0 using 5:6 with points pt 7 ps 1.5 title "1/sqrt(x) on [0,1]", \
	"error_quality.data" every ::1::1 using 5:6 with points pt 9 ps 1.5 title "log(x)/sqrt(x) on [0,1]", \
	"error_quality.data" every ::2::2 using 5:6 with points pt 5 ps 1.5 title "exp(-x) on [0,inf)", \
	"error_quality.data" every ::3::3 using 5:6 with points pt 13 ps 1.5 title "exp(-x*x) on [0,inf)", \


