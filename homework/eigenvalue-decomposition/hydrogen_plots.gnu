set datafile commentschars "#"
set terminal svg background rgb "#ffffff"
set grid

set output "convergence_dr.svg"
set title "Ground-state energy convergence vs dr"
set xlabel "dr [Bohr]"
set ylabel "Energy [Hartree]"
plot \
  "hydrogen_convergence_dr.txt" using 1:2 with linespoints lw 2 pt 7 title "E0 numerical", \
  "hydrogen_convergence_dr.txt" using 1:3 with lines lw 2 dt 2 title "E0 exact (-0.5)"

set output "convergence_rmax.svg"
set title "Ground-state energy convergence vs rmax"
set xlabel "rmax [Bohr]"
set ylabel "Energy [Hartree]"
plot \
  "hydrogen_convergence_rmax.txt" using 1:2 with linespoints lw 2 pt 7 title "E0 numerical", \
  "hydrogen_convergence_rmax.txt" using 1:3 with lines lw 2 dt 2 title "E0 exact (-0.5)"

set output "energies.svg"
set title "Lowest s-wave energies: numerical vs exact"
set xlabel "n"
set ylabel "Energy [Hartree]"
plot \
  "hydrogen_energies.txt" using 1:2 with points pt 7 ps 1.2 title "Numerical", \
  "hydrogen_energies.txt" using 1:3 with linespoints pt 5 lw 2 title "Exact"

set output "wavefunctions.svg"
set title "Reduced radial wavefunctions f_n(r)"
set xlabel "r [Bohr]"
set ylabel "f_n(r)"
plot \
  "hydrogen_wavefunctions.txt" using 1:2 with lines lw 2 title "n=1 numerical", \
  "hydrogen_wavefunctions.txt" using 1:3 with lines lw 2 dt 2 title "n=1 exact", \
  "hydrogen_wavefunctions.txt" using 1:4 with lines lw 2 title "n=2 numerical", \
  "hydrogen_wavefunctions.txt" using 1:5 with lines lw 2 dt 2 title "n=2 exact", \
  "hydrogen_wavefunctions.txt" using 1:6 with lines lw 2 title "n=3 numerical", \
  "hydrogen_wavefunctions.txt" using 1:7 with lines lw 2 dt 2 title "n=3 exact"
