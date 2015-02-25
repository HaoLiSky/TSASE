#!/usr/bin/bash

gnuplot << EOF
reset
set terminal postscript eps color enhanced "Helvetica" 20
set output "pdf.eps"
set title "Pair distribution"
set nokey
set xlabel "Angstroms"
set ylabel "g(r)"

plot "pdf.dat"  u 1:2 w p lt 3 lw 2.0 pt 7 ps 1.0 lc -1, \
    "interprdf.dat"  u 1:2 w l lt 1 lw 2.4 lc 3 

EOF
