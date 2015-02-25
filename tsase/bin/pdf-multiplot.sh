#!/usr/bin/bash

gnuplot << EOF
reset
set terminal postscript eps color enhanced "Helvetica" 20
set output "pdf.eps"
set title "Pair distribution of all atoms"
#set nokey
set xlabel "Angstroms"
set ylabel "g(r)"

set style line 1  lc rgb '#a6cee3' lt 1 lw 1.5 # --- blue
set style line 2  lc rgb '#1f78b4' lt 1 lw 1.5 #      .
set style line 3  lc rgb '#b2df8a' lt 1 lw 1.5 #      .
set style line 4  lc rgb '#33a02c' lt 1 lw 1.5 #      .
set style line 5  lc rgb '#fb9a99' lt 1 lw 1.5 #      .
set style line 6  lc rgb '#e31a1c' lt 1 lw 1.5 #      .
set style line 7  lc rgb '#fdbf6f' lt 1 lw 1.5 #      .
set style line 8  lc rgb '#ff7f00' lt 1 lw 1.5 #      .
set style line 9  lc rgb '#cab2d6' lt 1 lw 1.5 #      .
set style line 10 lc rgb '#6a3d9a' lt 1 lw 1.5 #      .
set style line 11 lc rgb '#ffff99' lt 1 lw 1.5 #      .
set style line 12 lc rgb '#b15928' lt 1 lw 1.5 # --- green
set style line 13 lc rgb '#66c2a5' lt 1 lw 1.5 # --- green
set style line 14 lc rgb '#fc8d62' lt 1 lw 1.5 # --- green
set style line 15 lc rgb '#8da0cb' lt 1 lw 1.5 # --- green
set style line 16 lc rgb '#e78ac3' lt 1 lw 1.5 # --- green
set style line 17 lc rgb '#a6d854' lt 1 lw 1.5 # --- green
set style line 18 lc rgb '#ffd92f' lt 1 lw 1.5 # --- green
set style line 19 lc rgb '#e5c494' lt 1 lw 1.5 # --- green
set style line 20 lc rgb '#b3b3b3' lt 1 lw 1.5 # --- green


plot for [i=1:20] 'interprdf-'.i.'.dat'  u 1:2 w l ls i title 'CON'.i

EOF
