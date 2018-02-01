set terminal png tiny size 800,800
set output "check.png"
set size 1,1
set grid
unset key
set border 15
set tics scale 0
set xlabel "chr1"
set ylabel "chr1"
set format "%.0f"
if(GPVAL_VERSION < 5) { set mouse clipboardformat "[%.0f, %.0f]" }
set xrange [1:1000000]
set yrange [1:986464]
set style line 1  lt 1 lw 3 pt 6 ps 1
set style line 2  lt 3 lw 3 pt 6 ps 1
set style line 3  lt 2 lw 3 pt 6 ps 1
plot \
 "check.fplot" title "FWD" w lp ls 1, \
 "check.rplot" title "REV" w lp ls 2
