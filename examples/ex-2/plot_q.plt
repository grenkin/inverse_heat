set term postscript eps enhanced font ', 18' 
set output "q.eps"
set ylabel "{/Helvetica-Italic q}({/Helvetica-Italic t})" norotate font ',20' offset 2,0
set xlabel '{/Helvetica-Italic t}' font ',20'
set yrange [0.22:0.37]
unset key
plot "output_q.txt" using 1:2 with lines linewidth 3