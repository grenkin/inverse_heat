set term postscript eps enhanced font ', 16' 
set output "theta.eps"
set ylabel "{/Symbol-Oblique q}" norotate font ',20' offset 1.5,0
set xlabel '{/Helvetica-Italic x}' font ',20' offset 0,0.5
set key top right spacing 1.2
plot "t_2.txt" using 1:2 with lines lw 2 lt 3 title "{/Helvetica-Italic t} = 2", \
  "t_5.txt" using 1:2 with lines lw 2 lt 2 title "{/Helvetica-Italic t} = 5", \
  "t_10.txt" using 1:2 with lines lw 2 lt 1 title "{/Helvetica-Italic t} = 10"