#!/usr/bin/gnuplot
#set key Left
#set key left top
#set log x
#set log y
set pointsize 2
FILES = system("ls -1 *brent_info*.dat")

plot for [i=1:words(FILES):1] word(FILES,i) u 2:3 title word(FILES,i) with linespoints lc i pointtype 2

pause -1
