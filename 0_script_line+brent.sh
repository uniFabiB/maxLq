#!/usr/bin/gnuplot
#set key Left
#set key left top
set log x
#set log y
FILES_brent = system("ls -1 *brent_info*.dat")
FILES_brak = system("ls -1 *lineMin_info*.dat")
do for [i=1:words(FILES_brent)] {
	set style line i lc i		# set lc = line color
}

plot \
for [i=1:words(FILES_brent):1] \
	word(FILES_brent,i) u 2:3 title word(FILES_brent,i)." brent" with linespoints ls i pointtype 2 pointsize 2, \
for [i=1:words(FILES_brak):1] \
	word(FILES_brak,i) u 3:6 title word(FILES_brak,i)." brak" with linespoints ls i pointtype 4 pointsize 1

pause -1
