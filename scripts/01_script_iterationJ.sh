#!/usr/bin/gnuplot
set key Left 
#set key left bottom

forceLog = 1			# 0: nothing		1: forces log y by setting y -> |y|

if(forceLog) {
	set log y
}
pointTypeIndex=0
set pointsize 2
set title system("pwd")
FILES = system("ls -1 iteration-info-B*.dat")
do for [i=1:words(FILES)] {
	set style line i lc i
}

#f(x) = (x);	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
g(x) = (x == 0 ? NaN : abs(x));	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
h(x) = (x > 0 ? 2 : 10)

if(forceLog) {
	plot for [i=1:words(FILES)] word(FILES,i) u 1:(g($4)):(h($4)):1 title word(FILES,i)."   × pos, ▿ neg" with linespoints ls i pt variable
}
else {
	plot\
		for [i=1:words(FILES)] word(FILES,i) u 1:4 title word(FILES,i) with linespoints ls i
}

pause -1

