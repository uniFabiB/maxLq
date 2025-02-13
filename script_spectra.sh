#!/usr/bin/gnuplot
subfiles=7			# number of files to be grouped

set key Left			# legend left aligned text
set key left bottom		# legend position
set log x
set log y
set xrange [5:*]		# x range starts at 5
set pointsize 2			# marker size

pointTypeIndex=0
FILES = system("ls -1 *spectrum*.dat")
do for [i=1:words(FILES)] {
	do for [j=1:subfiles] {
		if(i % subfiles == j % subfiles) {
			set style line i lc j		# set lc = line color
		}
	}
	if((i-1) % subfiles == 0){
		pointTypeIndex = pointTypeIndex+2	# +2 to get better differentiable markers
	}
	set style line i pt pointTypeIndex		# set pt = point type
}
plot for [i=1:words(FILES)] word(FILES,i) u 1:2 title word(FILES,i) with linespoints ls i
#plot for [i=1:words(FILES)] word(FILES,i) u 1:($2*($1**4.25)) title word(FILES,i) with linespoints ls i	#calculating example

pause -1			# keeps window open
