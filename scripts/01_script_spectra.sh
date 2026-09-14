#!/usr/bin/gnuplot
#FILES = system("ls -1 *spectrum*.dat")		# all
uFiles = system("ls -1 u*spectrum*.dat")	# just uvec
#FILES = system("ls -1 gradJ*spectrum*.dat") 	# just gradJ
dFiles = system("ls -1 d*spectrum*.dat") 	# just d


set terminal qt size 1600, 900
set key Left
set key left bottom
subfiles=7
#set log x
set xrange
set log y
pointTypeIndex=0
set pointsize 2
if ((words(uFiles)==0) && (words(dFiles)==0)) {
	system("echo 'no uFiles or dFiles found'")
	system("sleep 1")
} else if (words(dFiles)==0) {
	system("echo 'no dFiles found'")
	FILES=uFiles
	set title system("pwd")." u";
	plot for [i=1:words(FILES)] word(FILES,i) u 1:2 title word(FILES,i) with linespoints ls i
} else if (words(uFiles)==0) {
	system("echo 'no uFiles found'")
	FILES=dFiles
	set title system("pwd")." d";
	plot for [i=1:words(FILES)] word(FILES,i) u 1:2 title word(FILES,i) with linespoints ls i
} else {
	set multiplot layout 1, 2 ;
	set title system("pwd")." u";
	FILES=uFiles
	plot for [i=1:words(FILES)] word(FILES,i) u 1:2 title word(FILES,i) with linespoints ls i
	set title system("pwd")." d";
	FILES=dFiles
	plot for [i=1:words(FILES)] word(FILES,i) u 1:2 title word(FILES,i) with linespoints ls i
	unset multiplot
}

stats word(FILES,1) using 1 nooutput
maxN=STATS_records
dealiasingCut=(2.0/3.0*pi*maxN)
set arrow from dealiasingCut, graph 0 to dealiasingCut, graph 1 nohead

pause -1
