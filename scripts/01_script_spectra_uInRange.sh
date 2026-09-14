#!/usr/bin/gnuplot
#FILES = system("ls -1 */constraintDirs/B0{16..30}*/uvec*spectrum*iterend*.dat")
FILES="`bash -c 'ls -1 constraintDirs/*/spectra/uvec*spectrum*B0{17..22}*iterend*.dat'`"
print FILES
set key Left
set key left bottom
#set log x
set xrange
set log y
pointTypeIndex=0
set pointsize 2
set title system("pwd") 
plot for [i=1:words(FILES)] word(FILES,i) u 1:2 title word(FILES,i) with linespoints
file1=word(FILES,1)
stats word(FILES,1) using 1 nooutput
maxN=STATS_records
dealiasingCut=(2.0/3.0*pi*maxN)
set arrow from dealiasingCut, graph 0 to dealiasingCut, graph 1 nohead

pause -1
