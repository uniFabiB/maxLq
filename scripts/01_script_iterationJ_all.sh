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
FILES = system("ls -1 */constraintDirs/B*/iteration-info*.dat")


numberFilesPerColor = "0"		# number of files that should get the same color, if first <= 0, then all differently colored



do for [i=1:words(FILES)] {
	filesStart = 0
	filesEnd = 0
	if(word(numberFilesPerColor,1)<1){
		set style line i lc i		
	}
	else {
		do for [j=1:words(numberFilesPerColor)] {
			numberOfFiles = word(numberFilesPerColor,j)
			filesEnd = filesEnd + numberOfFiles
			if(i>filesStart && i<=filesEnd){
				set style line i lc j
			}
			filesStart = filesStart + numberOfFiles
		}
	}
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

