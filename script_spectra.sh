#!/usr/bin/gnuplot

subfiles=6
set log x
set log y
maxFileNumber=200
pointTypeIndex=0
set pointsize 2
do for [i=1:maxFileNumber] {
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
FILES = system("ls -1 *spectrum*.dat")
plot for [i=1:words(FILES)] word(FILES,i) u 1:2 title word(FILES,i) with linespoints ls i



pause -1

### recovery ###
#colorsArray = "red blue yellow green cyan"	# has to be number of subfiles
#set style line i lc rgb word(colorsArray,j)
