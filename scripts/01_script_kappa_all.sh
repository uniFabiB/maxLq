#!/usr/bin/gnuplot
set key Left
#set key left bottom

subfiles=5
set log x
pointTypeIndex=0
set pointsize 20
set title system("pwd")
FILES = system("ls -1 constraintDirs/B*/kappaTest/kappa*.dat")
#FILES = system("ls -1 output/constraintDirs/B*/*kappa*initial*.dat")
#FILES = system("ls -1 output/constraintDirs/B*/*kappa*end*.dat")
do for [i=1:words(FILES)] {
	do for [j=1:subfiles] {
		if(i % subfiles == j % subfiles) {
			set style line i lc j
		}
	}
	if((i-1) % subfiles == 0){
		pointTypeIndex = pointTypeIndex+2	# +2 to get better differentiable markers
	}
	set style line i pt pointTypeIndex ps 2
}
plot for [i=1:words(FILES)] word(FILES,i) u 1:6 title word(FILES,i) with linespoints ls i



pause -1

### recovery ###
#colorsArray = "red blue yellow green cyan"	# has to be number of subfiles
#set style line i lc rgb word(colorsArray,j)
