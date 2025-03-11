#!/usr/bin/gnuplot
#set key Left
#set key left top
#set log x
#set log y
maxFileNumber=1
pointTypeIndex=0
set pointsize 2
FILES = system("ls -1 *iteration_info*.dat")
do for [i=1:words(FILES)] {
	set style line i lc i		# set lc = line color
}

#f(x) = (x);	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
#f(x) = (x == 0 ? NaN : (x));	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
				# for example tc undefined and corresponding Fc = 0 so ignore those
#f(x) = (x == 0 ? NaN : abs(x));	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
#g(x) = (x);			# do nothing with column
#g(x) = (x+10**(-18));		# + 10**(-18) ~ mach_eps to avoid log(0)


plot \
for [i=1:words(FILES):1] \
	word(FILES,i) u 1:4 title word(FILES,i) with linespoints ls i pointtype 2

pause -1
