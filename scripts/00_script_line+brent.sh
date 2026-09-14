#!/usr/bin/gnuplot
#set key Left
set key left bottom
set log x
forceLog = 1			# 0: nothing		1: forces log y by setting y -> |y|
plot0 = 1			# 0: nothing		1: adds 10^-20 to 0 x values to plot these in log scale 
plotTitles = 1			# 0: don't plot		1: plot titles



if(forceLog) {
	set log y
}
FILES_all_J = system("ls -1 *all*J*info*.dat")
FILES_brent = system("ls -1 *brent*info*.dat")
FILES_brak = system("ls -1 *lineMin*info*.dat")
do for [i=1:words(FILES_brent)] {
	set style line i lc i		# set lc = line color
}

#f(x) = (x);	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
f(x) = (x == 0 ? NaN : (x));	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
				# for example tc undefined and corresponding Fc = 0 so ignore those
if(forceLog) {
	f(x) = (x == 0 ? NaN : abs(x));	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
}
else {
	f(x) = (x == 0 ? NaN : x);	# ($6 == 0 ? NaN : abs($6)) checks for 0 values and ignores them by setting it NaN
}

if(plot0) {
	g(x) = (x+10**(-20));		# + 10**(-20) ~ mach_eps to avoid log(0)
}
else {
	g(x) = (x);			# do nothing with column
}


set key noautotitle		# no titles in plot

if (plotTitles) {		# switch titles on/off
	plot \
	for [i=1:words(FILES_brent):1] \
		word(FILES_brent,i) u (g($2)):(f($3)) title word(FILES_brent,i)." brent" with points ls i pointtype 2 pointsize 2, \
	for [i=1:words(FILES_brak):1] \
		word(FILES_brak,i) u (g($2)):(f($5)) title word(FILES_brak,i)." brak fa" with points ls i pointtype 4 pointsize 2, \
	for [i=1:words(FILES_brak):1] \
		word(FILES_brak,i) u (g($3)):(f($6)) title word(FILES_brak,i)." brak fb" with points ls i pointtype 6 pointsize 2, \
	for [i=1:words(FILES_brak):1] \
		word(FILES_brak,i) u (g($4)):(f($7)) title word(FILES_brak,i)." brak fc" with points ls i pointtype 8 pointsize 2, \
	for [i=1:words(FILES_all_J):1] \
		word(FILES_all_J,i) u (g($3)):(f($6)) title word(FILES_all_J,i)." all" with points ls i pointtype 4 pointsize 0, \
}
else {
	plot \
	for [i=1:words(FILES_brent):1] \
		word(FILES_brent,i) u (g($2)):(f($3)) with points ls i pointtype 2 pointsize 2, \
	for [i=1:words(FILES_brak):1] \
		word(FILES_brak,i) u (g($2)):(f($5)) with points ls i pointtype 4 pointsize 2, \
	for [i=1:words(FILES_brak):1] \
		word(FILES_brak,i) u (g($3)):(f($6)) with points ls i pointtype 6 pointsize 2, \
	for [i=1:words(FILES_brak):1] \
		word(FILES_brak,i) u (g($4)):(f($7)) with points ls i pointtype 8 pointsize 2, \
	for [i=1:words(FILES_all_J):1] \
		word(FILES_all_J,i) u (g($3)):(f($6)) with points ls i pointtype 4 pointsize 0, \
}

pause -1
