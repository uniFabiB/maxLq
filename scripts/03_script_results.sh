#!/usr/bin/gnuplot

showLinesForData=1	# 1=true=default, 0=false	# show line for base data
showFit=1		# 1=true=default, 0=false	# show exponential fit
showAprioriScaling=0	# 1=true, 0=false=default	# shows the expected scaling from the a priori bound
dividBySol=0		# 1=true, 0=false=default	# divide by a priori bound
saveFitLog=0		# 1=true, 0=false=default	# save fit log file
plotEnstrophy=0		# 1=true=default, 0=false	# plots the enstrophy scaling
forceEnstrPos=1		# 1=true=default, 0=false	# force all enstrophy values to be positive to get exponent
ignoreQContByFile=1	# 1=true=default, 0=false	# ignore results_qCont_B* since they are sorted by filename
subdir=""		# default "", ex "*_used/"	# subdirectory to use, has to include "/" in the end



#showNegative=1		# 1=true=default, 0=false	# depricated




set key Right
#set key left bottom
set log x
set log y

if(saveFitLog==0){
	set fit logfile '/dev/null'
}

pointTypeIndex=0
set pointsize 2
set title system("pwd")
FILES=""
if(ignoreQContByFile==1){
	lsCommand="ls -1 ".subdir."*results*.dat | grep -v 'qCont_B'"
}
else {
	lsCommand="ls -1 ".subdir."*results*.dat"
}
FILES=system(lsCommand)


print sprintf("%g files found",words(FILES))
do for [i=1:words(FILES)] {
	print sprintf("   ".word(FILES,i))
}


array qValues[words(FILES)]
array expValues[words(FILES)]
array expErrorValues[words(FILES)]
array rValues[words(FILES)]
array rErrorValues[words(FILES)]
array titleFit[words(FILES)]
array minBpos[words(FILES)]
array maxBpos[words(FILES)]
array goalExp[words(FILES)]

posValues(x) = (x > 0 ? x : 1/0);

# f(u) = d u^a => log(f(u)) = log(d u^a) = log(d) + a log(u) = c + p x => d = exp(c)

do for [i=1:words(FILES)] {
	if(strstrt(word(FILES,i), 'q4')>0 || strstrt(word(FILES,i), 'q04')>0) {
		qValues[i]=4
		print "q=4"
		goalExp[i]=12
	}
	else if(strstrt(word(FILES,i), 'q5')>0 || strstrt(word(FILES,i), 'q05')>0) {
		qValues[i]=5
		print "q=5"
		goalExp[i]=10
	}
	else if(strstrt(word(FILES,i), 'q6')>0 || strstrt(word(FILES,i), 'q06')>0) {
		qValues[i]=6
		print "q=6"
		goalExp[i]=10
	}
	else if(strstrt(word(FILES,i), 'q7')>0 || strstrt(word(FILES,i), 'q07')>0) {
		qValues[i]=7
		print "q=7"
		goalExp[i]=10.5
	}
	else if(strstrt(word(FILES,i), 'q8')>0 || strstrt(word(FILES,i), 'q08')>0) {
		qValues[i]=8
		print "q=8"
		goalExp[i]=11.2
	}
	else if(strstrt(word(FILES,i), 'q9')>0 || strstrt(word(FILES,i), 'q09')>0) {
		qValues[i]=9
		print "q=9"
		goalExp[i]=12
	}
	else {
		qValues[i]=0
		goalExp[i]=0
	}
	
	
	stats word(FILES,i) using ($3 > 0 ? 1 : 0) nooutput
	positive_count = STATS_sum
	print "number positive: ", positive_count
	print word(FILES,i)
	if(positive_count > 2) {
		fitFunc(x) = x*p+c;
		fit fitFunc(x) word(FILES,i) u (log($2)):(log($3)) via p, c;
		expValues[i]=p;
		expErrorValues[i]=p_err;
		rValues[i]=exp(c);
		rErrorValues[i]=exp(c_err);
		titleFit[i]=word(FILES,i).sprintf(' fit: {%.2f}*B^{%.5f}', rValues[i],p);
		stats word(FILES,i) using 2:(posValues($3)) nooutput;
		minBpos[i]=STATS_min_x;
		maxBpos[i]=STATS_max_x;
		print sprintf("min %g, max %g",minBpos[i], maxBpos[i])
		print word(FILES,i).sprintf("\n\n\nq value %i",qValues[i])
		print word(FILES,i).sprintf("\texponent = %.5f +- %.5f",p,p_err)
		print word(FILES,i).sprintf("\tconstant = %.5e */ %.5e\n\n\n",rValues[i],rErrorValues[i])
	}
	else {
		stats word(FILES,i) using 2 nooutput;
		minBpos[i]=STATS_min;
		maxBpos[i]=STATS_max;
		rValues[i]=-1.0;
		expValues[i]=0.0;
		titleFit[i]='to few positive values to fit '.word(FILES,i);
	}
}

style(x) = (x > 0 ? 2 : 10);

#if(showNegative==0){
#	usingCommand=" u 2:3 lc i"
#}
#else {
#	usingCommand=" u 2:(abs($3)):(style($3)):1 title word(FILES,i).'   × pos, ▿ neg' pt variable lc i"
#}
if(dividBySol==0){
	usingCommand=" u 2:(abs($3)):(style($3)):1 title word(FILES,i).'   × pos, ▿ neg' pt variable lc i"	
}
else {
	usingCommand=" u 2:(abs($3)/(abs($2)**(goalExp[i]))):(style($3)):1 title word(FILES,i).'/x^{q(q-1)/(q-3)}     × pos, ▿ neg' pt variable lc i"
	showAprioriScaling=0
	showFit=0
}
if(showLinesForData==0){
	linesPointsCommand=" "
}
else {
	linesPointsCommand=" with linespoints"
}
baseCommand="plot for [i=1:words(FILES)] word(FILES,i)"
if(showFit==0){
	fitCommand=" "
}
else {
	fitCommand=", for [i=1:words(FILES)] [minBpos[i]/1.1:1.1*maxBpos[i]] + (rValues[i]*x**expValues[i]) t titleFit[i] lc i dashtype 2"
}
if(showAprioriScaling==0){
	goalCommand=" "
}
else {
	goalCommand=", for [i=1:words(FILES)] [minBpos[i]/1.1:1.1*maxBpos[i]] + (rValues[i]*x**(goalExp[i])*minBpos[i]**(expValues[i]-goalExp[i])) t word(FILES,i).' a pri: B^{'.goalExp[i].'}' lc i dashtype 3"
}
if(plotEnstrophy==0){
	enstrophyCommand=" "
}
else {

	array expValuesEnst[words(FILES)]
	array rValuesEnst[words(FILES)]
	array titleFitEnst[words(FILES)]
	array minBposEnst[words(FILES)]
	array maxBposEnst[words(FILES)]
	array goalExpEnst[words(FILES)]

	do for [i=1:words(FILES)] {
		stats word(FILES,i) using ($6 > 0 ? 1 : 0) nooutput
		positive_count = STATS_sum
		print ""
		print "-----------------------------"
		print "number positive: ", positive_count
		if(positive_count > 2) {
			print "count > 2"
			fitFunc(x) = x*p+c;
			fit fitFunc(x) word(FILES,i) u (log($5)):(log($6)) via p, c;
			expValuesEnst[i]=p;
			rValuesEnst[i]=exp(c);
			titleFitEnst[i]=word(FILES,i).sprintf(' fit: {%.2f}*enstrophy^{%.2f}', rValuesEnst[i],p);
			stats word(FILES,i) using 2:(posValues($6)) nooutput;
			minBposEnst[i]=STATS_min_x;
			maxBposEnst[i]=STATS_max_x;
		}
		else {
			if(forceEnstrPos) {
				print "forcing positive enstrophy"
				fitFunc(x) = x*p+c;
				fit fitFunc(x) word(FILES,i) u (log($5)):(log(abs($6))) via p, c;
				expValuesEnst[i]=p;
				rValuesEnst[i]=exp(c);
				titleFitEnst[i]=word(FILES,i).sprintf(' fit: -{%.2f}*enstrophy^{%.2f}', rValuesEnst[i],p);
				stats word(FILES,i) using 2:(posValues(abs($6))) nooutput;
				minBposEnst[i]=STATS_min_x;
				maxBposEnst[i]=STATS_max_x;
			}
			else {
				stats word(FILES,i) using 2 nooutput;
				minBposEnst[i]=STATS_min;
				maxBposEnst[i]=STATS_max;
				rValuesEnst[i]=-1.0;
				expValuesEnst[i]=0.0;
				titleFitEnst[i]='to few positive values to fit enstrophy for '.word(FILES,i);
			}
		}
	}
	enstrophyCommand=", for [i=1:words(FILES)] word(FILES,i) u 2:(abs($5)):(style($6)):1 title word(FILES,i).' enstrophy   × pos, ▿ neg' pt variable lc i dashtype 3 with linespoints"
	if(showFit==1){
		enstrophyCommand=enstrophyCommand.", for [i=1:words(FILES)] [minBposEnst[i]/1.1:1.1*maxBposEnst[i]] + (rValuesEnst[i]*x**expValuesEnst[i]) t titleFitEnst[i] lc i dashtype 4"
	}
}
command=baseCommand.usingCommand.linesPointsCommand.fitCommand.goalCommand.enstrophyCommand
eval command
pause -1
