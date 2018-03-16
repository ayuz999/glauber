all:
	loop for number 
	echo >> run.csh
	run.csh contains doScanK_my.csh XX XX XX XX
	echo >> sub.con conatins run.csh
	condor_submit sub.con
