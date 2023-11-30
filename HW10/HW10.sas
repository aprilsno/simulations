libname data "/home/u49497589/BIOS668"; 

data question1;
	set data.hw10_casecont;
run;

proc freq data=question1;
	table bfany*case / norow nocol nopercent;
run;

proc logistic data=question1 descending;
	model case=bfany age kcal;
run;

proc power;
	twosamplemeans test=diff
	nulldiff = 0
	meandiff = 0.34 
	stddev = .
	ntotal = 352
	alpha = 0.05
	power = 0.9;
run;
