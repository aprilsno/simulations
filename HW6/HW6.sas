*********************************************************************
*  Assignment:    BIOS668 HW6  
*  Name:          Sara O'Brien
*  Date:          3/29/23 
********************************************************************;

* Deng et al.;

* Sample size calculation; *Test of superiority;
proc power;
	twosamplemeans test=diff
	nulldiff = 0
	meandiff = 1.25 
	stddev = 1.88
	ntotal = .
	alpha = 0.05
	power = 0.9;
run;

* Power calculation;
proc power;
	twosamplemeans test=diff
	nulldiff = 0
	meandiff = 1.25 
	stddev = 1.8
	ntotal = 88
	alpha = 0.05
	power = .;
run;

* Simulation;
%macro diff(meandiff=, stddev=, alpha=, sampsize=, n_sim=);
	
data sim;
	call streaminit(730317945);
		
	do n_sim=1 to &n_sim;
		do n=1 to floor(&sampsize/2);
			* Treatment A;
			trt='A';
			y=rand('normal', 0, &stddev);
			output;
			* Treatment B;
			trt='B';
			y=rand('normal', &meandiff, &stddev);
			output;
		end;		
	end;
run;

ods select none;
	
proc ttest data = sim;
	ods output ttests=sup_rejection;
	by n_sim;
	class trt;
	var y;
run;
	
data tests_rate;
	set sup_rejection;
	where method='Pooled';
	i_sig=probt<&alpha;
run;
	
proc freq data = tests_rate;
	ods select binomial;
	tables i_sig / binomial(level='1');
run;
	
%mend;

* Check power;
%diff(meandiff=1.25,stddev=1.88,alpha=0.05,sampsize=88,n_sim=1000);

* Check type 1 error;
%diff(meandiff=0,stddev=1.88,alpha=0.05,sampsize=88,n_sim=1000);

/*The true acupuncture group improved 5.50 (SE, +/- 1.48) points 
on the Functional Assessment of Chronic Illness Therapy-Fatigue 
Subscale (FACIT-F), whereas the sham acupuncture group improved by 
3.73 (SE +/- 1.92) points. This difference was not statistically 
significant (p = .37).*/

* Turkova et al.;

* Sample size calculation; * Test of non-inferiority;
* Margin: 6 percentage points;

* Use two sample noninferiority;

* Std calculation; *Test of noninferiority;
proc power;
	twosamplemeans test=diff
	nulldiff = 0.06
	meandiff = .08 
	stddev = .
	sides = U
	ntotal = 1200
	alpha = 0.05
	power = 0.9;
run;

* Sample size calculation;
proc power;
	twosamplemeans test=diff
	nulldiff = 0.06
	meandiff = .08 
	stddev = .107

	ntotal = .
	alpha = 0.05
	power = 0.9;
run;

* Power calculation;
proc power;
	twosamplemeans test=diff
	nulldiff = 0.06
	meandiff = -0.4 
	stddev = 2.46
	ntotal = 1200
	alpha = 0.05
	power = .;
run;

* Simulation;
%macro diff(meandiff=, stddev=, alpha=, sampsize=, n_sim=);
	
data sim;
	call streaminit(730317945);
		
	do n_sim=1 to &n_sim;
		do n=1 to floor(&sampsize/2);
			* Treatment A;
			trt='A';
			y=rand('normal', 0, &stddev);
			output;
			* Treatment B;
			trt='B';
			y=rand('normal', &meandiff, &stddev);
			output;
		end;		
	end;
run;

ods select none;
	
proc ttest data = sim;
	ods output ttests=sup_rejection;
	by n_sim;
	class trt;
	var y;
run;
	
data tests_rate;
	set sup_rejection;
	where method='Pooled';
	i_sig=probt<&alpha;
run;
	
proc freq data = tests_rate;
	ods select binomial;
	tables i_sig / binomial(level='1');
run;
	
%mend;

* Check power;
%diff(meandiff=.08,stddev=.107,alpha=0.05,sampsize=1200,n_sim=1000);

* Check type 1 error;
%diff(meandiff=0,stddev=.107,alpha=0.05,sampsize=1200,n_sim=1000);