*********************************************************************
*  Assignment:    BIOS668 HW5  
*  Name:          Sara O'Brien
*  Date:          3/21/23 
********************************************************************;

/* Q1.1 (test for superiority)
	Design: parallel-group, randomized (1:1), controlled study
	Meandiff = 1.8
	Power = 0.80
	Alpha level = 0.05
	Std = 7.8 days for the two treatment strategies */

proc power;
	twosamplemeans test=diff
	nulldiff = 0
	meandiff = 1.8 
	stddev = 7.8
	ntotal = .
	alpha = 0.05
	power = 0.8;
run;

/* Q1.2 (test for equivalence)
	Design: parallel-group, randomized (1:1), controlled study
	Power = 0.80 (under the assumption that the true difference was 0.25)
	Alpha level = 0.05
	Limits = 2.2 days
	Meandiff = 0.25
	Std = 8.0 days for the two treatment strategies */
	
proc power;
	twosamplemeans test=equiv_diff 
	lower = -2.2
	upper = 2.2
	meandiff = 0.25 
	stddev = 8.0
	ntotal = .
	alpha = 0.05
	power = 0.8;
run;

/* Q2.1 (superiority)

SUPERIORITY:
drug A: x
drug B: y
Generate x and y with means mu(x)=mu1 and mu(x)=mu2
H0: mu1 = mu2
H1: mu1 ^= mu2, mu1-mu2=1.8
n total (_ x and _ y)

Use proc ttest to test if x and y are significantly different (test mu1=mu2),
if rejected count 1 if not count 0

For each iteration, run this proc ttest to get total num of rejections

Checking power: x~N(0,sigma^2) y~N(1.8,sigma^2), make sure H1 is true
Checking type 1 error: x~N(0, sigma^2), y~N(0, sigma^2), total number of rejections/1000 (make sure H0 is true) */

%macro superiority(meandiff=, stddev=, alpha=, sampsize=, n_sim=);
	
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
%superiority(meandiff=1.8,stddev=7.8,alpha=0.05,sampsize=592,n_sim=1000);

* Check type 1 error;
%superiority(meandiff=0,stddev=7.8,alpha=0.05,sampsize=592,n_sim=1000);

/* Q2.2 (equivalence) */

%macro equiv(meandiff=, stddev=, alpha=, sampsize=, n_sim=, l_equiv=, u_equiv=);
	
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
	
proc ttest data = sim tost(&l_equiv, &u_equiv);
	ods output EquivTests=equiv_rejection;
	by n_sim;
	class trt;
	var y;
run;
	
data tests_rate;
	set equiv_rejection;
	where method='Pooled' and test='Overall';
	i_sig=probt<&alpha;
run;
	
proc freq data = tests_rate;
	ods select binomial;
	tables i_sig / binomial(level='1');
run;
	
%mend;

* Check power;
%equiv(meandiff=0.25,stddev=8.0,alpha=0.05,sampsize=478,n_sim=1000,l_equiv=-2.2,u_equiv=2.2);

* Check type 1 error;
%equiv(meandiff=2.2,stddev=8.0,alpha=0.05,sampsize=478,n_sim=1000,l_equiv=-2.2,u_equiv=2.2);

/* Q3.1 (Urn Model): Simulate 100 trials each with two treatemnt arms, and total sample size 
60, and use parameter k=5. */

%macro urn(k=,sampsize=,n_sim=);
	
data sim;
	retain a 1;
	retain b 1;
	call streaminit(730317945);
	
	do n_sim=1 to &n_sim;
		do n=1 to &sampsize;
			u=rand('uniform');
			if u < (a/(a+b)) then do;
				trt='A';
				b=b+&k;
			end;
		
			else if u >= (a/(a+b)) then do;
				trt='B';
				a=a+&k;
			end;
			output;
		end;
	end;
run;

data sim_trtA;
	set sim;
	where trt='A';
run;

proc sql;
	create table sum_trtA as
	select freq(trt) as freq_trtA
	from sim_trtA
	group by n_sim;
quit;

data sum_trtA;
	set sum_trtA;
	freq_trtB=&sampsize - freq_trtA;
	diff=freq_trtA - freq_trtB;
run;

proc means data=sum_trtA;
	var freq_trtA freq_trtB diff;
run;

%mend;

%urn(k=5,sampsize=60,n_sim=100);
	
/* Q3.2 (Response Adaptive Urn Model): Simulate 100 trials each with two treatemnt arms, 
and total sample size 40, and use parameter r=3, and assuming that success rate for 
treatment A is pA=0.4 and for treatment B is pB=0.6. */

%macro respadapturn(r=,pA=,pB=,sampsize=,n_sim=);

data sim;
	retain a 1;
	retain b 1;
	call streaminit(730317945);
	
	do n_sim=1 to &n_sim;
		do n=1 to &sampsize;
			u=rand('uniform');
			if  u < (a/(a+b)) then do;
				trt='A';
				y=rand('binomial',&pA,1);
				if y=1 then do;
					a=a+&r;
				end;
				if y=0 then do;
					b=b+&r;
				end;
			end;
			else if u >= (a/(a+b)) then do;
				trt='B';
				y=rand('binomial',&pB,1);
				if y=1 then do;
					b=b+&r;
				end;
				if y=0 then do;
					a=a+&r;
				end;
			end;
			output;
		end;
	end;
run;

proc sql;
	create table sum as 
	select sum(y) as sumy
	from sim
	group by n_sim;
quit;

data sim_A;
	set sim;
	where trt='A';
run;

proc sql;
	create table sumA as
	select freq(trt) as freq_trtA
	from sim_A
	group by n_sim;
quit;

data sumA;
	set sumA;
	freq_trtB = &sampsize - freq_trtA;
	diff = freq_trtA - freq_trtB;
run;

data mergesum;
	merge sumA sum;
run;

proc means data=mergesum;
	var freq_trtA freq_trtB diff sumy;
run;

%mend;

%respadapturn(r=3,pA=0.40,pB=0.60,sampsize=40,n_sim=100);