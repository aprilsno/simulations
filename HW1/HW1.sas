*********************************************************************
*  Assignment:    HW1  
*  Name:          Sara O'Brien
*  Date:          1/19/23 
********************************************************************;

*********************************************************************
Q1:
Please write SAS code to generate 5 random numbers from the following
distribution (use your PID as seed)
(1) Binomial: with p=.25 n=25, where p is probability of success and n 
is the number of trials.
(2) Negative binomial: with p=.25 k=25, where p is probability of success 
and k is the number of successes.
(3) Exponential distribution: with rate parameter lamba=2.0 (i.e., 
mean = 0.5).
********************************************************************;

* 5 random variables with binomial distribution;
data randbinomial;
call streaminit(730317945);       /* set random number seed */
do i = 1 to 5;
   u = rand("Binomial",.25,25);     /* u ~ Binom(.25,25) */
   drop i;
   output;
end;
run;

* Print 5 vars;
proc print data=randbinomial noobs label;
	title1 '5 random numbers from';
	title2 'binomial distribution';
	label u='Random number';
run;

* 5 random variables with neg binomial distribution;
data randnbinomial;
call streaminit(730317945);       /* set random number seed */
do i = 1 to 5;
   u = rand("NEGBinomial",.25,25);     /* u ~ NBinom(.25,25) */
   drop i;
   output;
end;
run;

* Print 5 vars;
proc print data=randnbinomial noobs label;
	title1 '5 random numbers from';
	title2 'negative binomial distribution';
	label u='Random number';
run;

* 5 random variables with exponential distribution;
data randexp;
call streaminit(730317945);       /* set random number seed */
do i = 1 to 5;
   u = rand("EXPONENTIAL",2.0);     /* u ~ Exponential(2.0) */
   output;
   drop i;
end;
run;

* Print 5 vars;
proc print data=randexp noobs label;
	title1 '5 random numbers from';
	title2 'exponential distribution';
	label u='Random number';
run;

*********************************************************************
Q2:
Please write SAS/R code to conduct simulation studies as specified below. 
In this simulation, you simulate data from a linear regression model as we 
showed in the lecture notes. In one linear regression analysis, only Trt is 
included as the covariate, and in another analysis, other covariates are 
included in the linear regression model. In both analyses, the main interest 
is whether Trt has statistically significant (p<0.05) effect on the outcome.
Please try both options with parameters specified as follows:

Parameters 	Set 1 		Set 2
ssize 		300 		300
reps 		500 		500
beta0 		2.3 		2.3
beta1 		0.0 		1.5
beta2 		1.0 		1.0
beta3 		-3.5 		-3.5
beta4 		-1.5 		-1.5
Sd_e 		3.5 		3.5
seed0 		Your PID 	Your PID
********************************************************************;

options NOERRORABEND;

* Simulate data from linear model set 1;
data sim1;
	call streaminit(730317945);
	
	do rep=1 to 500;
		do i=1 to 300;
			Trt = rand('Binom', .5, 1);
			gender = rand('Binom', .3, 1);
			age = 20 + 50*rand('Uniform');
			BP = 75 + 20*rand('Uniform');
			epsilon = 3.5*rand('Normal');
			y = 2.3 + 0.0*Trt + 1.0*gender - 3.5*age - 1.5*BP + epsilon;
			output;
		end;
	end;
run;

* Create output dataset of reps with significant p values for F statistic;
proc glm data = sim1;
	class Trt;
	by rep;
	model y = Trt / solution;
	ods output OverallANOVA = sim1trt(keep=ProbF where=(0 <=ProbF< 0.05));
run;
* 22 reps with p-value<.05;
	
* Create output dataset of reps with significant p values for F statistic;
proc glm data = sim1;
	class Trt;
	by rep;
	model y = Trt gender age bp / solution;
	ods output OverallANOVA = sim1all(keep=ProbF where=(0 <=ProbF< 0.05));
run;
* 500 p-value<.05;

* Simulate data from linear model set 2;
data sim2;
	call streaminit(730317945);
	
	do rep=1 to 500;
		do i=1 to 300;
			Trt = rand('Binom', .5, 1);
			gender = rand('Binom', .3, 1);
			age = 20 + 50*rand('Uniform');
			BP = 75 + 20*rand('Uniform');
			epsilon = 3.5*rand('Normal');
			y = 2.3 + 1.5*Trt + 1.0*gender - 3.5*age - 1.5*BP + epsilon;
			output;
		end;
	end;
run;

* Create output dataset of reps with significant p values for F statistic;
proc glm data = sim2;
	class Trt;
	by rep;
	model y = Trt / solution;
	ods output OverallANOVA = sim2trt(keep=ProbF where=(0 <=ProbF< 0.05));
run;
* 28 p-value<.05;

* Create output dataset of reps with significant p values for F statistic;
proc glm data = sim2;
	class Trt;
	by rep;
	model y = Trt gender age bp / solution;
	ods output OverallANOVA = sim2all(keep=ProbF where=(0 <=ProbF< 0.05));
run;
* 500 p-value<.05;

*********************************************************************
Q3 (2'): Please write SAS/R code to simulate the gambling problem, as 
specified below.

A gambler plans to play in a gambling game.

In this game, each time, if he wins, he gets a reward of the same value 
of bet: if he loses, he loses his bet. For example, if he puts $100 bet 
to play, he gets $200 ($100 bet + $100 reward) back if he wins after a 
game, but he loses $100 (bet) if he loses this game.

Each time, his chance of winning is p (e.g., p=55%).

He has A0 = $2000 in the beginning.

He plans to plan this game for T times (e.g., T=15).

What does his final asset look like if he puts X% (e.g., X% = 50%) of what 
he has as his bet each time?

Please write SAS/R code to simulate (use your PID as the seed) the above 
games, and replicate the same games for 5000 times. After each replicate, 
please record the final asset of this gambler, and then (1) draw a histogram 
of distribution of the asset, and (2) calculate the mean, median, 25%, 75% 
quantiles of the final assets.
********************************************************************;

* Run a simulation of 5000 reps of the gambling game;
data gamble;
	call streaminit(730317945);
	retain asset;
	do rep=1 to 5000;
	asset = 2000;
		do i=1 to 15;
			win = rand('bernoulli', 0.55);
			bet = asset*0.5;
			if win = 0 then asset = asset - bet;
			else asset = asset + bet;		
		end;
		output;
	end;
run;

* Produce a histogram and univariate report of final assets;
proc univariate data=gamble;
	var asset;
	histogram asset;
run;