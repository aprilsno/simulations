*********************************************************************
*  Assignment:    BIOS668 HW7  
*  Name:          Sara O'Brien
*  Date:          4/9/23 
********************************************************************;

proc power;
	twosamplesurvival test=logrank
	curve("Standard") = 5 : 0.75 /* Standard survival changed to 75% */
	curve("Proposed") = 5 : 0.85 /* Proposed survival still 85% */
	groupsurvival = "Standard" | "Proposed"
	accrualtime = 4
	totaltime = 7
	npergroup = .
	power = 0.8;
run;

* Original .73 hr and sample size = 1832;
data sim_data;
  call streaminit(730317945);
  do n_sim = 1 to 1000;
      do i = 1 to floor(1832/2);
      	hrgroup='Standard';
        y = rand('exponential', 1/0.73);
        c = rand('uniform', 0,7);
        t = min(y,c);
        if y <= c then censor = 1; 
        else if y > c then censor = 0;
        output;
        hrgroup='Proposed';
        y = rand('exponential', 1/0.85);
        c = rand('uniform', 0,7);
        t = min(y,c);
        if y <= c then censor = 1; 
        else if y > c then censor = 0;
        output;
      end;
  end;
run;

ODS EXCLUDE ALL;
proc lifetest data=sim_data;  
	by n_sim;
	time t*censor(0);
	strata hrgroup;
	ods output HomTests = Table_logrank;
run;
ODS EXCLUDE NONE;

data tests_rate;
	set Table_logrank;
	where Test='Log-Rank';
	i_sig=ProbChiSq<0.05;
run;

proc freq data = tests_rate;
	ods select binomial;
	tables i_sig / binomial(level='1');
run;
	
* Simulation with hr 0.75 and proc power sample size 248;
data sim_data;
  call streaminit(730317945);
  do n_sim = 1 to 1000;
      do i = 1 to floor(248);
      	trt='Standard';
        y = rand('exponential', 1/0.75);
        c = rand('uniform', 0,7);
        t = min(y,c);
        if y <= c then censor = 1; 
        else if y > c then censor = 0;
        output;
        trt='Placebo';
        y = rand('exponential', 1/0.85);
        c = rand('uniform', 0,7);
        t = min(y,c);
        if y <= c then censor = 1; 
        else if y > c then censor = 0;
        output;
      end;
  end;
run;

ODS EXCLUDE ALL;
proc lifetest data=sim_data;  
	by n_sim;
	time t*censor(0);
	strata trt;
	ods output HomTests = Table_logrank;
run;
ODS EXCLUDE NONE;

data tests_rate;
	set Table_logrank;
	where Test='Log-Rank';
	i_sig=ProbChiSq<0.05;
run;

proc freq data = tests_rate;
	ods select binomial;
	tables i_sig / binomial(level='0');
run;
	

