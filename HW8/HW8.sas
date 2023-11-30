

/* Compute required sample size */
data _null_;
	std1 = 0.5;
	std2 = 1.0;
	corr = 0.15;
	sd_within = sqrt(((std1**2) + (std2**2) - (2*corr*std1*std2)) / 2);
	put sd_within=;
run;

proc power;
  twosamplemeans test=diff
  alpha=0.10
  power=0.80
  meandiff=0.80
  stddev=0.8043997
  npergroup = .
  sides=2;
run;

* Simulation;
/* Define input parameters */
%let alpha = 0.10; /* Significance level */
%let power = 0.80; /* Desired power */
%let diff = 0.8; /* True mean difference */
%let sd_within = 0.7416198487; /* Within-subject standard deviation */
%let n_periods = 2; /* Number of treatment periods */
%let n_washout = 3; /* Number of washout days */
%let n_measurements = 2; /* Number of biomarker measurements per subject */
%let n_simulations = 1000; /* Number of simulations */

/* Generate data for simulation */
data sim_data;
  do subject = 1 to 24; /* Number of subjects */
    do period = 1 to 2;
      if period = 1 then treatment = 1;
      else treatment = 2;
      biomarker = rannor(1234) * &sd_within + (treatment - 1) * &diff;
      output;
    end;
  end;
run;

/* Conduct simulation */
%let n_total = %sysevalf(%sysevalf(%sysevalf(&n_periods * &n_measurements) + &n_washout) * 40); /* Total sample size */
%let n_per_group = %sysevalf(&n_total / 2); /* Sample size per treatment arm */
%let results = sim_results; /* Output dataset name */
%macro simulate;
  %do i = 1 %to &n_simulations;
    /* Randomize treatment order */
    data sim_data;
      set sim_data;
      if mod(subject, 2) = 0 then treatment = 3 - treatment;
      if period = 2 then output; /* Exclude first period */
    run;

    /* Fit crossover model and test for treatment effect */
    proc mixed data=sim_data;
      class subject treatment period;
      model biomarker = treatment period treatment*period / solution;
      repeated / subject=subject type=cs;
      runquit;
      
      /* Store results */
      ods output TestsForFixedEffects=tests;
      data _null_;
        set tests;
        if Test = "t treatment" then do;
          p_value = 2 * cdf("t", -abs(TValue), DF);
          output;
        end;
      run;
      data &results;
        set &results;
        if _n_ = 1 then do;
          power = .;
          n_total = &n_total;
          n_per_group = &n_per_group;
        end;
        if p_value <= &alpha then power = &power;
        else power = sum(power, 1 / &n_simulations);
      run;
  %end;
%mend simulate;
%simulate;

/* Calculate power and sample size */
proc sql;
  select *, input(round(quantile(power, 0.95), 0.01), 5.) as power_95 format=5.2
  from &results
  where power > 0
  order by abs(power - &power);
quit;


proc seqdesign pss stopprob errspend;
    *TwoSidedPocock: design nstages=X method=poc INFO = CUM(0.67 1.00)
    alpha = 0.05 beta = 0.20;
    
    TwoSidedOBrienFleming: 
    design nstages=2 method=obf INFO = CUM(0.67 1.00)
    alpha = 0.025 beta = 0.20;
    
    samplesize model=TWOSAMPLEFREQ(NULLPROP=0.50 PROP=0.20 TEST=PROP REF=NULLPROP weight=1);
run;

proc seqdesign pss stopprob errspend;
 OneSidedPocock: design nstages=5 method=poc ALT=UPPER INFO = CUM(0.20 .40 .60 .80 1.000) alpha = 0.025 beta = 0.20;
 OneSidedOBrienFleming: design nstages=5 method=obf ALT=UPPER INFO = CUM(0.20 .40 .60 .80 1.000) alpha = 0.025 beta = 0.20;
 OneSidedHaybittle: design nstages=5 method=HP ALT=UPPER INFO = CUM(0.20 .40 .60 .80 1.000) alpha = 0.025 beta = 0.20;
 OneSidedOBFAlpha: design nstages=5 method=ERRFUNCOBF ALT=UPPER INFO = CUM(0.20 .40 .60 .80 1.000) alpha = 0.025 beta = 0.20;
 OneSidedPocockAlpha: design nstages=5 method=ERRFUNCPOC ALT=UPPER INFO = CUM(0.20 .40 .60 .80 1.000) alpha = 0.025 beta = 0.20;
run;

