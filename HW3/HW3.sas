*********************************************************************
*  Assignment:    BIOS668 HW3  
*  Name:          Sara O'Brien
*  Date:          2/14/23 
********************************************************************;

ODS PDF FILE="/home/u49497589/BIOS668/OBrien_HW3.pdf" STYLE=JOURNAL;

********************************************
Question 1: Code for simulating 3+3 design
*******************************************;

* Define update function;
proc iml;
  start updata(data_m, l_dose, n_pt, n_dlt); 
 
	if ( (l_dose < 1) | (l_dose >  nrow(data_m) ) ) then do;    
        msg = "dose do not treat";   
       	return(msg) ; 
      	end; 
       
    if (n_dlt > n_pt) then do;  
        msg = "ndlt > npt is impossible";   
        return(msg) ; 
      	end; 
 
    if ( (n_dlt < 0) | (n_pt < 0) ) then do;   
        msg = "ndlt & npt must be positive";   
        return(msg) ; 
      	end; 
        
    if (n_pt = 0) then  
       return(data_m); 
 
      	idx = l_dose;    
      	data_m[idx, 2] = data_m[idx, 2] + n_pt;    
      	data_m[idx, 3] = data_m[idx, 3] + n_dlt;   
    	return(data_m); 
	finish; 
 
* Define 3+3 decision;
  start TriPTri(data_m, l_dose); 
    idx =  l_dose; 
    ndlt = data_m[idx, 3]; 
    npt = data_m[idx, 2]; 
   
    mtd = .;  
 
    if (npt = 3) then do;  
      if (ndlt = 0) then nextdose = l_dose + 1; 
      if (ndlt = 1) then nextdose = l_dose; 
      if (ndlt > 1) then nextdose = .; 
    end;  
 
    if (npt = 6) then do;  
      if (ndlt = 1) then nextdose = l_dose + 1; 
      if (ndlt > 1) then nextdose = .; 
    end;  
    
      if ( ( (nextdose < 1) | (nextdose > nrow(data_m)) ) | ndlt > 1) then do;   
        nextdose = .;  
        if (l_dose = 1) then mtd = .; 
      if (l_dose > 1 & ndlt > 1) then mtd = l_dose - 1 ; 
      if (l_dose > 1 & ndlt <= 1) then mtd = l_dose; 
      end;    
     
      res = j(1, 2, 0); 
      res[1] = nextdose;
      res[2] = mtd; 
      return(res);  
  	finish; 
  
* Simulate 1 trial of 3+3; 
* Data structure: empty #pt and # DLT;

  start SimulOneTriPTri(TrueRates) ;  
 
  DLevel = NROW(TrueRates); 
 
   Data_mat = j(DLevel , 3, 0); /* dose level, npt, ndlt */ 
   do Dose = 1 to DLevel; 
   		Data_mat[Dose, 1] = Dose;    
   end; 
 
    nextdose = 1; 
   do until( (nextdose < 1) | (nextdose > DLevel) ); 
        lastdose = nextdose; 
    	call randgen(x, 'BINOM', TrueRates[lastdose], 3);   
    	ndlt = x; 
 
    	ttmp = j(1,2,.); 
    	ttmp[1] = lastdose; 
    	ttmp[2] = ndlt; 
 
      tmp = updata(data_mat, lastdose, 3, ndlt); 
    	if ( (nrow(tmp) = DLevel) & (ncol(tmp) = 3) ) THEN data_mat = tmp; 
 
    	myres = TriPTri(data_mat, lastdose); 
      nextdose = myres[1]; 
    end; 
 
  	myres = TriPTri(data_mat, lastdose); 
    mtd = myres[2]; 
 
    data_res = j(DLevel , 3+2, 0); 
    data_res[, 1:3] = data_mat; 
    data_res[, 4] = mtd; 
    data_res[, 5] = lastdose; 
 
    return(data_res);  
  finish; 
 
* Simulate multiple trials of 3+3; 
  start SimulTriPTri(TrueRates, Nrep, Seed) ;  
    MyTrueRates = TrueRates;   
    call randseed(seed);
  	MyDLevel = NROW(MyTrueRates); 
   
    Trip_Data_simul = j(Nrep*MyDLevel, 6, . ); 
    DO rep = 1 TO Nrep; 
      my_data_res =  SimulOneTriPTri(MyTrueRates); 
    	tmp_data = insert(my_data_res, j(MyDLevel,1 ,rep), 0, 1); 
      idxn = (rep-1)*MyDLevel  ; 
      Trip_Data_simul[(idxn+1):(idxn+MyDLevel), ] = tmp_data; 
    end;   
    return(Trip_Data_simul);  
  finish;  
 
* Call the main function for simulating 3+3 trials;
TrueRates = {.08, .12, .15, .21, .35, .40}; 
Nrep = 2000;  
Seed = 730317945;  
MySimul3p3 = SimulTriPTri(TrueRates, Nrep, Seed) ;  
 
* Export the data matrix to a SAS dataset; 
VarName = {"Rep" "DoseLevel" "Num_pt" "Num_DLT" "MTD" "LastDose" }; 
create TriPTriData from MySimul3p3[colname=VarName]; 
append from MySimul3p3; 
close TriPTriData; 
 
quit; 
 
proc sort data = Triptridata; 
  by DoseLevel; 
run; 

proc univariate data = Triptridata; 
  title1 'Mean number of patients treated and DLTs';
  title2 'at each dose level';
  var Num_pt Num_DLT MTD; 
  by DoseLevel; 
run; 
 
proc freq data = Triptridata (where = (DoseLevel = 1)); 
  title 'Distributions of estimated MTD';
  table MTD; 
  *by DoseLevel; 
run; 

*******************************************************************
Question 2.1: find the single-stage design with minimal sample size 
******************************************************************;

/* Single stage */ 
%macro one_stage( 
p0= , /* Unacceptable response */ 
p1= , /* Acceptable response */ 
alpha= , /* Type I Error (probability of accepting a poor drug) */ 
beta= , /* 1-Power (probability of rejecting a good drug) */ 
usern=  /* User-specified maximum sample size */ 
); 
 
DATA one_stage_output;  
  RETAIN low_n 0 p0 0 p1 0 alpha 0 beta 0 stat 0; 
 
  p0=&p0; p1=&p1; alpha=&alpha; beta=&beta; 
  low_n = INT( 0.25*(p0+p1)*(2-p0-p1)*( ( ( quantile('NORMAL', 1-alpha)+ quantile('NORMAL', 1-beta))/
  (p1-p0))**2)); 
   
  DO n = low_n to &usern;  /* Total sample size */ 
      DO r = 0 to n; /* Stage 1 cut-off */ 
        term1_p0 = cdf('BINOMIAL', r, &p0, n); 
        term1_p1 = cdf('BINOMIAL', r, &p1, n); 
        if (term1_p1 <=  &beta and 1-term1_p0<=&alpha) then output; 
      END; 
  END;   
RUN; 
%mend one_stage; 
 
%one_stage(p0=0.30,p1=0.60,alpha=0.10,beta=0.20,usern=80); 
 
DATA result_singlestage; 
  SET one_stage_output; 
  real_alpha = 1-(term1_p0); 
  real_beta = (term1_p1);  
  KEEP alpha beta p0 p1 n r real_alpha real_beta; 
RUN; 

PROC SORT DATA = result_singlestage; 
  BY n; 
RUN; 

proc print data = result_singlestage (obs=1);
	title 'Single-stage design with minimal sample size';
run;

*******************************************************
Question 2.1: find the Optimal Simon’s 2-stage design
******************************************************;

/* Simon's two-stage */
%macro simon_twostage( 
p0= , /* Unacceptable response */ 
p1= , /* Acceptable response */ 
alpha= , /* Type I Error (probability of accepting a poor drug) */ 
beta= , /* 1-Power (probability of rejecting a good drug) */ 
usern=  /* User-specified maximum sample size */ 
); 
 
DATA simon_twostage; 
  RETAIN low_n 0 p0 0 p1 0 alpha 0 beta 0 stat 0; 
 
  p0=&p0; p1=&p1; alpha=&alpha; beta=&beta; 
  low_n = INT( 0.25*(p0+p1)*(2-p0-p1)*( ( ( quantile('NORMAL', 1-alpha)+ quantile('NORMAL', 1-beta))/(p1-p0) 
  )**2)); 
   
  DO n = low_n to &usern;  /* Total sample size */ 
    DO n1 = 1 to n-1;  /* Stage 1 Sample size */ 
      DO r1 = 0 to n1; /* Stage 1 cut-off */ 
        term1_p0 = cdf('BINOMIAL', r1, &p0, n1); 
        term1_p1 = cdf('BINOMIAL', r1, &p1, n1); 
        if term1_p1=<&beta then DO; /*remove solution sets that do not meet the beta requirement*/ 
           stat = 0; 
           DO r = n to r1 BY - 1 WHILE (stat=0); 
              term2_p0=0; *initialize the summation terms for alpha & beta calculations; 
              term2_p1=0; 
              do x=r1+1 to min(r, n1); 
                dum0=pdf('BINOMIAL', x, &p0, n1)*cdf('BINOMIAL', r-x, &p0, n-n1); 
                dum1=pdf('BINOMIAL', x, &p1, n1)*cdf('BINOMIAL', r-x, &p1, n-n1); 
                term2_p0= term2_p0 + dum0;   
                term2_p1= term2_p1 + dum1; 
                end; 
              if 1-(term1_p0+term2_p0)=<&alpha and (term1_p1+term2_p1)=<&beta then DO; 
                              
output; stat = 1; END; 
                  END;     
                END;  
      END; 
    END; 
  END;   
RUN; 
%mend simon_twostage; 
 
%simon_twostage( 
p0=0.30,  /* Unacceptable response */ 
p1=0.60,  /* Acceptable response */ 
alpha=0.10,  /* Type I Error (probability of accepting a poor drug) */ 
beta=0.20,  /* 1-Power (probability of rejecting a good drug) */ 
usern=80  /* User-specified maximum sample size */ 
); 
 
DATA result_twostage; 
  SET simon_twostage; 
  real_alpha = 1-(term1_p0+term2_p0); 
  real_beta = (term1_p1+term2_p1);  
  PET = term1_p0; 
  EN = n1+(1-PET)*(n-n1); 
  KEEP alpha beta p0 p1 n1 r1 n r real_alpha real_beta PET EN; 
RUN; 
 
PROC SORT DATA = result_twostage; 
  BY EN; 
RUN; 

proc print data = result_twostage (obs=1);
	title 'Optimal Simon’s 2-stage design';
run;

********************************************************
Question 2.2: Simulate single-stage studies and use the 
code to verify the design in (2.1). 
*******************************************************;

data sim1;
	call streaminit(730317945);
	
	type1 = 0;
	power = 0;
	
	do rep=1 to 5000;
		null = rand('Binom', .3, 14);
		if null>6 then type1=type1+(1/5000);
		else type1=type1;
		alt = rand('Binom', .6, 14);
		if alt>6 then power=power+(1/5000);
		else power=power;
	end;
run;

proc print data=sim1;
	title 'Type I Error and Power of 1-stage simulation';
	var type1 power;
run;

********************************************************
Question 2.3: Simulate Simon's 2-stage studies and use 
the code to verify the design in (2.1). 
*******************************************************;
	
data sim2;
	call streaminit(730317945);
	
	type1 = 0;
	power = 0;
	
	do rep=1 to 5000;
		null = rand('Binom', .3, 5);
		alt = rand('Binom', .6, 5);
		
		if null <= 1 then null=null;
			else null = null+rand('Binom', .3, 9);
		if alt <= 1 then alt=alt;
			else alt = alt+rand('Binom', .6, 9);
		if null>6 then type1=type1+(1/5000);
		if alt>6 then power=power+(1/5000);
	end;
run;

proc print data=sim2;
	title 'Type I Error and Power of 2-stage simulation';
	var type1 power;
run;

ODS PDF CLOSE;