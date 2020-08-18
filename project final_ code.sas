/* Final Project - Klein Section 1.3 Bone Marrow Transplant Data Set */
option ls=79;
ods rtf file='final_proj_output.rtf';
data bmt;
input g $ T1 T2 d1 d2 d3 TA A TC Chr TP Plalt Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9$ Z10;

if d1=0 and d2=0 then d_all=0;
if d1=1 and d2=0 then d_all=1;
if d2=1 then d_all=2;
cards;
1  2081  2081  0  0  0    67  1   121  1    13  1  26  33  1  0  1  1    98  0  1  0
1  1602  1602  0  0  0  1602  0   139  1    18  1  21  37  1  1  0  0  1720  0  1  0
3   164   164  1  0  1   164  0   164  0   164  0  19  
3    16    16  1  0  1    16  0    16  0    16  0  27  36  
3   363   363  1  0  1   363  0   363  0    19  1  52  48  1  1  1  0  
;
run;

/*exploration*/
proc lifetest data=bmt method=km;
time t2*d2(0);
test A chr Plalt z1 z2 z3 z4 z5 z6 z7 z8 z10;
run;

proc lifetest data=bmt method=km;
time t1*d1(0);
test A chr Plalt z1 z2 z3 z4 z5 z6 z7 z8 z10;
run;

%macro km_plot(event_time, fail_indicator, pred);
proc lifetest data=bmt method=km conftype=loglog nelson plots=survival(cl)
	plots=(s,lls) graphics;
	time &event_time*&fail_indicator(0);
	strata &pred;
	symbol1 v=none color=black line=1;
	symbol2 v=none color=black line=2;
run;
%mend km_plot;
* test statment doesn't work for categorical variables with >2 levels;
%km_plot(T2, D2, g); * sig difference;
%km_plot(T2, D2, Z9);

%km_plot(T1, D1, g); * sig difference;
%km_plot(T1, D1, Z9);

/* competing risks */
DATA mort;
  SET bmt;
  event=(d1=1 and d2=0);
  type=1;
proc print data=mort;
run;

DATA relapse;
  SET bmt;
  event=(d2=1);
  type=2;

proc print data=relapse;
run;

DATA combine;
  SET relapse mort;

*proc print data=combine;
*run;

PROC LIFETEST DATA=COMBINE PLOTS=LLS;
  TIME T2*event(0);
  STRATA type;
RUN;

PROC LOGISTIC DATA=bmt;
   where d_all ne 0;
   MODEL d_all = T2 / LINK=LOGIT;
RUN;

proc phreg data=combine;
class g(ref='1');
model  T2*d_all(0) = g A Chr Plalt  Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10/ties=efron;
strata type;
run;
proc phreg data= bmt;
class g(ref='1');
model  T2*d_all(0) = g A Chr Plalt  Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10/ties=efron;
run;

proc phreg data= bmt;
class g(ref='1');
model  T2*d_all(0,1) = g A Chr Plalt  Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10/ties=efron;
run;

proc phreg data= bmt;
class g(ref='1');
model  T2*d_all(0,2) = g A Chr Plalt  Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10/ties=efron;
run;


/* model selection with time-dependent covariates */
%macro td_cox(event_time, covar_time, fail_indicator, td_indicator, td_covar);
proc phreg data=bmt;
    class g(ref='1');
    model &event_time*&fail_indicator(0)= g Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10 &td_covar/ties=efron selection=backward;
	if &covar_time >= &event_time OR &td_indicator=0 then &td_covar=0; else &td_covar=1;
run;
%mend td_cox;

%td_cox(T2, TA, D2, A, tda);
%td_cox(T2, TC, D2, Chr, tdc);
%td_cox(T2, TP, D2, Plalt, tdp);

proc phreg data=bmt;
class g(ref='1');
model T2*d2(0) = g Z8/ties=efron;
assess ph/resample;
run;
ods graphics off;

/* final model for relapse */
proc phreg data=bmt;
    class g(ref='1');
    model T2*d2(0)= g Z8/ties=efron;
	strata g;
	output out= res survival=s xbeta=risk_score;
	assess ph/resample;
run;
ods graphics off;

/* Cox-Snell residuals: check for goodness of fit */
data d;
   set res;
   e=-log(s);
run;

proc lifetest data=d plots=(ls) notable graphics;
   time e*d2(0);/* d2 is the status of survival data*/
   symbol1 v=none;
run;

/* time to death with time-dependent covariates; */
%td_cox(T1, TA, D1, A, tda);
%td_cox(T1, TC, D1, Chr, tdc);
%td_cox(T1, TP, D1, Plalt, tdp);

ods graphics off;
ods rtf close;


