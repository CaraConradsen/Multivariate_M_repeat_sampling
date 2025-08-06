
proc sort data = mrwings;
by Trait Gen;
run;
quit;

proc standard data=mrwings out=stdmrwings std=1 mean=0;
by Trait;
var Score;
run;
quit;




ods output covparms=Covparms1_1;
ods output FitStatistics=FitStat1_1;
ods output ConvergenceStatus=Converge1_1;
ods output asycov=Asycov1_1;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (1) and Gen in (1);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms1_1
OUTFILE= ".\un_sas_output\Covparms1_1.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat1_1
OUTFILE= ".\un_sas_output\FitStat1_1.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge1_1
OUTFILE= ".\un_sas_output\Converge1_1.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov1_1
OUTFILE= ".\un_sas_output\Asycov1_1.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms1_2;
ods output FitStatistics=FitStat1_2;
ods output ConvergenceStatus=Converge1_2;
ods output asycov=Asycov1_2;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (1) and Gen in (2);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms1_2
OUTFILE= ".\un_sas_output\Covparms1_2.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat1_2
OUTFILE= ".\un_sas_output\FitStat1_2.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge1_2
OUTFILE= ".\un_sas_output\Converge1_2.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov1_2
OUTFILE= ".\un_sas_output\Asycov1_2.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms1_3;
ods output FitStatistics=FitStat1_3;
ods output ConvergenceStatus=Converge1_3;
ods output asycov=Asycov1_3;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (1) and Gen in (3);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms1_3
OUTFILE= ".\un_sas_output\Covparms1_3.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat1_3
OUTFILE= ".\un_sas_output\FitStat1_3.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge1_3
OUTFILE= ".\un_sas_output\Converge1_3.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov1_3
OUTFILE= ".\un_sas_output\Asycov1_3.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms1_4;
ods output FitStatistics=FitStat1_4;
ods output ConvergenceStatus=Converge1_4;
ods output asycov=Asycov1_4;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (1) and Gen in (4);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms1_4
OUTFILE= ".\un_sas_output\Covparms1_4.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat1_4
OUTFILE= ".\un_sas_output\FitStat1_4.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge1_4
OUTFILE= ".\un_sas_output\Converge1_4.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov1_4
OUTFILE= ".\un_sas_output\Asycov1_4.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms1_5;
ods output FitStatistics=FitStat1_5;
ods output ConvergenceStatus=Converge1_5;
ods output asycov=Asycov1_5;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (1) and Gen in (5);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms1_5
OUTFILE= ".\un_sas_output\Covparms1_5.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat1_5
OUTFILE= ".\un_sas_output\FitStat1_5.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge1_5
OUTFILE= ".\un_sas_output\Converge1_5.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov1_5
OUTFILE= ".\un_sas_output\Asycov1_5.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms1_6;
ods output FitStatistics=FitStat1_6;
ods output ConvergenceStatus=Converge1_6;
ods output asycov=Asycov1_6;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (1) and Gen in (6);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms1_6
OUTFILE= ".\un_sas_output\Covparms1_6.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat1_6
OUTFILE= ".\un_sas_output\FitStat1_6.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge1_6
OUTFILE= ".\un_sas_output\Converge1_6.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov1_6
OUTFILE= ".\un_sas_output\Asycov1_6.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms3_1;
ods output FitStatistics=FitStat3_1;
ods output ConvergenceStatus=Converge3_1;
ods output asycov=Asycov3_1;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (3) and Gen in (1);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms3_1
OUTFILE= ".\un_sas_output\Covparms3_1.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat3_1
OUTFILE= ".\un_sas_output\FitStat3_1.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge3_1
OUTFILE= ".\un_sas_output\Converge3_1.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov3_1
OUTFILE= ".\un_sas_output\Asycov3_1.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms3_2;
ods output FitStatistics=FitStat3_2;
ods output ConvergenceStatus=Converge3_2;
ods output asycov=Asycov3_2;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (3) and Gen in (2);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms3_2
OUTFILE= ".\un_sas_output\Covparms3_2.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat3_2
OUTFILE= ".\un_sas_output\FitStat3_2.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge3_2
OUTFILE= ".\un_sas_output\Converge3_2.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov3_2
OUTFILE= ".\un_sas_output\Asycov3_2.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms3_3;
ods output FitStatistics=FitStat3_3;
ods output ConvergenceStatus=Converge3_3;
ods output asycov=Asycov3_3;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (3) and Gen in (3);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms3_3
OUTFILE= ".\un_sas_output\Covparms3_3.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat3_3
OUTFILE= ".\un_sas_output\FitStat3_3.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge3_3
OUTFILE= ".\un_sas_output\Converge3_3.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov3_3
OUTFILE= ".\un_sas_output\Asycov3_3.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms3_4;
ods output FitStatistics=FitStat3_4;
ods output ConvergenceStatus=Converge3_4;
ods output asycov=Asycov3_4;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (3) and Gen in (4);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms3_4
OUTFILE= ".\un_sas_output\Covparms3_4.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat3_4
OUTFILE= ".\un_sas_output\FitStat3_4.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge3_4
OUTFILE= ".\un_sas_output\Converge3_4.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov3_4
OUTFILE= ".\un_sas_output\Asycov3_4.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms3_5;
ods output FitStatistics=FitStat3_5;
ods output ConvergenceStatus=Converge3_5;
ods output asycov=Asycov3_5;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (3) and Gen in (5);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms3_5
OUTFILE= ".\un_sas_output\Covparms3_5.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat3_5
OUTFILE= ".\un_sas_output\FitStat3_5.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge3_5
OUTFILE= ".\un_sas_output\Converge3_5.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov3_5
OUTFILE= ".\un_sas_output\Asycov3_5.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




ods output covparms=Covparms3_6;
ods output FitStatistics=FitStat3_6;
ods output ConvergenceStatus=Converge3_6;
ods output asycov=Asycov3_6;
proc mixed data=stdmrwings method=reml asycov lognote scoring=5;
where multiout not in (1) and Treat in (3) and Gen in (6);
class Treat Line Animal Vial Trait;
model Score = Trait;
random Trait / subject= Line type=un;
random Trait / subject= Vial(Treat Line) type=un;
repeated / subject = Animal(Treat Line Vial) type=un;
run;
quit;

PROC EXPORT DATA= WORK.Covparms3_6
OUTFILE= ".\un_sas_output\Covparms3_6.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.FitStat3_6
OUTFILE= ".\un_sas_output\FitStat3_6.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Converge3_6
OUTFILE= ".\un_sas_output\Converge3_6.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;
PROC EXPORT DATA= WORK.Asycov3_6
OUTFILE= ".\un_sas_output\Asycov3_6.csv"
DBMS=CSV REPLACE;
PUTNAMES=YES;
RUN;




