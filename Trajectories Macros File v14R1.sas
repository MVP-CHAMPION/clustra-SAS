/*Trajectories macro version 14 -- 2-14-2022*/

/*Graphics Setup Program for Trajectories Macros

The following code should be run before running any of the trajectory macros as it
sets up graphic styles, devices and default fonts, line styles and colors.

The GRAPHSET macro below should be run with your system settings.
*/

/*Create Graphics Template Defaults*/
/*You only have to ever run PROC TEMPLATE once on a system to create the style
  if you have write permission to where styles are stored.  Otherwise, once per session */
proc template;
	define style Styles.splinecurve;
	parent = Styles.Listing;
	class GraphWalls /
		linethickness = 3px
		linestyle = 1
		frameborder = off
		contrastcolor = GraphColors('gaxis')
		backgroundcolor = white
		color = GraphColors('gwalls');
	class GraphAxisLines /
		tickdisplay = "outside"
		linethickness = 3px
		linestyle = 1
		contrastcolor = GraphColors('gaxis')
		color = GraphColors('gaxis'); 
	class GraphFonts "Fonts used in graph styles" /
		'NodeDetailFont' = ("<sans-serif>, <MTsans-serif>",14pt)
		'NodeLinkLabelFont' = ("<sans-serif>, <MTsans-serif>",18pt)
		'NodeInputLabelFont' = ("<sans-serif>, <MTsans-serif>",18pt)
		'NodeLabelFont' = ("<sans-serif>, <MTsans-serif>",18pt)
		'NodeTitleFont' = ("<sans-serif>, <MTsans-serif>",18pt)
		'GraphDataFont' = ("<sans-serif>, <MTsans-serif>",14pt)
		'GraphUnicodeFont' = ("<MTsans-serif-unicode>",18pt)
		'GraphValueFont' = ("<sans-serif>, <MTsans-serif>",14pt)
		'GraphLabel2Font' = ("<sans-serif>, <MTsans-serif>",18pt)
		'GraphLabelFont' = ("<sans-serif>, <MTsans-serif>",18pt,bold)
		'GraphFootnoteFont' = ("<sans-serif>, <MTsans-serif>",10pt)
		'GraphTitleFont' = ("<sans-serif>, <MTsans-serif>",15pt,bold)
		'GraphTitle1Font' = ("<sans-serif>, <MTsans-serif>",18pt,bold)
		'GraphAnnoFont' = ("<sans-serif>, <MTsans-serif>",18pt)
	;
	class body /  
		fontfamily = "Arial, sans-serif"
		color = #000000
		backgroundcolor = #FFFFFF
	;  
	class Header/
 		background=white
 		foreground=black
	;  
end; 
run; 

/***********************************************************************************************/
%macro graphset(gpath,dpi,style,format=emf,type=listing);
/*GRAPHICS OUTPUT OPTIONS
Not all combinations of format and type are supported
TYPE            FORMATS
Listing         PNG (default), GIF, BMP, DIB, EMF, EPSI, GIF, JFIF, JPEG, PBM, PS, TIFF, WMF
HTML            PNG (default), GIF, JPEG,
LATEX           PostScript (default), EPSI, GIF, JPEG, PDF, PNG
PDF             PDF

Parameters:
gpath       Location where graphic output will be stored.  This is a filename for PDF output
            and a directory for all other output types.
dpi         Dots per inch.  Not relevant for PDF or vector graphics [EMF]
style       Which style to use.  "splinecurve" is recommended
format      Which graphic format to use.  Default is "EMF'.
type        Which output destination to use.  Default is "Listing".
fname		Allowing for a filename to be used when specifying the output.

Tested combinations include Listing/EMF, Listing/JPEG, HTML/JPEG and PDF/PDF.  Note that
PDF files can be converted to PS or EPS output by Acrobat Pro or Inkscape
ODS1-ODS3 are needed in a macro before output is produced.
ODS4 is needed at the end of the output, especially if PDF output is desired.
     to close the PDF file.
*/
	%global ods1 ods2 ods3 ods4;

	%if &style eq 	%then %let style=default;
	%if &dpi eq  	%then %let dpi=300;

	%let 	ods1= ods &type gpath="&gpath" style=&style	image_dpi = &dpi; 
	%let 	ods2= ods graphics on / width=10in height=6.8in IMAGEFMT=&format;/*reset=all*/

	run;

	%if %upcase(&type) eq PDF %then %do;
		%let ods1= ods pdf file="&gpath" style=&style;
		%let ods3= options orientation=landscape nodate nonumber;
		%let ods4= ods pdf close;
	%end;
%mend graphset;

/*  Sample executions of %graphset:
%graphset(C:\Users\gagnon\Documents\figure output\test.pdf,,splinecurve,format=pdf,type=pdf);
%graphset(C:\Users\gagnon\Documents\figure output,,splinecurve,format=emf,type=listing);
%graphset(N:\Figure Output,,splinecurve,format=jpeg,type=listing);
*/
/***********************************************************************************************/

/*Set default graphic attributes for spline macros*/
data graphattr;
/*scatter and series use linecolor, band uses fillcolor*/
	length ID $ 7  linecolor $ 11 fillcolor $ 11  ;
	input value  linecolor  fillcolor linethickness ;

	markersymbol="squarefilled";
	markersize=40;
	filltransparency = .60;
	linepattern = 'Solid';

	if _n_ le 12 then ID = 'color';
	else ID='grays';

	datalines;
 1  black       black      2
 2  red         red        2
 3  blue        blue       2 
 4  green       green      2 
 5  cxD17800    cxD17800   2
 6  cxB26084    cxB26084   2
 7  cx2597FA    cx2597FA   2
 8  brown       brown      2
 9  yellow      yellow     2
 10 gray        gray       2
 11 lightred    lightred   2
 12 lightgreen  lightgreen 2
 1  black       black      2
 2  gray        cxEEEEEE   3
 3  black       black      5
 4  gray        cxEEEEEE   5
;
run;

/*****************************************************************************************/
/*
CLUSSMOOTH defines the form of smoothing used in TRAJLOOP and other macros.  The METHOD
is set in the TRAJSETUP macro for TRAJLOOP.

*/
%macro clussmooth(method,dasn,idvar,wtvar,timvar,rskvar,dfvar,outdsn,numvar);
/*

    METHOD      Define the smoothing method to use. See below.
    DASN        Input data set name for procedure
    IDVAR       Subject ID variable.  &GID for %TRAJLOOP.
    WTVAR       Weight variable for the regression, which determines the subgroup
                to be analyzed.  
    TIMVAR      Independent time variable for the spline.
    RSKVAR      Dependent variable for the spline
    DFVAR       For most methods, defines the degrees of freedom to use
    OUTDSN      Output data set name for predicted values
    NUMVAR      Numerical variable that identifies variable names for output data sets
                and predicted values.
*/
	%if &method eq 0 %then %do;
	    proc gampl data=&dasn ;
		    title "GAMPL with default optimization.";
		    weight &wtvar ;
		    id &idvar  &timvar &rskvar &wtvar;
		    model &rskvar =  spline(&timvar /maxdf=&dfvar);
		    output out= &outdsn p=pred&numvar r=resid;
		    performance threads=8;
	    run;quit;
	%end;
	%if &method eq 1 %then %do;
	    proc gampl data=&dasn smoothoptions(technique=quanew);
		    title "GAMPL with dual quasi-Newton optimization.";
		    weight &wtvar ;
		    id &idvar  &timvar &rskvar &wtvar;
		    model &rskvar =  spline(&timvar /maxdf=&dfvar);
		    output out= &outdsn p=pred&numvar r=resid;
		    performance threads=8;
	    run;quit;
	%end;
	%if &method eq 2 %then %do;
	    proc gampl data=&dasn  smoothoptions(technique=dbldog);
		    title "GAMPL with double-dogleg optimization";
		    weight &wtvar;
		    id &idvar  &timvar &rskvar &wtvar ;
		    model &rskvar =  spline(&timvar /maxdf=&dfvar);
		    output out= &outdsn p=pred&numvar r=resid;
		    performance threads=8;
	    run;quit;
	%end;
	%if &method eq 3 %then %do;
	    proc gampl data=&dasn smoothoptions(technique=congra);
		    title "GAMPL with conjugate-gradient optimization.";
		    weight &wtvar;
		    id &idvar  &timvar &rskvar &wtvar ;
		    model &rskvar =  spline(&timvar /maxdf=&dfvar);
		    output out= &outdsn p=pred&numvar r=resid;
		    performance threads=8;
	    run;quit;
	%end;
	%if &method eq 4 %then %do;
	    proc gampl data=&dasn smoothoptions(technique=nmsimp);
		    title "GAMPL with Nelder-Mead simplex optimization";
		    weight &wtvar;
		    id &idvar  &timvar &rskvar &wtvar;
		    model &rskvar =  spline(&timvar /maxdf=&dfvar);
		    output out= &outdsn p=pred&numvar r=resid;
		    performance threads=8;
	    run;quit;
	%end;
	%if &method eq 5 %then %do;
	    proc gam data=&dasn;
		    title "GAM";
		    where &wtvar eq 1;
		    model &rskvar =  spline(&timvar,df=&dfvar); 
		    score data=&dasn out=&outdsn ;
	    run;
	    data &outdsn; set &outdsn;
		    rename P_&rskvar = pred&numvar;
		    keep P_&rskvar &idvar  &timvar &rskvar &wtvar;
	    run;quit;
	%end;
	%if &method eq 6 %then %do;
	    proc glimmix data=&dasn method=laplace empirical;    
		    title "GLIMMIX with EFFECT";
		    effect spl = spline(&timvar /naturalcubic knotmethod=percentiles(&dfvar)) ; 
		    weight &wtvar;
		    id &idvar  &timvar &rskvar &wtvar;
		    model &rskvar = spl ;
		    output out= &outdsn predicted=pred&numvar residual=resid;
	    run;quit;
	%end;
	%if &method eq 7 %then %do;
	    proc tpspline data=&dasn ;    
		    title "TPSPLINE";
		    where &wtvar eq 1;
		    id &idvar  &timvar &rskvar ;
		    model &rskvar = (&timvar) /D=&dfvar;
		    score data=&dasn out=&outdsn ;
	    run;quit;
	    data &outdsn; set &outdsn;
		    rename P_&rskvar = pred&numvar;
		    keep P_&rskvar &idvar  &timvar &rskvar &wtvar;
	    run;
	%end;
	%if &method eq 8 %then %do;
	    proc transreg data=&dasn;
		    title 'TRANSREG';
		    weight &wtvar;
		    id &idvar;
		    model identity(&rskvar) = pbspline(&timvar/gcv);
		    output out= &outdsn pprefix=pred rprefix=resid;
	    run;quit;
	    data &outdsn; set &outdsn;
		    rename pred&rskvar=pred&numvar resid1=resid;
		    keep &idvar &timvar &rskvar pred&rskvar &wtvar resid1;
	    run;
	%end;
%mend clussmooth;
run;
/*****************************************************************************************/
/*
TRAJSETUP needs to be run first for any clustering attempts
It creates global macro variables that are used elsewhere as well
*/

%macro trajsetup(dsn,id,time,riskvar,ngroups,maxdf,ptrim,seed=1234,steps=1,method=0,random=YES);
/* This macro sets up the initial conditions for K-means clustering of trajectories.
   Execute this before the TRAJLOOP macro.  This macro can also be run with the RANDOM=NO
   option to re-initialize the clustering parameters and restart without the initial 
   randomization step.

Parameters:
   DSN          Input data set.  Requires ID, TIME and RISKVAR variables.
                One record per observation.  If restarting, use the CLUSOUT permanent
                data set created before shutdown.
   ID           Subject identifying variable.  Must be non-missing.
   TIME         X axis variable.  Can be a time or age variable. Must be non-missing.
   RISKVAR      Y axis variable.  It should be continuous. Must be non-missing.
   NGROUPS      Number of groups to create.
   MAXDF        Maximum degrees of freedom for the thin plate regression splines.
   PTRIM        Percentage of subjects to exclude from each iteration.  Good for
                reducing the effect of outlier observations. Trimmed subjects can
                return at later iteration if the group curves approach trimmed subjects.
   SEED         Random number seed used for initial cluster assignment.  Values < 1 will
                assign a random seed.
   STEPS        Is this a one-step [1] or two-step [2] clustering algorithm?  A two-step
                algorithm prepares for hierarchical clustering.
   METHOD       Selects spline technique to use.  See CLUSSMOOTH macro for details.
   RANDOM       If initiating cluster analysis, RANDOM=YES.  If restarting after a SAS 
                shutdown, reusuing the last produced iteration, RANDOM=NO.


*/
	options validvarname=V7; *check with David what this is?;

/* Setup the tracking dataset to see how things change between iterations */
	data tracking;
		dummy = 1;
	run;

/*Remove records with missing data*/
	proc sort data=&dsn; by &id &time &riskvar;
	data ranstart; set &dsn;
		by &id &time &riskvar;
		if cmiss(&id,&time,&riskvar) eq 0;
	run;

/*Create macro variable 'records' to count observations
  KEEPIT variable defines how many subjects to keep (based on trimming)*/
	data _null_;
		set ranstart (obs=1) nobs=nobs;
		call symput('records',nobs);
	run;
	%let keepit=%sysevalf(&records*(1-&ptrim/100));

/* Randomly assign the first group using the specified seed */
	%if &random eq YES %then %do;
		data ranstart; 
			set ranstart;
			by &id &time &riskvar;
			retain group;
			call streaminit(&seed);
			if first.&id then group = ceilz(&ngroups*rand('uniform'));
			%do ngrp = 1 %to &ngroups;
				grp&ngrp = group eq &ngrp;
			%end;
		run;
	%end;

/*Create global variables that will apply to TRAJLOOP macro*/
	%global gid gtime griskvar gngroups gmaxdf gptrim gkeepit gloop gsteps gmethod;
	%let gid=&id;
	%let gtime=&time;
	%let griskvar=&riskvar;
	%let gngroups=&ngroups;
	%let gmaxdf=&maxdf;
	%let gptrim=&ptrim;
	%let gkeepit=&keepit;
	%let gloop = 0;
	%let gsteps = &steps;
	%let gmethod = &method;
	run;
%mend trajsetup;

/**************************************************************************/

/* The trajloop macro does the K-means clustering algorithm */

%macro trajloop(outlib,iter,minchange,minsubs,min_x,max_x,by_x,min_y,max_y,by_y,showall=NO,showany=YES,AIC=YES);
/* Macro for iterating the K-means clustering of longitudinal trajectories.
   Please run TRAJSETUP first to initialize the system.  
   TRAJLOOP can be run multiple times to continue iterating.

Parameters:
   OUTLIB               Libname for directory storing output data sets.  If blank,
                        default is WORK
   ITER                 Maximum iterations of K-means algorithm
   MINCHANGE            Minimum percentage of subjects changing groups.  
                        When the percentage changing group at an iteration falls below
                        this threshold, the algorithm will stop iterating, regardless 
                        of the value of ITER.  Set to zero if you want to continue
                        iterating to the full value of ITER.  A minimum of one iteration
                        will occur with each invocation of the macro.  This is a percentage.
   MINSUBS              For two-stage clustering, indicate the minimum cluster size to pass 
                        on to the hierarchical clustering macro
   MIN_X, MAX_X, BY_X   minimum, maximum and interval values on X axis of plot
   MIN_Y, MAX_Y, BY_Y   minimum, maximum and interval values of Y axis of plot
   SHOWALL              Default is NO.  Any other value will display all output from
                        PROC GAMPL and will produce a graphic for each iteration.
                        If NO, only the final graphic will display along with the
                        change history.
   SHOWANY              Show any graphics?  Default is YES.  Anything else allows graphics
   AIC                  If YES, it calculates fit statistics and silhouette statistics. Also
                        generates FINALFIT and SILHOUT data sets used by the TRAJFIT macro.

Output Data Sets Stored in OUTLIB:
   CLUSOUT              Output from cluster analysis with group identifiers
   GRAPHOUT             Data set suitable for graphing trajectories
   MERGIT               Data set with last set of distances at observation level, 
                        useful for diagnostics
   RMS                  Data set with last set of mean squared distances at subject level
                        Used for silhouette plots
   WIDE                 Data set used for two-stage clustering -- input for TRAJHCLUS
   TRACKING             Store information on clustering history 


*/
	%if &outlib eq  %then %let outlib = work;

	%let stopit = NO;
	%let loopend = ;
	%let maxloop = %eval(&gloop+&iter);

	%do loopit = &gloop %to &maxloop;
  		%if (&stopit eq YES) AND (&loopend eq ) %then %do;
       		%let loopend = %eval(&loopit - 1);
       		title "Iteration &loopend";
  		%end;

  		%if &stopit eq NO %then %do; /*Start of the loop where groups are assigned */
 	  		%if &showall eq NO %then %do; 
				ods listing close; 
			%end;
        	title "Iteration &loopit";
			%do ngrp = 1 %to &gngroups;

			/*Trap groups with sample size less than max degrees of freedom*/
  				ods output onewayfreqs=minsamp;
  				proc freq data=ranstart;
				/*  title "Error capture iteration &loopit, group &ngrp";*/
				    tables grp&ngrp;
    			run;

				data minsamp; 
					set minsamp end=eof;
    				output;
  					if eof and (grp&ngrp ne 1) then do;
    					grp&ngrp=1; 
						frequency=0; 
						output; 
					end;
				run;
  				data minsamp; 
					set minsamp;
  					if grp&ngrp eq 1;
  					if Frequency lt &gmaxdf then call symput('singular',1);
    				else call symput('singular',0);;
    			run;

				%if &singular eq 0 %then %do;
    				%clussmooth(method=&gmethod, dasn=ranstart, idvar=&gid, wtvar=grp&ngrp,
								timvar=&gtime, rskvar=&griskvar, dfvar=&gmaxdf,
								outdsn=grpout&ngrp, numvar=&ngrp);
					run;
    			%end;
    			%else %do;
    				data grpout&ngrp; 
						set grpout&ngrp;
    					pred&ngrp = .;
    				run;
    			%end;
			%end;

			data mergit; 
				merge grpout1-grpout&gngroups; 
				by &gid &gtime &griskvar;
  				array preds (&gngroups) pred1-pred&gngroups;
  				array mss (&gngroups) mss1-mss&gngroups;
  				do i = 1 to &gngroups;
    				mss(i) = (preds(i) - &griskvar)**2;
  				end;
  				keep &gid &gtime &griskvar mss1-mss&gngroups pred1-pred&gngroups;
  			run;

  			proc means data=mergit noprint; 
				by &gid;
  				var mss1-mss&gngroups;
  				output out=rms mean=mean1-mean&gngroups;
  			run;

  			data rms; 
				set rms;
  				array mss (&gngroups) mean1-mean&gngroups;
  				minmean = min(of mean1-mean&gngroups);
  				do i = 1 to &gngroups;
    				if minmean = mss(i) then newgroup = i;
  				end;
  				keep &gid newgroup minmean mean1-mean&gngroups;
  			run;

  			data ranstart; 
				merge ranstart  rms; *do we need to sort the datasets beforehand?;
				by &gid;
  				%do ngrp = 1 %to &gngroups;  /*why is this in a %do loop and not regular do loop?;*/
     				grp&ngrp = newgroup eq &ngrp;
  				%end;
  				oldgrp = group;
  				if group = newgroup then change = 0;
  				else change = 1;
  				group=newgroup;
  				drop newgroup;
  			run;

			proc freq data= ranstart; /*Can we suppress this output? */
  				tables oldgrp*group ;
  			run;

 		 	/*Trimming*/
  			%if &gptrim ne 0 %then %do;
    			proc sort data=ranstart; by minmean;
    			data ranstart; 
					set ranstart; 
					by minmean;
    				array ngrps (&gngroups) grp1-grp&gngroups;

    				if _n_ gt &gkeepit then do i= 1 to &gngroups;
      					ngrps(i) = 0;
      					group = .;
    				end;
    			run;
    			proc sort data=ranstart; 
					by &gid &gtime &griskvar;
    			run;
  			%end;

		/* Calculate Changes From Last Iteration -- modified to one record/subject */
    		proc sort data=ranstart; 
				by &gid;
    		run;
    		data difcheck; 
				set ranstart; 
				by &gid;
    			if first.&gid;
    		run;

    		ods listing close;
    		ods output OneWayFreqs=owf;
   	 		proc freq data=difcheck;
    			where (group ne .) and (oldgrp ne .);
    			tables change;
    		run;
    		quit; *is this needed?;
    		ods listing;

			data owf; 
				set owf;

    			if 	((change eq 1) and (Percent le &minchange)) or 
       				((change eq 0) and (Percent eq 100)) 		then call symput('stopit','YES');

				if ((change eq 0) and (Percent eq 100)) then do;
       				change = 1;
       				Frequency = 0;
       				Percent = 0;
    			end;
			run;

			data tracking; 
				set tracking owf;
    			if change eq 1;
    		run;
		%end;

		/*Graphics section*/
		%if ((&showall ne NO) or (&loopit ge &maxloop)) and (&showany eq YES) %then %do;
    		data graphin; 
				set grpout1-grpout&gngroups;
    			array grps (&gngroups) grp1-grp&gngroups;
    			array preds (&gngroups) pred1-pred&gngroups;
    			do i = 1 to &gngroups;
      				if grps(i) eq 1 then do;
       					group=i;
       					pred=preds(i);
      				end;
    			end;
    			if group ne . then output;
    			keep pred group &gtime;
    		run;

    		proc sort data=graphin nodupkey out=graphout; 
				by group &gtime;
    		run;

			%trajplot(dsn=graphout, rvar=pred, groups=group, timevar=&gtime, attrid=color, 
						xlabel=Time, ylabel=Predicted, /*could consider parameterizing this*/
						min_x=&min_x, max_x=&max_x, by_x=&by_x,
						min_y=&min_y, max_y=&max_y, by_y=&by_y);
			run;

			data &outlib..graphout; 
				set graphout;
			run;
		%end;

	%end;

	%if &loopend eq  %then %let loopend = %eval(&loopit - 1);
	run;

	proc print data=tracking;
		title 'Summary of Group Changes by Iteration';
		where dummy ne 1;
		var Frequency Percent;
	run;

/*Setup predicted values for hierarchical clustering*/
	%if &gsteps eq 2 %then %do;
		proc means data=ranstart;
			var &gtime;
			output out=hdata min=minx max=maxx;
		run;
		data hdata; 
			set hdata;
			lrate = (maxx-minx)/100;
			array ggroups (&gngroups) grp1-grp&gngroups;
			do i = 1 to &gngroups;
			  	ggroups(i) = -1;
			end;
 			&griskvar = .;
			do &gtime = minx to maxx by lrate;
  				output;
			end;
		run;

		data ranoutt; 
			set ranstart hdata;
		run;

		%do ngrp = 1 %to &gngroups;
		/*Trap groups with sample size less than max degrees of freedom*/
  			ods output onewayfreqs=minsampf;
  			proc freq data=ranoutt;
    			title "Error capture final, group &ngrp";
    			tables grp&ngrp;
    		run;
  			data minsampf; 
				set minsampf end=eof;
    			output;
  				if eof and (grp&ngrp ne 1) then do;
    				grp&ngrp=1; 
					frequency=0; 
					output; 
				end;
    		run;
  			data minsampf; 
				set minsampf;
  				if grp&ngrp eq 1;
  				if Frequency lt &gmaxdf then call symput('singf',1);
    			else call symput('singf',0);;
    		run;

   			%if &singf eq 0 %then %do;
    			%clussmooth(method=&gmethod, dasn=ranoutt, idvar=&gid, wtvar=grp&ngrp,
							timvar=&gtime, rskvar=&griskvar, dfvar=&gmaxdf,
							outdsn=hpout&ngrp, numvar=&ngrp);
				run;
				data hpout&ngrp; 
					set hpout&ngrp;
  					group = &ngrp;
  					pred = pred&ngrp;
  				run;
				/*RMSSTD calculation*/
   				proc means noprint data=hpout&ngrp;
  	 				where grp&ngrp eq 1;
   					id group;
   					var resid;
   					output out=resid&ngrp std=std;
   				run;

   				data hpout&ngrp; 
					set hpout&ngrp;
     				where grp&ngrp eq -1;
     				group = &ngrp;
     				obs = _n_;
   				run;
			%end;
   			%else %do;
			/*RMSSTD calculation*/
   				data resid&ngrp;
   					std = .;
   					group = &ngrp;
				run;
   				data hpout&ngrp; 
    				grp&ngrp eq -1;
     				group = &ngrp;
     				obs = 0;
   				run;
    		%end;
		%end;

		/*FREQ calculation*/
  		proc sort data=ranoutt; by &gid;
		data unique; 
			set ranoutt; 
			by &gid;
			if first.&gid;
			if &gid ne .;
			keep &gid group &griskvar;
		run;
		proc means noprint data=unique;
			class group;
			var &griskvar;
			output out=nout n=obsnum;
		run;

  		data allhout;
			set hpout1-hpout&gngroups;
			by group obs;
  		run;
  		data allresid; 
			set resid1-resid&gngroups; 
			by group;
		run;

		proc transpose data=allhout out=wide1 prefix=&griskvar; 
			by group;
			id obs;
			var pred;
		run;
  
		data &outlib..wide; 
			merge wide1 allresid nout; 
			by group;
			/*Remove clusters with fewer than &minsubs*/
			if obsnum ge &minsubs;
			if group ne .;
			drop _freq_;
			rename std=_rmsstd_ obsnum=_freq_;
		run;

	options ls=80 ps=60; /*why is this here?*/
		proc means data=&outlib..wide; /*do we need to print this?*/
		run;
	run; /*why is this here?*/
	%end;

/* Create output data sets */
	data &outlib..clusout; set ranstart;
	data &outlib..graphout; set graphout;
	data &outlib..mergit; set mergit;
	data &outlib..rms; set rms;
	data &outlib..tracking; set tracking;
	run;

	%if &AIC eq YES %then %do; /* Wouldn't this make more sense to say you want to do the plots regardless of AIC?*/
		/* Run Silhouette Macro for Elbow and Silhouette Statistics */
		%silhouette(lib=&outlib,sid=&gid,groupsn=&gngroups);
		run;
		data silout; 
			set silout;
			loops = &loopend;
			groups = &gngroups;
		run;

/* Run Joint Model for AIC*/
		data ranfinish; 
			set ranstart;
  			%do ngrp = 1 %to &gngroups;
  				time&ngrp = &gtime*grp&ngrp;
  			%end;
		run;

		ods output FitStatistics=finalfit;
    	proc gampl data=ranfinish;
    		title "global model for iteration &gloop with &gngroups groups AIC";
    		id &gid  &griskvar ;
    		model &griskvar =  
				%do ngrp = 1 %to &gngroups;
    				spline(time&ngrp /maxdf=&gmaxdf)
  				%end;
			;
    		output out= final p=predfinal;
    		performance threads=8;
    	run;
		quit;

		data finalfit; 
			set finalfit;
			loops = &loopend;
			groups = &gngroups;
		run;

		proc print data=finalfit;
			Title "fit statistics for joint model, &loopend iterations";
		run;

		proc print data=silout;
			title "Silhouette and Elbow values for clusters, &gngroups groups";
		run;
	%end;

%mend trajloop;


/***********************************************************************************************/

%macro trajplot(dsn,rvar,groups,timevar,attrid,xlabel,ylabel,min_x,max_x,by_x,min_y,max_y,by_y);
/*  Macro for printing cluster analysis results.  May be used for additional copies of 
    trajectory output without re-running the TRAJLOOP macro.  This macro is called within TRAJLOOP.

DSN                     Data set name of printable data set.  Use full name [libname.dataset]
                        Usually this is the data set GRAPHOUT produced by TRAJLOOP.
RVAR                    Risk factor variable used for clustering
GROUPS                  Grouping variable for clustering.  Usually GROUP [one-step models] or CLUSTER
                        [two-step models].
TIMEVAR                 Time variable used in clustering
ATTRID                  Either "color" for color or "grays" for gray-scale plots
XLABEL                  Label for time variable
YLABEL                  Label for risk factor variable
MIN_X, MAX_X, BY_X      These are minimum, maximum and by variables for the X axis [time]
MIN_Y, MAX_Y, BY_Y      These are minimum, maximum and by variables for the Y axis [risk factor]

*/
	&ods1; &ods2; &ods3;
	proc sort data=&dsn; 
		by &groups &timevar &rvar;
	run;
	ods graphics / imagename="TrajPlotSAS"; /*adding image name 9/11/23... may want to add reset=index to stop accumulating unwanted graphics */
	proc sgplot data=&dsn /*noautolegend*/ dattrmap=graphattr; /*commented out the noautolegend for ease of graphics interpretability*/
		title ;
		series y=&rvar x=&timevar / group=&groups attrid=&attrid;
		xaxis  values=(&min_x to &max_x by &by_x) label=  "&xlabel";
		yaxis  values=(&min_y to &max_y by &by_y) label= "&ylabel";
	run;
	quit;

	/* commenting the section out due to un-needed graphics at the moment */
	/*
	ods graphics / imagename="Histogram_&timevar.";
	proc sgpanel data=&dsn dattrmap=graphattr;
		panelby &groups;
		histogram &timevar / scale=count fillattrs=(color=gray transparency=0.6) binwidth=&by_x ;
	run;quit;
	*/

	&ods4;
	run;
%mend trajplot;
/******************************************************************************/

%macro silhouette(lib,sid,groupsn,cgroup=group);
/* Macro for producing silhouette plots.  This macro uses output from the TRAJLOOP macro

    LIB         Libname identifying where TRAJLOOP output data sets are stored
    SID         Subjects unique identifier
    GROUPSN     Number of clusters produced.  Should match TRAJLOOP macro call
    CGROUP      Cluster group variable name.  Should be "group"


*/
	data membership &lib.clusprime; 
		set &lib..clusout (keep=&sid &cgroup);
	run;
	proc sort data=membership; by &sid;
	data membership; 
		set membership; 
		by &sid;
		if first.&sid;
		keep &sid &cgroup;
		rename &cgroup=primegroup;
	run;
	proc sort data=&lib..rms; by &sid;
	data newrms; 
		merge &lib..rms membership; 
		by &sid;
		array dist (&groupsn) mean1-mean&groupsn;
		ingrpmean = dist(primegroup);
		dist(primegroup)=.;
		second = min(of dist(*));
		do i = 1 to &groupsn;
			if dist(i) = second then group2 = i;
		end;
		silhouette = (second-ingrpmean)/(max(ingrpmean,second));
		start = 0;
		if silhouette eq . then loser = .;
		else if silhouette lt 0 then loser = 1;
		else loser = 0;
		label group2="Nearest Group";
	run;
	proc means data=newrms noprint; /*added noprint option 9-11-23;*/
		var ingrpmean;
		output out=meandist mean=meandist sum=sumdist;
	run;
	data _null_; set meandist;
		call symputx("elbow&GROUPSN",meandist,'G');
	run;
	proc sort data=newrms; by primegroup descending silhouette;
	run;
	options ls=80 ps=55;
	proc means data=newrms min mean median max noprint;  /*added noprint option 9-11-23;*/
		var silhouette;
		output out=tmp1x mean=msil median=mdsil;
	run;
	data tmp1x; 
		set tmp1x;
		call symputx("meansil&GROUPSN",msil,'G');
		call symputx("medsil&GROUPSN",mdsil,'G');
		call symputx("meansil_p",round(msil,0.001)); /*hg added for graphics */
		call symputx("medsil_p",round(mdsil,0.001)); /*hg added for graphics */
	run;
	data plotdat; 
		set newrms; 
		by primegroup descending silhouette;
		xvar = _n_;
		label primegroup='cluster' xvar='observation #';
	run;
	data silout; 
		meansil = &&meansil&GROUPSN;
		medsil = &&medsil&GROUPSN;
		elbow = &&elbow&GROUPSN;
	run;
	&ods1;
	&ods2;
	&ods3;
	ods graphics / imagename="Silhouette" reset=index; /*adding reset=index to clean up output file names*/
	proc sgplot data=plotdat;
		title1 'Silhouette Plot';
		title2 h=1 "Mean=&meansil_p (Black), Median=&medsil_p (Red)"; /*changed printing output*/
		band x=xvar upper=silhouette lower=start /group=primegroup;
		refline &&meansil&GROUPSN /axis=y lineattrs=(color=black);
		refline &&medsil&GROUPSN /axis=y lineattrs=(color=red);
	run;
	&ods4;
	run;

	/* adding output of data to streamline analytic workflow (hg 9/11/23) */
	data &lib..sil_plotdat; set plotdat;
	data &lib..silout; set silout;
	run;

	/* commenting out this code because it doesn't seem as helpful for the manuscript */
/*
	proc freq data=newrms;
		title;
		tables loser*primegroup;
		tables primegroup*group2;
	run;
	proc means data=plotdat min mean median max std; 
		class primegroup;
		var silhouette;
	run;
*/
%mend silhouette;

/*****************************************************************************************/


%macro trajfit(nstart,nstop,restart,tlibname,dsn,id,time,riskvar,maxdf,ptrim,iter,minchange);

/*This macro will run the trajectories macros multiple times, saving the AIC*/
/*values for comparisons.

nstart          Starting value for number of groups.  This should be >= 2.  If restarting
                an interupted series, simply start where the last run left off.
nstop           Ending value for the number of groups.  Must be >= nstart.
restart         if YES, delete previous results and start a new set.
tlibname        Libname for storing results.  If blank, it defaults to WORK.

TRAJSETUP macro variables:
DSN             Input data set.  Requires ID, TIME and RISKVAR variables.
                One record per observation
ID              Subject identifying variable.  Must be non-missing.
TIME            X axis variable.  Can be a time or age variable. Must be non-missing.
RISKVAR         Y axis variable.  Currently, it should be continuous. Must be non-missing.
MAXDF           Maximum degrees of freedom for the thin plate regression splines.
PTRIM           Percentage of subjects to exclude from each iteration.  Good for
                reducing the effect of outlier observations. Trimmed subjects can
                return at later iteration if the group curves approach trimmed subjects.

TRAJLOOP macro variables:
ITER            Maximum iterations of K-means algorithm
MINCHANGE       Minimum percentage of subjects changing groups.  
                When the percentage changing group at an iteration falls below
                this threshold, the algorithm will stop iterating, regardless 
                of the value of ITER.  Set to zero if you want to continue
                iterating to the full value of ITER.  A minimum of one iteration
                will occur with each invocation of the macro.


*/
	%if &tlibname eq  %then %let tlibname = work;

   	%if &restart eq YES %then %do;
     	data &tlibname..results;
     		loops = 0;
     		groups = 0;
     	run;
     	data &tlibname..silout;
     		loops = 0;
     		groups = 0;
     	run;
	%end;

	%do groupsw = &nstart %to &nstop;
		%trajsetup(dsn=&dsn, id=&id, time=&time, riskvar=&riskvar, ngroups=&groupsw,
					maxdf=&maxdf, ptrim&ptrim, steps=1, method=0, random=YES);
		run;

		%trajloop(outlib=&tlibname, iter=&iter, minchange=&minchange, minsubs=0,
          			min_x=20, max_x=110, by_x=20,
					min_y=100, max_y=180, by_y=20,
					showall=NO, showany=NO, AIC=YES);
		run;

		data &tlibname..results; 
			set &tlibname..results finalfit;
		run;
		data &tlibname..silout; 
			set &tlibname..silout silout;
		run;
	%end;

	proc sort data=&tlibname..results; by groups;
	data result; 
		set &tlibname..results; 
		by groups;
		if groups ne 0;
		retain AIC AICC BIC GCV PLL RP EDF EDFE;

		array tnames (8) AIC AICC BIC GCV PLL RP EDF EDFE;
		if first.groups then do i = 1 to 8;
		  	tnames(i)=.;
		end;

		if Description =: 'AIC ' then AIC = Value;
		else if Description =: 'AICC' then AICC = Value;
		else if Description =: 'BIC' then BIC = Value;
		else if Description =: 'GCV' then GCV = Value;
		else if Description =: 'Pen' then PLL = Value;
		else if Description =: 'Roug' then RP = Value;
		else if Description = 'Effective Degrees of Freedom' then EDF = Value;
		else if Description =: 'Effective Degrees of Freedom for Error' then EDFE = Value;

		keep AIC AICC BIC GCV PLL RP EDF EDFE groups;
		if last.groups;
	run;
	proc sort data=result; by groups;
	proc sort data=&tlibname..silout; by groups;
	data &tlibname..silout; 
		set &tlibname..silout; 
		by groups;
		if groups ne 0;
	run;
	data result; merge 
		result &tlibname..silout; 
		by groups;
	run;
	options ls=80 ps=60;
	proc print data=result;
		id groups;
	run;

	ods graphics / imagename="Silmean";
	proc sgplot data=result noautolegend dattrmap=graphattr;
		title "Mean Silhouette Plot";
		series y=meansil x=groups /  attrid=color lineattrs=(thickness=5);
		xaxis   label=  "Groups";
		yaxis   label= "Mean Silhouette Value";
	run;
	quit;

	ods graphics / imagename="Elbow";
	proc sgplot data=result noautolegend dattrmap=graphattr;
		title "Elbow Plot";
		series y=elbow x=groups /  attrid=color lineattrs=(thickness=5);
		xaxis   label=  "Groups";
		yaxis   label= "Mean Distance";
	run;
	quit;

	ods graphics / imagename="AIC";
	proc sgplot data=result noautolegend dattrmap=graphattr;
		title "AIC Plot";
		series y=AIC x=groups /  attrid=color  lineattrs=(thickness=5);
		xaxis   label=  "Groups";
		yaxis   label= "AIC";
	run;
	quit;

%mend trajfit;


/***************************** EXPERIMENTAL **************************************/

%macro trajhclus(lib,numclus,res,rfbase); 
/* Macro for using output of trajloop with too many groups to do hierarchical clustering

LIB         Libname for trajloop output data sets
NUMCLUS     Number of clusters to output from PROC TREE.  Note that it is best to pick 
            an arbitrary value [2-10] for a first pass and then use the diagnostics to 
            pick a final value.  Run this macro twice to do this.
res         Resolution of final predicted values.
rfbase      This is the variable name for the risk factor used in TRAJLOOP. For example,
            if sbp is your risk factor, rfbase=sbp will reflect the sbp1-sbp100 variables
            in the WIDE data set that is used for clustering

*/

	%let numclus = &numclus;
	ods graphics on;
	proc cluster data=&lib..wide outtree=tree1 ccc pseudo method=ward plots=all;
		id group;
		var &rfbase.1-&rfbase.101;
	run;
	proc tree data=tree1 noprint ncl=&numclus out=treeout;
	run;
	data &lib..treeout; 
		set treeout;
		group = _NAME_ /1;;
		drop _NAME_;
	run;
	proc sort data=&lib..treeout; by group;
	proc sort data=&lib..clusout; by group;
	data &lib..finalcluster; 
		merge &lib..treeout (in=in1) &lib..clusout; 
		by group;
		if in1;
	run;
	proc print data=treeout;
	run;

	proc means data=&lib..finalcluster;
	run;

	ods output OneWayFreqs=owf;
	proc freq data=&lib..finalcluster;
		tables cluster;
	run;
	quit;
	proc print data=owf;
	run;
	data owf; set owf end=eof;
		if eof then call symput("nclus",cluster);
	run;
	%let nclus = &nclus;
	proc means data=&lib..finalcluster;
		var &gtime;
		output out=hdata min=minx max=maxx;
	run;
	data hdata; 
		set hdata;

		array ggroups (&nclus) cluster1-cluster&nclus;
		do i = 1 to &nclus;
		  	ggroups(i) = -1;
		end;
		&griskvar = .;
		do &gtime = minx to maxx by &res;
		  	output;
		end;
	run;
	data &lib..finalcluster; 
		set &lib..finalcluster;
		array ggroups (&nclus) cluster1-cluster&nclus;
		do i = 1 to &nclus;
		  	ggroups(i) = -1;
		end;
		ggroups(cluster) = 1;
	run;
	data cdata; 
		set &lib..finalcluster hdata;
	run;
	proc sort data=cdata; by cluster; run;

	%do ngrp = 1 %to &nclus;
		%clussmooth(method=&gmethod, dasn=cdata, idvar=&gid, wtvar=cluster&ngrp,
					timvar=&gtime, rskvar=&griskvar, dfvar=&gmaxdf,
					outdsn=hcout&ngrp, numvar=&ngrp);
		run;
	/* predicted values for fit dataset*/
	 	data xhpout&ngrp; 
			set hcout&ngrp;
	  		where cluster&ngrp eq -1;
	     	cluster = &ngrp;
	     	pred = pred&ngrp;
	 	run;

	/* Current data plot data set*/
		data plotdata&ngrp; 
			set hcout&ngrp;
		  	where cluster&ngrp eq 1;
		  	cluster=&ngrp;
		  	pred = pred&ngrp;
		  	keep cluster &gtime pred &gid;
	  	run;
	%end;

	data &lib..predout;
		set xhpout1-xhpout&nclus;
		by cluster;
	run;
	data &lib..allplot; 
		set plotdata1-plotdata&nclus; 
		by cluster;
	run;

%mend trajhclus;
