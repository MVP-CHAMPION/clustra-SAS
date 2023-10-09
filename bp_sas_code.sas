/* Load the relevant macros */
%include "Trajectories Macros File V14R1.sas";
libname final "sas_results";

/* Import simulated data */
proc import	datafile = 'bp_data/simulated_data_27June2023.csv'
	dbms=csv
	out=trajdata
	replace;
 		getnames = yes;
   		delimiter=",";
run;
data trajdata; 
	set trajdata;
	rename group=true_group;
run;

/* checking dataset */
proc contents data=trajdata; run;
proc sort data=trajdata out=final.trajdata; by id time; run;
proc print data=trajdata (obs=20); run;

/* Macro for iterating through cluster size = n */
%macro doit(gnum=5);
	/* Title */
	title "Fitting &gnum. Clusters to the Data";

	/*location of data output*/
	libname traj&gnum.  "sas_results/&gnum.clusters";

	/*Setting up graphs for temporary output*/
	%graphset(gpath=sas_results/&gnum.clusters, 
			dpi=300, 
			style=splinecurve, 
			format=jpeg, 
			type=html);

	/* running with 2 trajectories, and response as systolic blood pressure */
	%trajsetup(dsn=final.trajdata, id=id, time=time, riskvar=response, ngroups=&gnum., 
				maxdf=30, ptrim=0, 
				seed=1256, steps=1, method=0, random=YES);

	ods exclude all; *options to suppress all the intermediate output from this process, while still retaining relevant output datasets;
	%trajloop(outlib=traj&gnum., iter=50, minchange=0.5, minsubs=0, /*specifying up to 50 iterations for consistency with 10 cluster approach */
				min_x=-365, max_x=730, by_x=50,
				min_y=110, max_y=180, by_y=10,
				showall=NO, showany=YES); *, AIC=YES);
	run;
	ods select all; *turn the results back on;

	data final.graphout_cl&gnum; set traj&gnum..graphout; run; /* save the dataset for George */
%mend;

/* macro for plotting the results following algorithms with cluster size = n specified 
		Note: depends on output from %doit(gnum=) macro*/
%macro plot_it(gnum=5);
	/* Title */
	title "Fitting &gnum. Clusters to the Data";

	/*location of data output*/
	libname traj&gnum.  "sas_results/&gnum.clusters";

	/* Generate the plots using the output data */
	ods html close; 
	/* Trajectories plot */
	ods html gpath="sas_results"
			style=splinecurve	image_dpi = 300;
	ods graphics on / reset width=10in height=6.8in IMAGEFMT=jpeg
			imagename="TrajPlotSAS_cl&gnum" ;
	proc sort data=traj&gnum..graphout; by group time pred;
	run;
	proc sgplot data=traj&gnum..graphout dattrmap=graphattr;
		title "Fitting &gnum. Clusters to the Data";
		series y=pred x=time / group=group attrid=color;
		xaxis  values=(-365 to 730 by 50) label=  "Time from Treatment Initiation (in days)";
		yaxis  values=(110 to 180 by 10) label= "Predicted Systolic Blood Pressure";
	run;
	quit;
	ods graphics off;
	ods html close;

	/* Silhouette plot */
	* Run the silhouette code to get the output data needed;
	/*%silhouette(lib=traj&gnum., sid=id, groupsn=&gnum, cgroup=group); *hg added specifications for outputting the plot data;*/
	/*run;*/

	* Extract the relevant silhouette statistics;
	data _null_; 
		set traj&gnum..silout; *output from %silhouette;
		call symputx("meansil_lab",round(meansil,0.001)); *grab the value for the mean silhouette from the output data;
		call symputx("medsil_lab",round(medsil,0.001)); *grab the value for the median silhouette from the output data;
	run;
	data plotdat; *grab the plot data that was output from the silhouette function;
		set traj&gnum..sil_plotdat; 
		label primegroup='group' xvar='observation #';
	run;

	/* Do the plotting */
	ods html gpath="sas_results" 
			style=splinecurve	image_dpi = 300;
	ods graphics on / width=10in height=6.8in IMAGEFMT=jpeg
			imagename="SilPlotSAS_cl&gnum";
	proc sgplot data=plotdat;
		title1 'Silhouette Plot';
		title2 h=1 "Mean=&meansil_lab (Black), Median=&medsil_lab (Red)"; /*changed printing output*/
		band x=xvar upper=silhouette lower=start /group=primegroup;
		refline &meansil_lab /axis=y lineattrs=(color=black);
		refline &medsil_lab /axis=y lineattrs=(color=red);
	run;
	ods graphics off;
%mend; 

*calculate the trajectories;
%doit(gnum=2); *2 clusters (1.5 min on VINCI SAS 9x);
%doit(gnum=5); *5 clusters (5.5 min on VINCI SAS 9x);
%doit(gnum=10); *10 clusters (45 min on VINCI SAS 9x);
run;

* generate the plots for the manuscript;
%plot_it(gnum=2);
%plot_it(gnum=5);
%plot_it(gnum=10);

* Save the final ID-grouping results for R comparisons;
data final.rms_cl2; set traj2.rms; run;
data final.rms_cl5; set traj5.rms; run;
data final.rms_cl10; set traj10.rms; run;

/* Look at the iterations and number of persons switching groups between iterations */
proc print data=traj2.tracking;
	title1 'Summary of changes per iteration';
	title2 "With 2 Clusters Specified";
	var Frequency Percent;
run;
proc print data=traj5.tracking;
	title1 'Summary of changes per iteration';
	title2 "With 5 Clusters Specified";
	var Frequency Percent;
run;
proc print data=traj10.tracking;
	title1 'Summary of changes per iteration';
	title2 "With 10 Clusters Specified";
	var Frequency Percent;
run;

/***********************************************************************************************************/
title 'Iteration Plots for Fitting 5 Clusters';

/* Iteration plots */
libname itertraj  "sas_results/iteration_plots";
%graphset(gpath=sas_results/iteration_plots, 
			dpi=300, style=splinecurve, format=jpeg,type=html);
/*ods html close; ods html;*/

%macro trajit(maxi=13);

%do iterate=0 %to &maxi;
/* running with 5 trajectories, and response as systolic blood pressure */
	data td&iterate.; set trajdata; keep id time response; run;
	proc sort data=td&iterate.; by id time; run;
	%trajsetup(dsn=td&iterate., id=id, time=time, riskvar=response, ngroups=5, maxdf=30, 
				ptrim=0, seed=1256, steps=1, method=0, random=YES);

	ods select none;
	%trajloop(outlib=itertraj,iter=&iterate.,minchange=0.5,minsubs=0,
				min_x=-365, max_x=730, by_x=50,
				min_y=110, max_y=180, by_y=10,
				showall=NO, showany=YES);
	run;
	data itertraj.graphout_iter&iterate.;
		set itertraj.graphout;
	run;
	ods select all;
%end;

%mend;

%macro trajit_plot(maxi=13);
%do iterate=0 %to &maxi;
	title 'Iteration Plots for Fitting 5 Clusters';
	ods html close;
	ods html gpath="sas_results/iteration_plots" 
	style=splinecurve image_dpi = 300; 
	ods graphics on / reset width=10in height=6.8in IMAGEFMT=jpeg
		imagename="TrajPlotSAS_cl5_iter&iterate.";
	proc sort data=itertraj.graphout_iter&iterate.; by group time pred;
	run;
	proc sgplot data=itertraj.graphout_iter&iterate. noautolegend dattrmap=graphattr;
		title2 "Iteration &iterate.";
		series y=pred x=time / group=group attrid=color;
		xaxis  values=(-365 to 730 by 50) label=  "Time from treatment initiation (in days)";
		yaxis  values=(110 to 180 by 10) label= "Predicted Systolic Blood Pressure";
	run;
	quit;
	ods html close;
%end;

%mend;

/* converged after 10 iterations using this seed, so look at the ten plots */
%trajit(maxi=10);

%trajit_plot(maxi=10);
