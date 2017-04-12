Carmichael et al. (submitted)

Code used for processing, analysis and visualization in Carmichael,
Gmaz, van der Meer (submitted),"Gamma oscillations in the rat ventral
striatum originate in the piriform cortex" ([preprint](http://www.biorxiv.org/content/early/2017/04/12/127126)).

Makes extensive use of the
[vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab/);
please use
[this commit](https://github.com/vandermeerlab/vandermeerlab/commit/245431a4e0e1f6f344dcc27cd82c561ba9664aef)
if you want to be sure you are using the same code that generated the
results in the paper.

To obtain the data, e-mail mvdm at dartmouth dot edu to get access to
the lab server. Then, to reproduce the results in the paper, run the
following:

LFP mapping: To load the data and generate intermediate files
(all_data_pre, all_data_post, all_data_task) run
[AMPX_LFP_MASTER.m](https://github.com/vandermeerlab/papers/blob/master/Carmichael_etal_submitted/Naris_Paper/Naris/AMPX_LFP_MASTER.m),
one cell at a time. Then, to generate Figures 1-6, run
[AMPX_figures_MASTER.m](https://github.com/vandermeerlab/papers/blob/master/Carmichael_etal_submitted/Naris_Paper/Naris/AMPX_figures_MASTER.m);
this requires the intermediate files generated in the previous step.

Naris occlusion: To load data, generate Figure 7, and get statistics
run
[Naris_master.m](https://github.com/vandermeerlab/papers/blob/master/Carmichael_etal_submitted/Naris_Paper/Naris/Naris_MASTER.m). This
can be done independently from the LFP mapping workflow.

Toolboxes used in this analysis:
 - FieldTrip - Used for LFP mapping experiment (Oostenveld, R., Fries,
   P., Maris, E., Schoffelen, JM (2011) FieldTrip: Open Source
   Software for Advanced Analysis of MEG, EEG, and Invasive
   Electrophysiological Data. Computational Intelligence and
   Neuroscience Volume 2011 (2011), Article ID 156869,
   doi:10.1155/2011/156869; http://www.ru.nl/neuroimaging/fieldtrip)

 - Chronux - used for part of the CSD analysis (http://chronux.org/)
   see reference material in: "Observed Brain Dynamics" by Partha
   Mitra and Hemant Bokil, Oxford University Press, New York, 2008
 
 - Circular Statistics Toolbox: P. Berens, CircStat: A Matlab Toolbox
   for Circular Statistics, Journal of Statistical Software, Volume
   31, Issue 10, 2009 http://www.jstatsoft.org/v31/i10

We used MATLAB R2014b running on 64-bit Windows 7.

Wiki-based tutorials introducing the codebase architecture, and
step-by-step explanations of specific analyses, are
[here](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:course-w16).
