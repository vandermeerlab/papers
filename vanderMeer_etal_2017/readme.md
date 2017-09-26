# van der Meer et al. (2017)

Code used for processing, analysis and visualization in van der Meer,
Carey, Tanaka (2017) "Optimizing for generalization in the decoding of
internally generated activity in the hippocampus", _Hippocampus_
27(5):580-595
([link](http://onlinelibrary.wiley.com/doi/10.1002/hipo.22714/full)).

Makes extensive use of the
[vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab);
please use
[this commit](https://github.com/vandermeerlab/vandermeerlab/commit/44a3547f059bfeb828cb3f5aaedb74e5e644f92d)
if you want to be sure you are using the same code that generated the
results in the paper.

Once you have checked out the above code, set up your MATLAB path
using [this script](https://github.com/mvdm/papers/blob/master/vanderMeer_etal_submitted/misc/MASTER_path.m).

To obtain the data, e-mail mvdm at dartmouth dot edu to get access to
the lab server, or get them from
[DataLad](http://datasets.datalad.org/?dir=/workshops/mind-2017/MotivationalT). Then,
to reproduce the results in the paper, run the following:

- For the parameter sweep results in Figures 4 and 5, first use [this
  script](https://github.com/mvdm/papers/blob/master/vanderMeer_etal_submitted/ParameterSweep/GENERATE_decErr_paramSweep.m) to generate a struct (`ALL_decErr`) containing decoding
  accuracy data across all sessions and parameter values for a given
  data split. Obtain and save this struct for each split you want to test.
- Plotting the results is done [here](https://github.com/mvdm/papers/blob/master/vanderMeer_etal_submitted/ParameterSweep/PLOT_decErr_paramSweep.m). Make sure that the filenames
  specified at the top of this script match what you saved them as in
  the previous step. This same script generates Figure 8.
- Figure 6 requires generating an output data struct for each
  parameter combination of interest using [this script](https://github.com/mvdm/papers/blob/master/vanderMeer_etal_submitted/NumberOfTrials/GENERATE_loocv_byLap.m) and then
  plotting the results with [this](https://github.com/mvdm/papers/blob/master/vanderMeer_etal_submitted/NumberOfTrials/PLOT_loocv_byLap_multi.m).
- Figure 7 follows the same structure: generate the output variables
  with [this](https://github.com/mvdm/papers/blob/master/vanderMeer_etal_submitted/TrialDistance/GENERATE_singleLapDist.m) and then [plot](https://github.com/mvdm/papers/blob/master/vanderMeer_etal_submitted/TrialDistance/PLOT_singleLapDist.m). As above, remember to update the
  filenames used.

Extras:

- To produce the cross-validation schematic in Figure 3, run [this](https://github.com/mvdm/papers/blob/master/vanderMeer_etal_submitted/misc/decSchematic.m).
- Counting the number of units recorded (Table 1) is done by
  [this script](https://github.com/vandermeerlab/vandermeerlab/blob/44a3547f059bfeb828cb3f5aaedb74e5e644f92d/code-matlab/tasks/Alyssa_Tmaze/paper/PAPER_Collect_nUnitsRecorded.m).

Miscellaneous:

We used MATLAB R2014b running on 64-bit Windows 7.

Wiki-based tutorials introducing the codebase architecture, and
step-by-step explanations of specific analyses, are
[here](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:course-w16). We
also have a
[Jupyter notebook](http://nbviewer.jupyter.org/github/summer-mind/mind_2017/blob/master/Tutorials/SpikeDecoding/spike_decoding_matlab.ipynb)
with example analyses using this data set.
