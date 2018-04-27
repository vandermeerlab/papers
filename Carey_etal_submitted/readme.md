# Carey et al. (submitted)

Code used for processing, analysis and visualization in Carey, Tanaka
and van der Meer, "Reward revaluation biases hippocampal sequence
content away from the preferred outcome" (submitted)

Makes extensive use of the
[vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab);
please use
[this commit](https://github.com/vandermeerlab/vandermeerlab/commit/08395c4648c986283f6b4b882c4b77708fdfb306)
if you want to be sure you are using the same code that generated the
results in the paper.

Once you have checked out the above code, set up your MATLAB path as follows:

```
restoredefaultpath; % start with clean slate
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('\GitHub\vandermeerlab\code-matlab\tasks\Alyssa_Tmaze'));
addpath(genpath('\GitHub\papers\Carey_etal_submitted'));
```

To obtain the data, e-mail mvdm at dartmouth dot edu to get access to
the lab server, or get them from
[datalad](http://datasets.datalad.org/?dir=/workshops/mind-2017/MotivationalT). Then,
to reproduce the results in the paper, run the following:

- For analysis of **behavior** (Figure 1), run [Behavior_GenData.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/behavior/Behavior_GenData.m)
  followed by [Behavior_MultiPlotData.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/behavior/Behavior_MultiPlotData.m) (for making the figure) and
  [BehavChi.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/behavior/BehavChi.m) (for generating the stats).
- For the main **left vs. right decoding log odds analysis** (Figures XX), first generate SWR candidate events using [MASTER_Generate_Tmaze_Candidates.m](https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/tasks/Alyssa_Tmaze/MASTER_Generate_Tmaze_Candidates.m). For
  the results in the main text, use the `amSWR` detector
  (default). For the results in Figure S4, use the following settings:

```
gen.SWRmethod = 'HT';
gen.MUAmethod = 'none';
gen.ThetaThreshold = [];
```

(remember to specify a different suffix, such as `suffix =
  '_HT1-3_noTheta';` to distinguish the resulting candidate files.)
- Then, run [ALL_Generate_DecSeqCombined.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/decoding_noSeq/ALL_Generate_DecSeqCombined.m), which is a batch script to generate decoding data for each session.
- Once that finishes, [PLOT_DecSeqCombined.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/decoding_noSeq/PLOT_DecSeqCombined.m) will output the statistics and figures.
- For the main **decoded sequence analysis** (Figure 2e, Figure 3a-b,
  supplementary figures), generate SWR candidate events using
  [MASTER_Generate_Tmaze_Candidates.m](https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/tasks/Alyssa_Tmaze/MASTER_Generate_Tmaze_Candidates.m) if you haven't already done so.
- Then, decode each session using [ALL_GenerateDecSeq.m](https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/tasks/Alyssa_Tmaze/decoding/ALL_Generate_DecSeq.m).
- Collect data across sessions with [ALL_Collect_DecSeq.m](https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/tasks/Alyssa_Tmaze/decoding/ALL_Collect_DecSeq.m), and
  generate the figures with [ALL_Plot_DecSeq.m](https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/tasks/Alyssa_Tmaze/decoding/ALL_Plot_DecSeq.m). For these last two
  steps, the default parameters generate the results in the main
  text. For the supplementary figures, modify the parameters accordingly.

Other items:

- To generate the simulated data in Figure 3d, run
  [MASTER_scenarios.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/simulations/MASTER_scenarios.m).
- For the decoding accuracy null hypothesis (Figure 3c), run
  [MASTER_xvalDecErr.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/nullHypothesis/MASTER_xvalDecErr.m).
- To run the suite of tests for sequence detection:
  `runtests('TEST_DecSeqDetectZ')`.

We used MATLAB R2017a running on 64-bit Windows 7.

Wiki-based tutorials introducing the codebase architecture, data and
metadata formats, and step-by-step explanations of specific analyses,
are
[here](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:course-w16). We
also have a [Jupyter notebook](http://nbviewer.jupyter.org/github/summer-mind/mind_2017/blob/master/Tutorials/SpikeDecoding/spike_decoding_matlab.ipynb) with example analyses using this data set.

The data set is the same as that described in van der Meer, Carey and
Tanaka (2017) Optimizing for generalization in the decoding of
internally generated activity in the hippocampus, _Hippocampus_
27(5):580-595
([link](http://onlinelibrary.wiley.com/doi/10.1002/hipo.22714/full)).
