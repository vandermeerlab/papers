# Carey et al. (submitted)

Code used for processing, analysis and visualization in Carey, Tanaka
and van der Meer, "Reward revaluation biases hippocampal replay
content away from the preferred outcome" (submitted)

Makes extensive use of the
[vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab);
please use
[this commit](https://github.com/vandermeerlab/vandermeerlab/commit/ad0bbd4d01726a436b36671c0a8b2db81476e946)
if you want to be sure you are using the same code that generated the
results in the paper.

Once you have checked out the above code, set up your MATLAB path as follows:

```
restoredefaultpath; % start with clean slate
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('\GitHub\vandermeerlab\code-matlab\tasks\Alyssa_Tmaze'));
addpath(genpath('\GitHub\papers\Carey_etal_submitted'));
```

The data is available from [datalad](http://datasets.datalad.org/?dir=/workshops/mind-2017/MotivationalT) or
from the lab server (email mvdm at dartmouth dot edu for credentials). Next, edit the 
`getTmazeDataPath.m` function to include the location where you put the data. 

Then, to reproduce the results in the paper, run the following:

- For analysis of **behavior** (Figure 1), run [Behavior_GenData.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/behavior/Behavior_GenData.m)
  followed by [PLOT_behaviorNewStyle.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/behavior/PLOT_behaviorNewStyle.m).
- For the main **sequenceless decoding log odds analysis** (Figures 2 and 3), place the [SWR candidate event files](https://github.com/vandermeerlab/papers/tree/master/Carey_etal_submitted/SWRcandidates) in their corresponding data folders.
  (or, use [MASTER_Generate_Tmaze_Candidates.m](https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/tasks/Alyssa_Tmaze/MASTER_Generate_Tmaze_Candidates.m) to generate your own.)
- Then, run [ALL_Generate_DecSeqCombined.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/decoding_noSeq/ALL_Generate_DecSeqCombined.m), which is a batch script to generate decoding data for each session.
- Once that finishes, [PLOT_DecSeqCombinedShuf.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/decoding_noSeq/PLOT_DecSeqCombinedShuf.m) will output the statistics and figures.
- [RESUME HERE] For the main **decoded sequence analysis** (Figure 5) decode each session using [ALL_GenerateDecSeq.m](https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/tasks/Alyssa_Tmaze/decoding/ALL_Generate_DecSeq.m).
- Collect data across sessions with [ALL_Collect_DecSeq_eligibleShuf.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/decoding_Seq/ALL_Collect_DecSeq_eligibleShuf.m), and
  generate the figures with [ALL_Plot_DecSeq_eligible.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/decoding_Seq/ALL_Plot_DecSeq_eligible.m). For these last two
  steps, the default parameters generate the results in the main text. For the supplementary figures, modify the parameters accordingly.

Other items:

- To generate the simulated data in Figure 6b, run
  [MASTER_scenarios.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/simulations/MASTER_scenarios.m).
- For the decoding accuracy null hypothesis (Figure 6a), run
  [MASTER_xvalDecErr.m](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/nullHypothesis/MASTER_xvalDecErr.m).
- For the tunning curve correlations between trials (Figure S2c-d) run [GENERATE_singleLapTCs](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/TCcorrelations/GENERATE_singleLapTCs.m) followed by [COLLECT_singleLapTCs](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/TCcorrelations/COLLECT_singleLapTCs.m).
- The regression models linking behavior and SWR content are [here](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/regression/).
- Example sequence plots (Figure 4) are generated with [this](https://github.com/vandermeerlab/papers/blob/master/Carey_etal_submitted/examples/GenerateExample.m).
- To run the suite of tests for sequence detection: `runtests('TEST_DecSeqDetectZ')`.

We used MATLAB R2017a running on 64-bit Windows 7.

Wiki-based tutorials introducing the codebase architecture, data and
metadata formats, and step-by-step explanations of specific analyses,
are
[here](http://discovery.dartmouth.edu/~mvdm/wiki/doku.php?id=analysis:nsb2018). We
also have a [Jupyter notebook](http://nbviewer.jupyter.org/github/summer-mind/mind_2017/blob/master/Tutorials/SpikeDecoding/spike_decoding_matlab.ipynb) with example analyses using this data set.

The data set is the same as that described in van der Meer, Carey and
Tanaka (2017) Optimizing for generalization in the decoding of
internally generated activity in the hippocampus, _Hippocampus_
27(5):580-595
([link](http://onlinelibrary.wiley.com/doi/10.1002/hipo.22714/full)).
