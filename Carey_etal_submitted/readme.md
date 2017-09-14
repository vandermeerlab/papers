# Carey et al. (submitted)

Code used for processing, analysis and visualization in Carey, Tanaka
and van der Meer, "Motivational shifts bias hippocampal sequence
content away from the preferred outcome" (submitted, [preprint]())

Makes extensive use of the
[vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab);
please use
[this commit](https://github.com/vandermeerlab/vandermeerlab/commit/932406314574b3e2d295aeb9aed303e66d94601b)
if you want to be sure you are using the same code that generated the
results in the paper.

Once you have checked out the above code, set up your MATLAB path as follows:

```
restoredefaultpath; % start with clean slate
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('\GitHub\vandermeerlab\code-matlab\tasks\Alyssa_Tmaze'));
```

To obtain the data, e-mail mvdm at dartmouth dot edu to get access to
the lab server, or get them from [datalad](http://datasets.datalad.org/?dir=/workshops/mind-2017/MotivationalT). Then, to reproduce the
results in the paper, run the following:

- For behavior (Figure 1), run [PAPER_Collect_Behavior.m]().
- For the main decoded sequence analysis (Figure 2e, Figure 3a-b,
  supplementary figures), first generate SWR candidate events using
  [MASTER_Generate_Tmaze_Candidates.m](). For the results in the main
  text, use the `amSWR` detector (default). For the results in Figure
  S4, use the following settings:
  ```
  gen.SWRmethod = 'HT'; 
  gen.MUAmethod = 'none';
  gen.ThetaThreshold = [];
  ```
  (remember to select a different suffix, such as `suffix =
  '_HT1-3_noTheta';` to distinguish the resulting candidate files.)
- Then, decode each session using [ALL_GenerateDecSeq.m]().
- Collect data across sessions with [ALL_Collect_DecSeq.m](), and
  generate the figures with [ALL_Plot_DecSeq.m](). For these last two
  steps, the default parameters generate the results in the main
  text. For the supplementary figures, modify the parameters accordingly.

Other items:

- To generate the simulated data in Figure 3d, run
  [MASTER_scenarios.m]().
- For the decoding accuracy null hypothesis (Figure 3c), run
  [MASTER_xvalDecErr.m]().
- To run the suite of tests for sequence detection:
  `runtests('TEST_DecSeqDetectZ')`.

We used MATLAB R2017a running on 64-bit Windows 7.

Wiki-based tutorials introducing the codebase architecture, and
step-by-step explanations of specific analyses, are
[here](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:course-w16).
