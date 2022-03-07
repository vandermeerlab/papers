# Catanese et al. (2016)

Code used for processing, analysis and visualization in Catanese,
Carmichael & van der Meer, "Low and high gamma oscillations deviate in
opposite directions from zero-phase synchrony in the limbic
corticostriatal loop" (2016)
[Journal of Neurophysiology](http://jn.physiology.org/content/early/2016/03/04/jn.00914.2015) ([preprint](http://www.vandermeerlab.org/JC_MvdM_gamma_accepted.pdf))

Makes extensive use of the
[vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab);
please use
[this release](https://github.com/mvdm/vandermeerlab/releases/tag/v1.1)
if you want to be sure you are using the same code that generated the
results in the paper.

You will also need the [FieldTrip toolbox](http://www.fieldtriptoolbox.org/); we used [this commit](https://github.com/fieldtrip/fieldtrip/commit/a93aa21f4f65f933da5254f20265f8b8489668fe). 

Once you have checked out the above code, set up your MATLAB path as follows:

```
restoredefaultpath; % start with clean slate
cd('\GitHub\fieldtrip'); % remember to replace paths with yours
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared'));
addpath('\GitHub\papers\Catanese_vanderMeer2016\utility');
ft_defaults; % warning can be ignored
```

To obtain the data, e-mail mvdm at dartmouth dot edu to get access to
the lab server. Then, to reproduce the results in the paper, run the
following:

- For all event-based analyses (Figures 4, 5, 6 and 9) first run the
  [event detection](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/master/MASTER_CollectGammaEvents.m);
  subsequent analyses require the resulting `ALL_evt` variable to
  exist in the workspace.
- Amplitude correlations: [MASTER_AmplCorr.m](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/master/MASTER_AmplCorr.m) followed by [PLOT_AmplCorr.m](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/plotting/PLOT_AmplCorr.m).
- Phase slopes:
  [MASTER_SpectralConnectivity.m](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/master/MASTER_SpectralConnectivity.m)
  followed by
  [PLOT_SpectralConnectivity.m](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/plotting/PLOT_SpectralConnectivity.m).
- Ensemble classification: [MASTER_ClassifyGammaEvents.m](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/master/MASTER_ClassifyGammaEvents.m) followed by [PLOT_Classify.m](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/plotting/PLOT_Classify.m).
- Analysis of phase lags and slopes for signal generator inputs:
  [MASTER_probe.m](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/master/MASTER_probe.m).
- For the epoch-based analyses (Figure 3), run [this script](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/master/MASTER_meanPSD_COH_Script_JC.m).

Other items:

- Phase slope index schematic (Figure 5a):
  [PLOT_PSIexample.m](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/plotting/PLOT_PSIexample.m).
- Code for checking the units in FieldTrip's output for the
  [phase angle](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/simulations/icoh_unit_check.m)
  and
  [phase slopes](https://github.com/mvdm/papers/blob/master/Catanese_etal_2016/simulations/psi_unit_check.m)
  on simulated data. 

We used MATLAB R2014b running on 64-bit Windows 7.

Wiki-based tutorials introducing the codebase architecture, and
step-by-step explanations of specific analyses, are
[here](https://rcweb.dartmouth.edu/~mvdm/wiki).
