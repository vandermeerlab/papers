# Butler et al. (submitted)
Code used for the head direction model fitting procedure in Butler et al. (submitted)

Makes extensive use of the
[vandermeerlab codebase](https://github.com/mvdm/vandermeerlab);
please use
[this commit](https://github.com/mvdm/vandermeerlab/commit/4977b3efd380fb81f29cb949195f1c7693050416)
if you want to be sure you are using the same code that generated the
results in the paper.

Once you have checked out the above code, set up your MATLAB path as follows:

```
restoredefaultpath; % start with clean slate
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('\GitHub\papers\Butler_etal'));
```

The workflow for generating the results is as follows:

- For the parameter recovery analysis, run MASTER_paramRecovery.m, followed by PLOT_paramRecovery.m

- For the model fits to HD cell data, run MASTER_hdfit.m, followed by PLOT_hdfit.m

Also:

- PLOT_examples.m to generate the insets in the model schematic
