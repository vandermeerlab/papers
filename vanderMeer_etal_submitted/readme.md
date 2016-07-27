# van der Meer et al. (submitted) Code used for processing, analysis
and visualization in van der Meer, Carey, Tanaka (submitted), "Towards
the principled decoding of internally generated sequences in the
hippocampus"

Makes extensive use of the
[vandermeerlab codebase](https://github.com/mvdm/vandermeerlab);
please use
[this commit](TBC)
if you want to be sure you are using the same code that generated the
results in the paper.

Once you have checked out the above code, set up your MATLAB path as follows:

```
restoredefaultpath; % start with clean slate
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared'));
addpath('');
```

To obtain the data, e-mail mvdm at dartmouth dot edu to get access to
the lab server. Then, to reproduce the results in the paper, run the
following:

- Step 1
- Step 2

Miscellaneous:

We used MATLAB R2014b running on 64-bit Windows 7.

Wiki-based tutorials introducing the codebase architecture, and
step-by-step explanations of specific analyses, are
[here](http://ctnsrv.uwaterloo.ca/vandermeerlab/doku.php?id=analysis:course-w16)
