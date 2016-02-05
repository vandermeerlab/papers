# Catanese et al. (submitted)
Code used for processing, analysis and visualization in Catanese, Carmichael et
al. (submitted)

Makes extensive use of the
[vandermeerlab codebase](https://github.com/mvdm/vandermeerlab);
please use
[this release](https://github.com/mvdm/vandermeerlab/releases/tag/v1.1)
if you want to be sure you are using the same code that generated the
results in the paper.

You will also need the [FieldTrip toolbox](http://www.fieldtriptoolbox.org/); we used [this commit](https://github.com/fieldtrip/fieldtrip/commit/a93aa21f4f65f933da5254f20265f8b8489668fe). 

Once you have checked out the above code, set up your MATLAB path as follows:

```
restoredefaultpath;
addpath(genpath('\GitHub\vandermeerlab\code-matlab\shared'));
cd('\GitHub\fieldtrip');
ft_defaults;
```

We used MATLAB R2014b running on 64-bit Windows 7.

The original data sets, as well as wiki-based tutorials introducing
the codebase architecture, along with various hands-on examples using
real data sets, are available on request.
