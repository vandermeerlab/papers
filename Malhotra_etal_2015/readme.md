# Malhotra et al. (2015)

MATLAB and R code to generate the results in Malhotra, Cross, Zhang
and van der Meer (2015),
[European Journal of Neuroscience](http://onlinelibrary.wiley.com/doi/10.1111/ejn.13069/full) ([preprint](http://www.vandermeerlab.org/SM_MvdM_gamma_accepted.pdf))

We believe sharing code to promote reproducibility, collaboration, and openness is perferable over not sharing code "because it's messy". we ask for your understanding, and feel free to get in touch if anything is unclear.

Raw/pre-processed data files and metadata are available on request from the authors, contact mvdm at dartmouth dot edu.

Makes extensive use of the
[vandermeerlab codebase](https://github.com/vandermeerlab/vandermeerlab);
please use
[this release](https://github.com/vandermeerlab/vandermeerlab/releases/tag/v1.1)
if you want to be sure you are using the same code that generated the
results in the paper.

You will also need the [FieldTrip toolbox](http://www.fieldtriptoolbox.org/); we used [this commit](https://github.com/fieldtrip/fieldtrip/commit/a93aa21f4f65f933da5254f20265f8b8489668fe). 

To run this code, use the following example to set your path:

```
restoredefaultpath;
addpath(genpath('\github\vandermeerlab\code-matlab\toolboxes\MClust3.5\'));
addpath(genpath('\github\papers\Malhotra_etal_2015\shared'));
cd('\github\fieldtrip');
ft_defaults;
```

All R scripts require the file `ALL_trialinfo_fixedbaseline.mat` which
can be generated with the `BatchTrialInfo_create.m` script.

Please keep in mind that our lab codebase is always evolving and the
current version may have much better implementations than those in the
release snapshot.

Tutorial-style introductions to the data formats, data types, loaders,
and analyses used here can be found on our
[lab wiki](https://rcweb.dartmouth.edu/~mvdm/wiki).

