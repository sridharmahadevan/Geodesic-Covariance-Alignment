This repository contains code for  the ECML 2018 paper titled
"A Unified Framework for Domain Adaptation using Metric Learning on Manifolds".

The paper can be downloaded from https://arxiv.org/abs/1804.10834

The code is in MATLAB, and requires several toolboxes, specifically
the Parallel Distributed Computing Toolbox.

To run the code, ensure that all the folders are in MATLAB's path. 

In the GCA directory, the main routine is gca123_vs_coral.m

As a sample experiment, try

gca123_vs_coral('dslr', 'amazon', 20, 0.9, 0.2, 0.1)

The first argument is the source domain, the second argument is the
target domain.  The third argument is the number of trials.  The last
three arguments are hyperparameters that are described more fully in
the paper.

The result should look similar to the PDF plot in the repository. The results will vary depending on source and target domains, and hyper parameter values.

It should be easy to extend the code to other domain adaptation
problems, besides the Office dataset.


