# splicedice
Splice DICE (Divergent Interval Co-Exclusion): A tool for RNA splicing quantification and analysis.

Currently contains the signature subprogram, which can perform differential splicing analysis on tables of percent-spliced values, and create a splicing signature for a particular sample group.

# Installation
All programs run as basic python programs. The current requirements are:
*python 3.7+
*numpy
*scipy
*matplotlib

# Usage
As input for the signature analysis, spliceDICE requires a table of percent-spliced values, previously output from spliceDICE or MESA, and a manifest (a tab-separated file with lists of samples and their group labels).

| sample 1 | group A
| -------- | --------
| sample 1 | group A
| sample 2 | group A
| sample 3 | group B
| sample 4 | group B
| sample 5 | control
| sample 6 | control

The recommended workflow is to first test for differential splicing, then use that set of significant splice intervals to generate a splicing signature. New samples can be queried against the signature to determine if they are a statistically significant match.

python ~/splicedice/code/signature.py compare -p project.ps.tsv -m manifest.tsv -o project
python ~/splicedice/code/signature.py fit_beta -s project.sig.tsv -p project.ps.tsv -m manifest.tsv -o project
python ~/splicedice/code/signature.py query -b project.beta.tsv -p new_samples.ps.tsv -o new_samples
python ~/splicedice/code/plot.py -q new_samples.pvals.tsv -m new_manifest.tsv
