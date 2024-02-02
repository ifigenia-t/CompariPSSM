# CompariPSSM

CompariPSSM quantifies the similarity between short linear motifs (SLiMs) using sliding window PSSM-PSSM (Position Specific Scoring Matrix) comparison and scores PSSM similarity using a randomisation-based probabilistic framework. CompariPSSM can be used for peptide classification against lists of validated motif classes, peptide clustering to group functionally related conserved disordered regions, and benchmarking of experimental motif discovery methods. The CompariPSSM pipeline can also screen motif binding determinant models against protein alignments using PSSM-PSSM comparison thereby including residue conservation, a strong discriminator of motif functionality, directly into in silico motif discovery.

## Table of contents
* [Description](#PSSM-PSSM-Comparison-Approach)
* [CompariPSSM Web server](#CompariPSSM-Web-server)
* [Software prerequisites](#software-prerequisites)
* [Installation methods](#install)
* [Usage](#usage)
* [Input files](#input-files)
* [Output files](#output-files)
* [License](#license)
* [Reference](#reference)

## PSSM-PSSM Comparison Approach
The CompariPSSM framework is a pipeline for comparing motif binding determinants encoded as PSSMs. The framework integrates information on PSSM column similarity, scored using Pearson's Correlation Coefficient, and the PSSM column importance, scored using the Gini Coefficient, to align the PSSMs and calculate a similarity score for the best alignment. The CompariPSSM pipeline includes two major steps (Figure 1): (i) the Background Probability Pipeline that calculates the likelihood of pairwise PSSM column Importance-Weighted Similarity (IWS) score using a randomisation-based probabilistic framework, and (ii) the Comparison Pipeline which aligns PSSM pairs using a sliding window approach to define the optimal PSSM-PSSM alignment and applies information from the Background Probability Pipeline to calculate the likelihood of such an alignment occurring by chance. An additional Importance-weighted Dissimilarity (IWD) score that identifies important positions in the query PSSM that are not present in the corresponding position of the comparison PSSM is also calculated in parallel to the IWS score calculation, to discriminate between similar but non-identical motif binding determinants. CompariPSSM takes as input a query PSSM and a set of one or more comparison PSSMs and outputs the aligned best-matching comparison PSSM to the query PSSM with a statistical measure of the likelihood of that match happening by chance.

<p align="center">
  <img src="./docs/img/comparipssm_pipeline_figure_1.png" width="90%" height="100%" title="CompariPSSM Pipeline">
</p>

## CompariPSSM Web server
The CompariPSSM pipeline and interactive visualisations have been made available as a web server at https://slim.icr.ac.uk/projects/comparipssm. The CompariPSSM server has numerous input options: (i) input PSSMs, which can be copied and pasted directly in or upload a PSSM in a Tab-Delimited Table or JSON format; (ii) sets of aligned or unaligned peptides; and (iii) protein regions defined UniProt accession, gene name or protein name, and region and start and stop offsets. The query PSSM can be compared against a user-defined input PSSM or the Eukaryotic Linear Motif (ELM) resource-curated PSSM dataset. The output is the best match to the query PSSM, the comparison similarity (IWSsigwin) and dissimilarity score (max IWD score), and logos visualising the PSSM-PSSM comparison. If there is more than one significant match, additional hits are shown in a tabular form.

## Software prerequisites



## Install


## Usage

## Input files

## Output files

## License
This source code is licensed under the MIT license found in the `LICENSE` file in the root directory of this source tree.

## Reference
If you find the pipeline useful in your research, we ask that you cite our paper:
```
In preparation: 
"CompariPSSM: a PSSM-PSSM comparison tool for motif specificity determinant analysis" I. Tsitsa, I. Krystkowiak, N.E. Davey (2024) PMID:TBD
```