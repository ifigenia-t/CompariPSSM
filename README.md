# CompariPSSM

CompariPSSM quantifies the similarity between short linear motifs (SLiMs) using sliding window PSSM-PSSM (Position Specific Scoring Matrix) comparison and scores PSSM similarity using a randomisation-based probabilistic framework. CompariPSSM can be used for peptide classification against lists of validated motif classes, peptide clustering to group functionally related conserved disordered regions, and benchmarking of experimental motif discovery methods. The CompariPSSM pipeline can also screen motif binding determinant models against protein alignments using PSSM-PSSM comparison thereby including residue conservation, a strong discriminator of motif functionality, directly into in silico motif discovery.

## Table of contents
* [Description](#PSSM-PSSM-Comparison-Approach)
* [CompariPSSM Web server](#CompariPSSM-Web-server)
* [Usage](#usage)
* [Input files](#input-files)
* [Output files](#output-files)
* [License](#license)
* [Reference](#reference)

## PSSM-PSSM Comparison Approach
The CompariPSSM framework is a pipeline for comparing motif binding determinants encoded as PSSMs. The framework integrates information on PSSM column similarity, scored using Pearson's Correlation Coefficient, and the PSSM column importance, scored using the Gini Coefficient, to align the PSSMs and calculate a similarity score for the best alignment. The CompariPSSM pipeline includes two major steps (Figure 1): 
(i) the Background Probability Pipeline that calculates the likelihood of pairwise PSSM column Importance-Weighted Similarity (IWS) score using a randomisation-based probabilistic framework, and 
(ii) the Comparison Pipeline which aligns PSSM pairs using a sliding window approach to define the optimal PSSM-PSSM alignment and applies information from the Background Probability Pipeline to calculate the likelihood of such an alignment occurring by chance. An additional Importance-weighted Dissimilarity (IWD) score that identifies important positions in the query PSSM that are not present in the corresponding position of the comparison PSSM is also calculated in parallel to the IWS score calculation, to discriminate between similar but non-identical motif binding determinants. CompariPSSM takes as input a query PSSM and a set of one or more comparison PSSMs and outputs the aligned best-matching comparison PSSM to the query PSSM with a statistical measure of the likelihood of that match happening by chance.

<p align="center">
  <img src="./docs/img/comparipssm_pipeline_figure_1.png" width="90%" height="100%" title="CompariPSSM Pipeline">
</p>
Figure 1. Overview of the CompariPSSM framework.

## CompariPSSM Web server
The CompariPSSM pipeline and interactive visualisations have been made available as a web server at https://slim.icr.ac.uk/projects/comparipssm. The CompariPSSM server has numerous input options: (i) input PSSMs, which can be copied and pasted directly in or upload a PSSM in a Tab-Delimited Table or JSON format; (ii) sets of aligned or unaligned peptides; and (iii) protein regions defined UniProt accession, gene name or protein name, and region and start and stop offsets. The query PSSM can be compared against a user-defined input PSSM or the Eukaryotic Linear Motif (ELM) resource-curated PSSM dataset. The output is the best match to the query PSSM, the comparison similarity (IWSsigwin) and dissimilarity score (max IWD score), and logos visualising the PSSM-PSSM comparison. If there is more than one significant match, additional hits are shown in a tabular form.

## Usage
### CLI
Call CompariPSSM directly from the command line:

```shell
python comparipssm.py --query_pssm_file ./pssm_sets/query_pssm.json --compare_pssm_file ./pssm_sets/elm_pssm.json --output_file ./pssm_sets/test.out.json 
```
### Import Library
Use CompariPSSM as a library:

```python
import comparipssm

comparipssm_runner = comparipssm.CompariPSSM() 
query_pssm = json.loads(open('./pssm_sets/query_pssm.json').read())
compare_pssm = json.loads(open('./pssm_sets/elm_pssm.json').read())
comparipssm_runner.options["query_pssm"] = query_pssm
comparipssm_runner.options["compare_pssm"] = compare_pssm
comparipssm_runner.options["significance_cutoff"] = 0.0001
pssm_comparison_response = comparipssm_runner.run_pssm_comparison_analysis()
```
Parameters:
- `significance_cutoff`: The PSSM-PSSM comparison p-value significance cutoff.

## Input files
Position-specific scoring matrices (PSSMs) store models of motif binding determinants as preference scores for individual amino acids at each position within the motif. PSSMs are represented as a matrix, with L columns, where L is the length of the motif peptide and 20 rows, one for each of the standard amino acids. PSSM columns generally refer to the position of the motif and the rows to the amino acid. PSSM JSON files are a dictionary of PSSMs where the key is the name of the PSSM. Each PSSM should have standard amino acids as keys with lists of the PSSM values for each position (see _PSSM file example format bellow_).

- `query_pssm_file`: .JSON file containing the query PSSM.
- `compare_pssm_file`: .JSON file containing the compare PSSM or multiple compare PSSMs.
- `elm_pssm.json`: .JSON file containing the PSSMs created from motif instances from the [ELM database](http://elm.eu.org/ "ELM database"). For each ELM class, peptides were extracted and aligned using the ELM-defined class consensus. The PSSMs were calculated from peptide alignments using the [PSSMSearch](https://slim.icr.ac.uk/pssmsearch/) web application with default parameters and the frequency PSSM construction method [(Krystkowiak et al., 2018)](https://academic.oup.com/nar/article/46/W1/W235/5033155?login=false).


#### PSSM file example format
```
{
    {"PPxY": 
        "A": [-0.178, -0.924, -0.178, 0.474, -0.924, -0.924, 0.474, -0.924, 0.841, -0.178, -0.178], 
        "C": [-0.094, 1.707, -0.094, -0.094, -0.094, -0.094, -0.094, -0.094, -0.094, -0.094, -0.094], 
        "D": [0.403, -0.219, 0.819, -0.604, -0.604, -0.604, -0.604, -0.604, -0.219, 0.403, 1.355], 
        .
        .
        .
        "Y": [-0.137, 1.408, -0.137, -0.137, -0.137, -0.137, -0.137, 20.0, -0.137, -0.137, -0.137]}
}
```

## Output files
The pipeline will produce a file:
- `output_file`: .JSON file containing the results of the PSSM-PSSM comparison.

Every PSSM-PSSM comparison dictionary contains:

- `best`: The best result of the PSSM-PSSM comparison.
- `significant`: All the significant hits are reported along with the respective motif class, consensus and similarity p-value and dissimilarity score.

#### Output scores

- `raw_score`: The Importance Weighted Similarity score for PSSM - PSSM comparison.
- `p_score_adj`: Final PSSM-PSSM comparison probability. The IWS sig win probability is the final p-value for a PSSM-PSSM comparison, and it is defined as the corrected Importance Weighted Similarity window probability (IWSpwin) score. 
- `dissimilarity_scores_max`: Importance Weighted Dissimilarity (IWD) score, used to highlight the difference between motif-PSSM and the model PSSM. Calculated as the maximum normalised product of the mean absolute error and Gini coefficient across the PSSM. A lower score indicates that the identified hit is more similar to the PSSM model.
- `query_motif_re`: The consensus of the query motif.
- `compare_motif_re`: The consensus of the compare motif.

#### Output file example format
```
{
    "query_PPxY": 
        {"best": 
            {"raw_score": 0.49891657444188475, 
            "p_score": 3.3682128176079197e-10, 
            "p_score_adj": 3.7050340993687115e-09,  
            "dissimilarity_scores_max": 0.394781248366091,  
            "query_motif_re": "[CG]TPPPPY[NV][PT][LD]", "compare_motif_re": "GTPPPPY.PL", 
            "query_offset_start": 1, "query_offset_end": 11, 
            .
            .
            .
            "compare_offset_start": 0, "compare_offset_end": 10,   
            "query_motif": "query_PPxY", "compare_motif": "compare_LIG_WW_1_compare", 
            "query_name": "PPxY", "compare_name": "LIG_WW_1_compare", 
            "compare_pssm_score_scheme": "frequency", "query_pssm_score_scheme": "frequency"}, 

        "significant": 
            {"compare_LIG_WW_1_compare": 
                {"p_score": 3.3682128176079197e-10, 
                "p_score_adj": 3.7050340993687115e-09, 
                "dissimilarity_scores_max": 0.394781248366091, 
                "query_name": "PPxY", 
                "compare_name": "LIG_WW_1_compare", 
                "compare_motif_re": "GTPPPPY.PL", 
                "query_motif_re": "[CG]TPPPPY[NV][PT][LD]"}}
        }
}
```
## Useful code snippets
```python
# FaSTPACE peptide alignment:

peptides = ['TSPDGGTTFEHLWSSL', 'SPEVFQHIWDFLEQPI', 'CPVDPFEAQWAALENK', 'EPPLSQETFSDLWKLL', 'APELDPFEAQWAALEG']

def align_peptides(peptides):
  # Requires fastpace installation - pip install fastpace
  from fastpace import run_motif_discovery
  results = run_motif_discovery(peptides, refine=0)
  return results["alignment"]["aligned_sequences"].values()

# PSSM Frequency Construction:

def get_frequency_pssm(peptides):
    aas =  list("ACDEFGHIKLMNPQRSTVWY")

    # check peptides
    if not peptides:
        return "Missing peptides"

    peptides_lens = [len(x) for x in peptides]
    peptides_lens = list(set(peptides_lens))
    if len(peptides_lens) != 1:
        return "Peptides not aligned."

    pssm_len = len(peptides[0])

    # pssm columns
    aaColumn = {}
    for i in range(0, pssm_len):
        aaColumn[i] = []

    for peptide in peptides:
        for i in range(0,pssm_len):
            aaColumn[i].append(peptide[i])

    ## pssm frequency
    pssm = {}
    for aa in aas:
        if aa not in pssm: pssm[aa] = []

        for i in range(0, pssm_len):
            pssm[aa].append(float(aaColumn[i].count(aa))/(len(aaColumn[i])))
    return pssm

# Gini Coefficient Calculation:

def get_pssm_gini(pssm):
    pssm_columns = get_columns(pssm)
    gini_scores = {}
    for column in pssm_columns:
        gini_scores[column] = gini_coefficient(pssm_columns[column])
    return gini_scores

def get_columns(pssm,aas=[]):
    if len(aas) == 0: aas = list("CPQNTSGAVILMFYWHKRDE")
    pssm_columns = {}
    for i in range(0,len(pssm['A'])):
        pssm_columns[i] = []
        for aa in aas:
            pssm_columns[i].append(pssm[aa][i])
    return pssm_columns

def gini_coefficient(column):
    diffsum = 0
    try:
        for i, xi in enumerate(column[:-1], 1):
            diffsum += sum([abs(xi -y) for y in column[i:]])
        return  diffsum / (len(column)**2 * sum(column)/len(column))
    except:
        return 0

```

## License
This source code is licensed under the MIT license found in the `LICENSE` file in the root directory of this source tree.

## Reference
If you find the pipeline useful in your research, we ask that you cite our paper:
```
"CompariPSSM: a PSSM-PSSM comparison tool for motif specificity determinant analysis" I. Tsitsa, I. Krystkowiak, N.E. Davey (2024)
PMID:39471470 [DOI: 10.1093/bioinformatics/btae644](https://doi.org/10.1093/bioinformatics/btae644)
```
