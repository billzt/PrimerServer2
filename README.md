# PrimerServer2
PrimerServer2: a high-throughput primer design and specificity-checking platform

![screenshot]( https://raw.githubusercontent.com/billzt/figure/master/PrimerServer-UI-typeA.png ) 

## Description
PrimerServer was proposed to design genome-wide specific PCR primers. It uses candidate primers produced by Primer3, uses BLAST and nucleotide thermodynamics to search for possible amplicons and filters out specific primers for each site. By using multiple threads, it runs very fast, ~0.4s per site in our case study for more than 10000 sites.

This repository is based on Python3 and acts as the successor of legacy [PrimerServer](https://github.com/billzt/PrimerServer).

## External Dependencies
**Add these two softwares to your system PATH**
* [Samtools](https://www.htslib.org/) (>=1.9).
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (>=2.2.18)

## Install

Don't use Python 3.9 or above since the [primer3-py](https://pypi.org/project/primer3-py/) module hasn't supported Python 3.9 yet.
```
conda create -n primer python=3.8
conda activate primer
```

### Via PIP (release only)
```
$ pip3 install primerserver2
```

### Via Github
```
$ git clone https://github.com/billzt/PrimerServer2.git
$ cd PrimerServer2
$ python3 setup.py install
```

## Run testing commands
```
** (if installed from pip,) tests/query_design_multiple and tests/example.fa can be obtained from this github repository.

** full mode: design primers and check specificity
$ primertool full tests/query_design_multiple tests/example.fa

** design mode: design primers only
$ primertool design tests/query_design_multiple tests/example.fa

** check mode: check specificity only
$ primertool check tests/query_check_multiple tests/example.fa

```

## Need to run the Web UI?
Please refer to the [wiki](https://github.com/billzt/PrimerServer2/wiki).

## Warning: About the reference genome
If you use reference genomes with many unplaced scaffolds, be caution since such scaffolds with great homology with main chromosomes might influence your results.
If possible, delete (some or all of ) these unplaced scaffolds.
For the human genome, we recommend the [no_alt_analysis_set](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/), which has all the PAR regions marked with N, to be used.

## Comparison of the CLI and Web version
|        |    CLI    |    Web UI    |
|    ----    |    ----    |    ----    |
|    Design primers    |     :heavy_check_mark:    |     :heavy_check_mark:    |
|    Checking specificity    |     :heavy_check_mark:    |     :heavy_check_mark:    |
|    Progress monitor    |     :heavy_check_mark:    |     :heavy_check_mark:    |
|    Number of tasks    |    High    |    Low    |
|    Alternative isoforms    |     :heavy_check_mark:    |     :x:    |
|    Exon-exon junction    |     :heavy_check_mark:    |     :x:    |
|    Pick internal oligos    |     :heavy_check_mark:    |     :x:    |
|    Custom Tm temperature    |     :heavy_check_mark:    |     :x:    |
|    Custom max amplicons    |     :heavy_check_mark:    |     :x:    |
|    Visualization    |     :x:    |     :heavy_check_mark:    |




