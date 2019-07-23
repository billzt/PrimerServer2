# PrimerServer2
PrimerServer2: a high-throughput primer design and specificity-checking platform

# Description
PrimerServer was proposed to design genome-wide specific PCR primers. It uses candidate primers produced by Primer3, uses BLAST and nucleotide thermodynamics to search for possible amplicons and filters out specific primers for each site. By using multiple threads, it runs very fast, ~0.4s per site in our case study for more than 10000 sites.

This repository is based on Python3 and acts as the successor of legacy [PrimerServer](https://github.com/billzt/PrimerServer2).

# External Dependencies
**Add these two softwares to your system PATH**
* [Samtools](https://www.htslib.org/) (>=1.9).
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (>=2.2.18)

# Install
## Via PIP
```
    pip3 install primerserver2
```

## Via Github
```
    git clone https://github.com/billzt/PrimerServer2.git
    cd PrimerServer2
    python3 setup.py install
```

# Run testing commands
```
    # (if installed from pip,) tests/query_design_multiple and tests/example.fa are included in this github repository.

    # design primers and check specificity (the default mode)
    primertool tests/query_design_multiple tests/example.fa -o example_design_check.json -t example_design_check.tsv

    # check specificity only
    primertool tests/query_check_multiple tests/example.fa --only-specificity -o example_check.json -t example_check.tsv
```