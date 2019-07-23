'''PrimerServer2: a high-throughput primer design and specificity-checking platform

Github: https://github.com/billzt/PrimerServer2

PrimerServer was proposed to design genome(or transcriptome)-wide specific PCR primers. It uses candidate primers 
produced by Primer3, uses BLAST and nucleotide thermodynamics to search for possible amplicons and filters out 
specific primers for each site. By using multiple CPUs, it runs very fast, ~0.4s per site in our case study for 
more than 10000 sites.

External Dependencies:

blastn and makeblastdb (from NCBI BLAST+)

samtools (>=v1.9)

'''