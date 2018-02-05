# Condamage
Score post-mortem deamination patterns. Also score patterns conditional on
deamination at the most 5' position, and patterns conditional on deamination
at the most 3' position.

# Prerequisites
Condamage uses *htslib* to parse bam and indexed fasta files.  The plotting
script requires *python* (tested with versions *2.7.14* and *3.6.1*) and
*matplotlib* (tested with version *2.1.0*)

# Installation
Check out the git repository, then build with `make`.

# Usage
1) Ensure the reference assembly has been indexed with `samtools` to create an
`fai` file. Ie. run
```
samtools faidx ref.fasta
```

2) Score the post-mortem damage patterns in `file.bam`, that was aligned to the
reference assembly `ref.fasta`.
```
condamage file.bam ref.fasta > mismatches.txt
```

3) Plot the damage patterns.
```
plot_condamage.py -o mismatches.pdf mismatches.txt
```
