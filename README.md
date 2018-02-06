# Condamage
Score post-mortem deamination patterns. Also score patterns conditional on
deamination at the most 5' position, and patterns conditional on deamination
at the most 3' position.

# Prerequisites
Condamage uses **htslib** to parse bam and indexed fasta files.  The plotting
script requires **python** (tested with versions **2.7.14** and **3.6.1**) and
**matplotlib** (tested with version **2.1.0**)

# Installation
Clone the git repository, then build with `make`.

# Usage
* Ensure the reference assembly has been indexed with `samtools` to create an
`fai` file. Ie. run
```
samtools faidx ref.fasta
```

* Score the post-mortem damage patterns in `file.bam`, that was aligned to the
reference assembly `ref.fasta`.
```
condamage file.bam ref.fasta > mismatches.txt
```

* Plot the damage patterns (double stranded library).
```
plot_condamage.py -o mismatches.pdf mismatches.txt
```

* Plot the damage patterns (single stranded library).
```
plot_condamage.py -s -o mismatches.pdf mismatches.txt
```
