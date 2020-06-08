### PPGfinder
Discover processed pseudogene from long structural variations. 

usage: 

`python findPPG.py [OPTIONS] gencode.bed.gz ref.fa splice_alignment.paf`

options:

-a  Output all annotations if set. [False]

-d  Store both coverages for exons and introns. [False]

-t  Threads. [1]

-e  At least INT exons should have coverage for an entry to be considered a PPG candidate.[3]

-c  Whether the output will be written to the current dir. Affects file naming. [False]
