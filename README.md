#### PPGfinder

Discover processed pseudogene from long structural variations.

#### usage

```sh 
# see help message of each script for more. 
python findPPG.py [OPTIONS] gencode.bed.gz ref.fa splice_alignment.paf   # append gene structure info to paf lines (PPGpaf); alignment is longSVs to a gene reference.
python findPPG_annotate.py [OPTIONS] -f dir_fasta -l dir_fasta_flk -p dir_PPGpaf -o output_name --nb_exons 2  # recover TSD, polyA and etc info and dump into a tsv. This script is shared for record purposes.
```
