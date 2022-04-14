Discover processed pseudogene from long-read assemblies, 
given a reference assembly and a gene reference. 
Code is reorganized based on the scripts we wrote for manuscript
"Higher Rates of Processed Pseudogene Acquisition in Humans 
and Three Great Apes Revealed by Long-Read Assemblies".

## Usage

### python dependency and external tools

- python librarise: numpy

- need minimap2 with paftools.js for alignment and variant calling.
They are not called by the scripts.

Tested with python 3.7.6 on linux server, although earlier python 3.x should be fine as well. 

### gzip and piping

All input files can be gzipped or plain text.

For step 2 described below, please avoid pipes when feeding in the variant calling files (.var files).
The sample names correlated with each structural variant, and eventually the PPG candidates,
are drew from the file names.
If losing this information does not matter, then pipes are probably fine.

### workflow

Input: 

- long-read assembly/ies, fasta/q(.gz)

- a reference genome, fasta/q(.gz)

- gene reference, BED and fasta. When converting BED to fasta, use full length and do not remove introns.
Each sequence is expected to be named like 'id|type|symbol', e.g.'ENST00000376576.3|protein_coding|KIAA2013'.

Output:

- 2 intermediate files: prefix.PPG.paf and prefix.PPG.intermediate_results.
These two can be ignore.

- 1 final result file: prefix.PPG.final_results.
This contains the identified processed pseudogenes, equivalent to table S2 in the manuscript.

Steps:

1. Align long-read assembly to the reference genome without splice-awareness.
Then call variants to find the long structural variants (SVs).

2. Convert the SVs to fasta format, with location information preserved 
in the sequence names. 
Use `python findPPG_fileprep.py` (no argument or -h for usage).

3. Align the long structural variants to the gene reference, 
this time splice-aware.

4. Find processed pseudogenes candidates from the splice alignment.
Use `python findPPG_parsealn.py`.

5. Annotate the candidates.
Use `python findPPG_annotate.py`.

### toy demo

Download the demo input from github pre-release and extract with `tar zxf demo_input.tar.gz`.
This is contains some contigs from a human assembly (AK1), truncated hs38 and gencode's gene references. 

```
#step1
minimap2 -xasm5 -c --cs -z10000,200 -r5k -t32 ref.fa.gz contigs.fa.gz  > aln.paf
sort -k6,6 -k8,8n aln.paf > aln_sorted.paf
paftools call -q0 aln_sorted.paf 1>call.var 2>call.vst

#step2
python findPPG_fileprep.py -o call.tsv call.var 1> SV.fa

#step3
minimap2 -c --cs -x splice -f10000 -N100 -p0.1 -t32 refgene.fa.gz SV.fa > aln_SV.paf

#step4
python findPPG_parsealn.py -o PPGFparse -t32 refgene.bed.gz ref.fa.gz aln_SV.paf

#step5
python findPPG_annotate.py -o PPGFannotate -t32 --assembly contigs.fa.gz \
ref.fa.gz \
refgene.bed.gz \
call.tsv \
PPGFparse.PPG.paf

# the final result would be in PPGFannotate.PPG.final_results
```



