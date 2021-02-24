import os, sys, gzip
import numpy as np
from collections import Counter
from multiprocessing import Pool

############### NOTE ###################
# This is a collection of example commands and might not
#   suit for reuse. Modify /path/to/* lines and print/os.system lines.
# A work directory is structured as the following:
#  /path
#     |
#     --to
#       |
#       ----0_aln
#       ----1_SV
#       ----1_SV_fa
#       ----1_SV_flk
#       ----2_splicealn
#       ----3_annotation
########## end of NOTE #################


#### (1) contig to reference alignment
# // minimap2 with -xasm5 -c --cs -z10000,200 -r5k for humans, -xasm20 -c --cs -z10000,200 -r5k for great apes
# // sort the paf with `sort -k6,6 -k8,8n`
fs = [_ for _ in os.listdir('/path/to/0_aln/') if 'sorted' in _]
cmds = []
for f in fs:
    cmd = 'cd /path/to/0_aln; paftools call -q0 ../0_aln/{0} 1>{0}.var 2>{0}.vst'.format(f)
    cmds.append(cmd)
print(' &\n'.join(cmds))

#### (2) get SV sequences
def var2fa(stream):
    """convert variant calling's .var file to fasta"""
    for line in stream:
        if line[0]!='V': continue
        line = line.strip().split('\t')
        _1, chrom, start, end, _2, _3, ref, alt, queryname, q_start, q_end, strand = line
        if abs(len(ref)-len(alt))<50: continue  # not long enough
        if len(ref)>len(alt):
            newname = 'DEL_'+'_'.join([queryname, chrom+strand, start+'-'+end, q_start+'-'+q_end])
            seq = ref.upper()
        else:
            newname = 'INS_'+'_'.join([queryname, chrom+strand, start+'-'+end, q_start+'-'+q_end])
            seq = alt.upper()
        sys.stdout.write('>{0}\n{1}\n'.format(newname, seq))
    return 0
dir_base = '/path/to/'
fs = [_ for _ in os.listdir(dir_base+'1_SV/') if _.endswith('.var')]
for f in fs:
    os.system('python filePrep.py {0}1_SV/{1} > {0}1_SV_fa/{1}.fa'.format(f))


#### (3) splice alignment 
# (write to is /path/to/2_splicealn)
# // minimap2 -c --cs -x splice -f10000 -N100 -p0.1 -t nb_cpu -d hs38_proteinCodingGenes-splice.mmi hs38_proteinCodingGenes.fa

#### (4) findPPG_parsealn.py