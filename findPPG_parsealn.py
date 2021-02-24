import os, sys, gzip, argparse, re, binascii
from math import ceil
from binascii import hexlify
from multiprocessing import Pool

import numpy as np

__VERSION__ = "v0.1.5"
__LAST_UPDATE__ = "Jun 28, 2020"

####################
#    file utils    #
####################

def is_gz(filename):
    "Tells if the given filename leads to invalid file (-1), gzipped file (0), not non-gzipped file(1)."
    try:
        file = open(filename, 'rb')
    except FileNotFoundError:
        print('[warning::is_gz] File ({0}) not found. Continue anyway.'.format(filename), file=sys.stderr)
        return -1
    if binascii.hexlify(file.read(2)) == b'1f8b':
        file.close()
        return 0
    else:
        file.close()
        return 1

def opener(filename):
    "Determines if the file is gzipped and deal with it."
    status = is_gz(filename)
    if status==-1:
        print('[error::myopen] File ({0}) not found. Abort.'.format(filename), file=sys.stderr)
        exit(1)
    elif status==0: gzipped = True
    elif status==1: gzipped = False
    if gzipped:
        with gzip.open(filename, 'rb') as file:
            for line in file:
                yield line.decode()
    else:
        with open(filename) as file:
            for line in file:
                yield line

def var2fa(stream, gzipped=False):
    """convert variant calling's .var file to fasta"""
    for line in stream:
        if gzipped: line = line.decode()
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

def vcf2fa(stream, gzipped=False):
    """convert sniffle's .vcf to fasta"""
    for line in stream:
        if gzipped: line = line.decode()
        if line[0]=='#': continue
        chrom, pos, ID, ref, alt = line[:5]
        refl = len(ref)
        altl = len(alt)
        if max(refl, altl)>50:
            if refl>altl:
                sys.stdout.write('>DEL_{0}:{1}|{2}\n{3}\n'.format(chrom, pos, ID, ref))
            else:
                sys.stdout.write('>INS_{0}:{1}|{2}\n{3}\n'.format(chrom, pos, ID, alt))
    sys.stdout.flush()
    return 0


##############
# other utils#
##############
def ori(seq, l=70):
    "for better fasta formating: wrap long lines"
    output = [seq[i*l:(i+1)*l] for i in range(ceil(len(seq)/l))]
    return '\n'.join(output)
    
def overlap(s1, s2):
    span1, span2 = s1[:2], s2[:2]  # ignore tailings
    return max(0, min(span1[1], span2[1]) - max(span1[0], span2[0]))

def group_overlaps(intervals, queries):
    """Intervals and queries shall not be swapped. (normally `intervals` is for 
    exon structures and `queries` for alignments from cigar.) intervals do not need 
    to be sorted. However, non-overlapping queries are assumed."""
    output = []
    intervals.sort()
    queries.sort()
    for interval in intervals:
        tmp = []
        for query in queries:
            tmp.append(overlap(interval, query))
        output.append(np.sum(tmp))
    return output

def cigar2queries(cigar, close_short_gaps=False):
    cigar = re.findall('[0-9]+[NMD]', cigar)
    queries = []
    offset = 0
    if close_short_gaps:
        cigar_modi = []
        for c in cigar:
            if c[-1]=='M': cigar_modi.append(c)
            elif c[-1] in 'ND': 
                if int(c[:-1])<close_short_gaps:
                    cigar_modi.append(c[:-1]+'M')
                else:
                    cigar_modi.append(c)
    for c in cigar:
        l = int(c[:-1])
        if c[-1] in 'ND': offset+=l
        elif c[-1]=='M':
            if len(queries)==0: 
                queries.append([offset, offset+l])
                offset+=l
                continue
            if queries[-1][-1]==offset:
                queries[-1][-1] = offset+l
            else:
                queries.append([offset, offset+l])
            offset+=l
    return queries

def merge_intervals(intervals):
    output = []
    intervals.sort()
    for i in intervals:
        if len(output)==0:output.append(i); continue
        if i[0]<=output[-1][-1]: output[-1][-1] = i[1]
        else: output.append(i)
    return output

def load_gencode(filename, with_structure=False, only_pc=False, select_types=[], select_inverse=False):
    gencode = {}
    file = opener(filename)  # hs38.gencode.v31.bed.gz
    for line in file:
        if only_pc and ('protein_coding' not in line): continue
        if len(select_types)!=0:
            if line.split('\t')[3].split('|')[1] not in select_types: continue
        line = line.strip().split('\t')
        chrom, start, end, name = line[0], int(line[1]), int(line[2]), line[3]
        if with_structure:
            lengths = [float(_) for _ in line[-2].split(',') if _!='']
            exon_starts = [float(_) for _ in line[-1].split(',') if _!='']
            if select_inverse:
                introns = []
                for i in range(len(exon_starts)-1):
                    introns.append([exon_starts[i]+lengths[i], exon_starts[i+1]])
                gencode[name] = [chrom, start, end, introns]
            else:
                exons = [[s, s+l] for s, l in zip(exon_starts, lengths)]
                gencode[name] = [chrom, start, end, exons]
        else:
            gencode[name] = [chrom, start, end]
    return gencode
            
def load_hs38_seq(filename):
    hs38 = {}
    file = opener(filename)  # hs38.fa
    name, l = '', ''
    for line in file:
        if line.startswith('>'):
            if name!='': hs38[name] = l
            name = line[1:].split(' ')[0]
            l = ''
        else:
            l+=line.strip()
    hs38[name] = l
    return hs38
    
def load_hs38_chrlength(filename):
    hs38_chrlen = {}
    file = opener(filename)  # hs38.fa
    name, l = '', 0
    for line in file:
        if line.startswith('>'):
            if name!='': hs38_chrlen[name] = l
            name = line[1:].split(' ')[0]
            l = 0
        else:
            l+=len(line)-1
    hs38_chrlen[name] = l
    return hs38_chrlen

##############
#  main func #
##############
def worker(bundle):
    f, file_out_name, qnames, workerID = bundle
    #print('[debug] worker{0}, {1}'.format(workerID, ', '.join(qnames[:10])), file=sys.stderr)
    if qnames!=None: qnames = set(qnames)
    file_out = open(file_out_name+'.tmp{0}'.format(workerID), 'w')
    file = opener(f)
    for line in file:
        line = line.strip().split('\t')
        if qnames!=None and line[0] not in qnames: 
            continue
        #print('[debug] process a line', file=sys.stderr)
        cigar = [_ for _ in line if 'cg:Z' in _][0][5:]
        if 'N' not in cigar: continue  # discard alignments without splicing
        cigar = '{0}N'.format(line[7])+cigar  # adjust offset (Target start on original strand (0-based))
        name_fragment, name_gene = line[0], line[5]
        intervals = np.array(gencode[name_gene][-1])
        queries = cigar2queries(cigar)
        cov = group_overlaps(intervals, queries)
        cov_relative = np.array([c/l for c, l in zip(cov, [e-s for s, e in gencode[name_gene][-1]])])
        if is_storeBothCov:
            intervals_inverse = np.array(gencode_inverse[name_gene][-1])
            covq = group_overlaps(intervals_inverse, queries)
            covq_relative = np.array([c/l for c, l in zip(covq, [e-s for s, e in intervals_inverse])])
        if is_outputall: condition = True
        else: 
            if is_storeBothCov: condition = np.sum(cov_relative>0.0)>=threshold_nb_exons and np.max(cov_relative)>0.7 and np.sum(covq_relative>0.5)<5  # intron restriction is loose
            else: condition = np.sum(cov_relative>0.0)>1 and np.max(cov_relative)>0.7
        if condition:
            #print(name_fragment, file=sys.stderr)
            name_gene_short = name_gene.split('|')[-1]
            # modify paf from gene coordinates to genomic coordinates
            transcript_chrom, transcript_start, transcript_end = gencode[name_gene][:3]
            line[5] = transcript_chrom
            line[6] = str(hs38_chrlen[line[5]])
            line[7] = str(int(line[7]) + transcript_start)
            line[8] = str(int(line[8]) + transcript_start)
            info = '|'.join(line[12:])
            if is_storeBothCov:
                line = line + [name_gene, str(transcript_start)+','+str(transcript_end),
                           ','.join([str(_) for _ in cov]),
                           ','.join(['%.2f'%_ for _ in cov_relative]),
                           ','.join([str(_) for _ in covq]),
                           ','.join(['%.2f'%_ for _ in covq_relative])]
            else:
                line = line + [name_gene, str(transcript_start)+','+str(transcript_end),
                           ','.join([str(_) for _ in cov]),
                           ','.join(['%.2f'%_ for _ in cov_relative])]
            file_out.write('\t'.join(line)+'\n')  # write to std
            file_out.flush()
#            print('[log] FOUND', file=sys.stderr)
        else:
#            print('[log] -', file=sys.stderr)
            pass
    file_out.close()


if __name__=="__main__":
    print("findPPG.py Version {0}, {1}".format(__VERSION__, __LAST_UPDATE__), file=sys.stderr)
    if len(sys.argv)<2:
        print("usage: python findPPG.py [OPTIONS] gencode.bed.gz ref.fa splice_alignment.paf", file=sys.stderr)
        print("       (All three files could be either gzipped or not.)", file=sys.stderr)
        print("options:", file=sys.stderr)
        print("   -a  Output all annotations if set. [False]", file=sys.stderr)
        print("   -d  Store both coverages for exons and introns. [False]", file=sys.stderr)
        print("   -t  Threads. [1]", file=sys.stderr)
        print("   -e  At least INT exons should have coverage for an entry to be considered a PPG candidate.[3]", file=sys.stderr)
        print("   -c  Whether the output will be written to the current dir. Affects file naming. [False]", file=sys.stderr)
        quit(0)

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', action="store_true", default=False)
    parser.add_argument('-d', action="store_true", default=False)
    parser.add_argument('-t', action="store", type=int, default=1)
    parser.add_argument('-e', action="store", type=int, default=3)  # least exon cov
    parser.add_argument('-c', action="store_true", default=False)
    parser.add_argument('gencode', action="store", type=str)
    parser.add_argument('ref', action="store", type=str)
    parser.add_argument('aln', action="store", type=str)
    args, unkown_args = parser.parse_known_args(sys.argv[1:])
    is_outputall = args.a
    is_storeBothCov = args.d
    is_writeInCurrentDir = args.c
    nb_cpu = args.t
    threshold_nb_exons = args.e
    filename_gencode = args.gencode
    filename_ref = args.ref
    filename_input = args.aln

    print('[log::main] requested nb_cpu=={0}'.format(nb_cpu), file=sys.stderr)
    print('[log::main] is_outputall:', is_outputall, file=sys.stderr)
    print('[log::main] loading references...', file=sys.stderr)
    gencode = load_gencode(filename_gencode, with_structure=True, only_pc=True)
    if is_storeBothCov: gencode_inverse = load_gencode(filename_gencode, with_structure=True, only_pc=True, select_inverse=True)  # load intron structures instead of exons
    hs38_chrlen = load_hs38_chrlength(filename_ref)
    print('[log::main] Got references.', file=sys.stderr)

    if nb_cpu==1:
        print('[log::main] single thread', file=sys.stderr)
        worker([filename_input, filename_input, None, 0])
    else:  # bad multip but ok
        print('[log::main] multi thread', file=sys.stderr)
        # probe
        total = 0
        qnames = []
        file = opener(filename_input)
        try:
            for line in file:
                qnames.append(line.split('\t')[0])
        except EOFError: 
            # This is to allow annotated a running alignment for a quick peek.
            print('[warning::main] EOF when preparing for multithread; continue anyway.', file=sys.stderr)
        file.close()

        # prepare
        qnames = sorted(list(set(qnames)))
        total = len(qnames)
        print('[log_multi::main] total {0} querie(s)'.format(total), file=sys.stderr)
        step = ceil(total/nb_cpu)
        print('[log_multi::main] each thread will have {0} querie(s) to process'.format(step), file=sys.stderr)
        qnames_bundled = [qnames[i*step:(i+1)*step] for i in range(nb_cpu)]
        if not is_writeInCurrentDir: 
            bundles = [[filename_input, filename_input, qnames_bundled[i], i] for i in range(nb_cpu)]
        else:
            bundles = [[filename_input, filename_input.split('/')[-1], qnames_bundled[i], i] for i in range(nb_cpu)]

        # run
        with Pool(nb_cpu) as p:
            p.map(worker, bundles)
    os.system('cat {0}.tmp* > {0}.PPGpaf; rm -f {0}.tmp*'.format(filename_input))
    quit(0)
