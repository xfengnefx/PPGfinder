import os, sys, gzip, argparse, re, binascii
from math import ceil
from binascii import hexlify
from multiprocessing import Pool
from readfq import readfq
from fileutil import opener
import numpy as np

    
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

def load_refgene_BED(filename, with_structure=False, only_pc=False, select_inverse=False):
    refgene = {}
    file = opener(filename)
    for line in file:
        if only_pc and ('protein_coding' not in line): continue
        line = line.strip().split('\t')
        chrom, start, end, name = line[0], int(line[1]), int(line[2]), line[3]
        if with_structure:
            lengths = [float(_) for _ in line[-2].split(',') if _!='']
            exon_starts = [float(_) for _ in line[-1].split(',') if _!='']
            if select_inverse:
                introns = []
                for i in range(len(exon_starts)-1):
                    introns.append([exon_starts[i]+lengths[i], exon_starts[i+1]])
                refgene[name] = [chrom, start, end, introns]
            else:
                exons = [[s, s+l] for s, l in zip(exon_starts, lengths)]
                refgene[name] = [chrom, start, end, exons]
        else:
            refgene[name] = [chrom, start, end]
    return refgene
            
    
def load_refasm_chrlength(filename):
    refasm_chrlen = {}
    file = opener(filename) 
    name, l = '', 0
    for qname, seq, qual in readfq(file):
        refasm_chrlen[qname] = len(seq)
    return refasm_chrlen


def worker(bundle):
    f, file_out_name, qnames, workerID = bundle
    file_out = open(file_out_name+'.tmp{0}'.format(workerID), 'w')
    if qnames!=None: 
        qnames = set(qnames)
    else:
        file_out.close()
        return
    file = opener(f)
    for line in file:
        line = line.strip().split('\t')
        if line[0] not in qnames: 
            continue

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
            if is_storeBothCov: 
                condition = np.sum(cov_relative>0.0)>=threshold_nb_exons and \
                                    np.max(cov_relative)>0.7 and \
                                    np.sum(covq_relative>0.5)<5  # intron restriction is loose
            else: 
                condition = np.sum(cov_relative>0.0)>1 and np.max(cov_relative)>0.7
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
                           ','.join([str(_) for _ in cov]),  # exon absolute coverage (bp)
                           ','.join(['%.2f'%_ for _ in cov_relative]),  # exon relative coverage
                           ','.join([str(_) for _ in covq]),  # intron absolute coverage (bp)
                           ','.join(['%.2f'%_ for _ in covq_relative])  # intron relative coverage
                           ]
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
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', action="store_true", default=False,
                        help='whether to output all locations')
    parser.add_argument('-t', action="store", type=int, default=1, help='number of threads')
    parser.add_argument('-e', action="store", type=int, default=3, 
                        help='require PPG candidate to cover at least INT exons')
    parser.add_argument('-o', type=str, default='PPGFparse', required=True,
                        help="prefix of intermediate and output files")
    parser.add_argument('gencode', action="store", type=str, help='gene reference in BED format (e.g. from Gencode).')
    parser.add_argument('ref', action="store", type=str, help='reference genome')
    parser.add_argument('aln', action="store", type=str, help='splice-aware alignment of long structural variants')

    if len(sys.argv)<2:
        parser.print_help()
        exit(1)


    args = parser.parse_args()
    is_outputall = args.a
    is_storeBothCov = True
    prefix_output = args.o
    nb_cpu = args.t
    threshold_nb_exons = args.e
    filename_gencode = args.gencode
    filename_ref = args.ref
    filename_input = args.aln

    print('[M] loading references...', file=sys.stderr)
    gencode = load_refgene_BED(filename_gencode, with_structure=True, only_pc=False)
    if is_storeBothCov: 
        # load intron structures
        gencode_inverse = load_refgene_BED(filename_gencode, with_structure=True, only_pc=False, select_inverse=True)
    hs38_chrlen = load_refasm_chrlength(filename_ref)

    if nb_cpu==1:
        worker([filename_input, filename_input, None, 0])
    else: 
        # probe
        total = 0
        qnames = []
        file = opener(filename_input)
        try:
            for line in file:
                qnames.append(line.split('\t')[0])
        except EOFError: 
            # dirty: allow annotation for a ongoing splice-aware alignment for a quick peek.
            #        Should not happen in normal runs.
            print('[W] EOF reading input file {0}; continue anyway.'.format(filename_input), file=sys.stderr)
        file.close()
        qnames = sorted(list(set(qnames)))
        print("[M] total {0} SVs aligned".format(len(qnames)), file=sys.stderr)

        # prepare
        step = ceil(len(qnames)/nb_cpu)
        print('[M] each thread will have {0} queries to process'.format(step), file=sys.stderr)
        qnames_bundled = [qnames[i*step:(i+1)*step] for i in range(nb_cpu)]
        bundles = [[filename_input, prefix_output, qnames_bundled[i], i] for i in range(nb_cpu)]

        # run
        with Pool(nb_cpu) as p:
            p.map(worker, bundles)
        print('[M] finished, concatenating tmp files.', file=sys.stderr)
    
    counter = 0
    with open(prefix_output+'.PPG.paf', 'w') as file_out:
        for i in range(nb_cpu):
            with open(prefix_output+'.tmp'+str(i)) as file:
                for line in file:
                    file_out.write(line)
                    counter+=1
            os.remove(prefix_output+'.tmp'+str(i))

    print('[M] done, {0} lines'.format(counter), file=sys.stderr)
    quit(0)
