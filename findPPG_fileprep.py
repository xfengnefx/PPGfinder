import os, sys, gzip, binascii
import argparse

TOTAL = 0

def is_gz(filename):
    "Tells if the given filename leads to invalid file (-1), gzipped file (0), not non-gzipped file(1)."
    try:
        file = open(filename, 'rb')
    except:
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


def var2fa(stream, input_name, stream_out):
    global TOTAL
    for line in stream:
        if line[0]!='V': continue
        line = line.strip().split('\t')
        _1, chrom, start, end, _2, _3, ref, alt, queryname, q_start, q_end, strand = line
        if abs(len(ref)-len(alt))<50: continue  # not long enough
        if len(ref)>len(alt):
            # newname = 'DEL_'+'_'.join([queryname, chrom+strand, start+'-'+end, q_start+'-'+q_end])
            newname = ['DEL', queryname, chrom, strand, start, end, q_start, q_end]
            seq = ref.upper()
        else:
            # newname = 'INS_'+'_'.join([queryname, chrom+strand, start+'-'+end, q_start+'-'+q_end])
            newname = ['INS', queryname, chrom, strand, start, end, q_start, q_end]
            seq = alt.upper()
        sys.stdout.write('>SV{0}\n{1}\n'.format(TOTAL, seq))
        stream_out.write('\t'.join(['SV'+str(TOTAL)]+[str(_) for _ in newname]+[input_name])+'\n')
        TOTAL+=1
    return 0


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Write long SV seqs in fasta foramt to stdout.")
    parser.add_argument('-o', type=str, required=True, default='',
                        help='Output file containing location info of the SVs. Tab-delimited.')
    parser.add_argument('filename_variant_calling', type=str, nargs='+')
    
    if (len(sys.argv)==1):
        parser.print_help()
        exit(1)
    args = parser.parse_args()


    assert(args.o.strip()!='')
    fs = args.filename_variant_calling

    with open(args.o, 'w') as file_out:
        for fn in fs:
            sys.stderr.write('[M] at {0}...'.format(fn))
            file = opener(fn)
            var2fa(file, fn, file_out)
            sys.stderr.write(' total SVs now: {0}.\n'.format(TOTAL))