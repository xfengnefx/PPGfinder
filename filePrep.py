import sys, gzip, argparse
import binascii

def is_gz(f):
    with open(f, 'rb') as file:
        return binascii.hexlify(file.read(2)) == b'1f8b'

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

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', type=bool, default=False)  # is_vcf
    parser.add_argument('filename', type=str)
    args = parser.parse_args()
    is_vcf = args.v
    f = args.filename
    gzipped = is_gz(f)
    if gzipped: file = gzip.open(f, 'rb')
    else: file = open(f)
    if is_vcf: vcf2fa(file, gzipped)
    else: var2fa(file, gzipped)
    file.close()
