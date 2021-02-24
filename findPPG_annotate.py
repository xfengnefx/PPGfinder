import sys, os
from readfq import readfq
from fileutil import opener
trans_revcomp = str.maketrans('ATCGatcg', 'TAGCtagc')
trans_nothing = str.maketrans('ATCGatcg', 'ATCGatcg')
import numpy as np
from multiprocessing import Pool

   
def find_TSD_r3e_core(FLKleft, FLKright):
    verbose = 0
    i1=0 
    i2 = 0
    i_miss=-1
    i_miss2 = -1
    TSD1 = ''
    TSD2 = ''
    nb_mismatch = 0
    shorter_length = min(len(FLKright), len(FLKleft))
    
    while True:
        if i1>=shorter_length or i2>=shorter_length:break
        base1 = FLKleft[i1]
        base2 = FLKright[i2]
        if verbose: print('aln at {0} {1} {2} {3}, nb mismatch {4}, {5}, {6}'.format(i1, i2, base1, base2, nb_mismatch, TSD1, TSD2))
        if base1!=base2:
            if nb_mismatch<=0:  # allow {} base off
                if i1+1>=shorter_length or i2+1>=shorter_length:
                    break
                nb_mismatch+=1
                if FLKleft[i1+1]==FLKright[i2+1]:  # substitute
                    i_miss = i1
                    i_miss2 = i2
                    i1+=1
                    i2+=1
                    TSD1+=base1
                    TSD2+=base2
                    continue
                elif FLKleft[i1]==FLKright[i2+1]:  # 1bp gap on left
                    i_miss = i1+1
                    i_miss2 = i2+1
                    i2+=1
                    TSD1+=base1
                    continue
                elif FLKleft[i1+1]==FLKright[i2]:  # 1bp gap on right
                    i_miss = i1+1
                    i_miss2 = i2+1
                    i1+=1
                    TSD2+=base2
                    continue
                else:  # next base is not a match, terminate
                    break
            else:
#                 print(TSD1, TSD2)
                break
        else: 
            TSD1+=base1
            TSD2+=base1
        i1+=1
        i2+=1
    # do not allow mismatch at the last base or the second last
    if i_miss!=-1:
        if i_miss+2>=i1 or i_miss+2>=i2:
            TSD1 = TSD1[:i_miss]
            TSD2 = TSD2[:i_miss2]
    if verbose:
        print('rough TSD1,', TSD1)
        print('      i_miss:', i_miss, 'i', i1)
        print('rough TSD2', TSD2, 'i', i2)
        print('      i_miss2:', i_miss2)
    return TSD1, TSD2, nb_mismatch
def find_TSD_r3e(FLKleft, FLKright):
    TSD1_, TSD2_, nb_mismatch_, code_ = '', '', 0, -1
    offset_1 = ''
    offset_2 = ''
    
    TSD1, TSD2, nb_mismatch = find_TSD_r3e_core(FLKleft, FLKright)
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2; offset_1=''; offset_2=''; nb_mismatch_=nb_mismatch; code_=0
    
    # examine cases with at most 2bp mis/indel at the begining
    TSD1, TSD2, nb_mismatch = find_TSD_r3e_core(FLKleft[1:], FLKright)
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2; offset_1=FLKleft[:1]; offset_2=''; nb_mismatch_=nb_mismatch; code_=1
    TSD1, TSD2, nb_mismatch = find_TSD_r3e_core(FLKleft[2:], FLKright)
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2; offset_1=FLKleft[:2]; offset_2=''; nb_mismatch_=nb_mismatch; code_=2
    
    TSD1, TSD2, nb_mismatch = find_TSD_r3e_core(FLKleft, FLKright[1:])
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2; offset_1=''; offset_2=FLKright[1]; nb_mismatch_=nb_mismatch; code_=3
    TSD1, TSD2, nb_mismatch = find_TSD_r3e_core(FLKleft, FLKright[2:])
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2;offset_1=''; offset_2=FLKright[:2]; nb_mismatch_=nb_mismatch; code_=4
    
    TSD1, TSD2,nb_mismatch = find_TSD_r3e_core(FLKleft[1:], FLKright[1:])
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2;offset_1=FLKleft[:1]; offset_2=FLKright[:1]; nb_mismatch_=nb_mismatch; code_=5
    TSD1, TSD2,nb_mismatch = find_TSD_r3e_core(FLKleft[1:], FLKright[2:])
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2;offset_1=FLKleft[:1]; offset_2=FLKright[:2]; nb_mismatch_=nb_mismatch; code_=5
    TSD1, TSD2,nb_mismatch = find_TSD_r3e_core(FLKleft[2:], FLKright[1:])
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2;offset_1=FLKleft[:2]; offset_2=FLKright[:1]; nb_mismatch_=nb_mismatch; code_=5
    TSD1, TSD2, nb_mismatch = find_TSD_r3e_core(FLKleft[2:], FLKright[2:])
    if min(len(TSD1), len(TSD2))>min(len(TSD1_), len(TSD2_)):
        TSD1_=TSD1; TSD2_=TSD2; offset_1=FLKleft[:2]; offset_2=FLKright[:2]; nb_mismatch_=nb_mismatch; code_=6
        
    return TSD1_, TSD2_, offset_1, offset_2, nb_mismatch_, code_


def match_length_from_cigar(cigar):
    ret = 0
    buf = ''
    for c in cigar:
        if c.isdigit():
            buf+=c
        else:
            if c=='M': ret+=int(buf)
            buf = ''
    return ret


def longest_conseq_polymer(seq, char, return_seq=False):
    cnt = 0
    max_cnt = 0
    buf = ''
    max_buf = ''
    for c in seq:
        if c!=char:
            if cnt>max_cnt:
                max_cnt = cnt
                max_buf = buf
            cnt = 0
            buf = ''
        else:
            cnt+=1
            buf+=c
    if cnt>max_cnt:
        max_cnt = cnt
        max_buf = buf
    if return_seq:
        return max_cnt, max_buf
    return max_cnt


def samplename_from_contigname(tigname):
    hap = ''
    tigname_short = tigname.split('#')[0]
    if 'h1' in tigname and 'h1' not in tigname_short:
        hap='H1'
    if 'H1' in tigname and 'H1' not in tigname_short:
        hap='H1'
    if 'h2' in tigname and 'h2' not in tigname_short:
        hap='H2'
    if 'H2' in tigname and 'H2' not in tigname_short:
        hap='H2'
    if hap!='':
        return tigname_short+'_'+hap
    return tigname_short


def samplename_from_SVname(SVname):
    return samplename_from_contigname(SVname.split('_')[1])


def longest_conseq_polymer_allowmismatch(seq, char, return_seq=False):
    cnt = 0
    max_cnt = 0
    buf = ''
    max_buf = ''
    is_mismatch = False
    for c in seq:
        if c!=char:
            if not is_mismatch:
                is_mismatch = True
                cnt+=1
                buf+=c
            else:  # 2 conseq miss, terminate
                is_mismatch = False
                cnt-=1
                buf = buf[:-1]
                if cnt>max_cnt:
                    max_cnt = cnt
                    max_buf = buf
                cnt = 0
                buf = ''
        else:
            is_mismatch = False
            cnt+=1
            buf+=c
    if cnt>max_cnt:
        max_cnt = cnt
        max_buf = buf
    if return_seq:
        return max_cnt, max_buf
    return max_cnt


def is_start_with_polyA_allowmismatch(seq0, chars='AT'):
    cnt = 0
    max_cnt = 0
    buf = ''
    max_buf = ''
    
    seqs = [seq0, seq0[::-1]]
    
    for char in chars:
        for seq in seqs:
#             print(char, seq)
            cnt = 0
            buf = ''
            is_mismatch = False
            for c in seq:
                if c!=char:
                    if not is_mismatch:
                        is_mismatch = True
                        cnt+=1
                        buf+=c
                    else:  # 2 conseq miss, terminate
                        is_mismatch = False
                        cnt-=1
                        buf = buf[:-1]
                        if cnt>max_cnt:
                            max_cnt = cnt
                            max_buf = buf
                        cnt = 0
                        buf = ''
                        break
                else:
                    is_mismatch = False
                    cnt+=1
                    buf+=c
            if cnt>max_cnt:
                max_cnt = cnt
                max_buf = buf
#             print(max_cnt, max_buf)
    if cnt>max_cnt:
        max_cnt = cnt
        max_buf = buf
    return max_cnt, max_buf













###########################
#         main            #
###########################
nb_cpu = 28
species = 'placeholder'
dir1 = '/path/to/3_annotation/'  # generated by step1
fs = [dir1+_ for _ in os.listdir(dir1) if _.endswith('PPGpaf')]


entries = []
special_retain = set()
for f in fs:
    print(f)
    load = 0
    with open(f) as file:
        for line in file:
            tmp = line.strip().split('\t')
            hit_exons = np.array([float(_) for _ in tmp[-3].split(',') if _!=''])
            hit_introns = np.array([float(_) for _ in tmp[-1].split(',') if _!=''])
            hit_introns_bases = np.array([float(_) for _ in tmp[-2].split(',') if _!=''])
            de_f = float([_ for _ in tmp if 'de:f' in _][0][5:])
            if (np.sum(hit_exons>0)==2 and np.sum(hit_introns>0.2)==0) or \
               (np.sum(hit_exons>0)>2 and np.sum(hit_introns>0.2)<=1 and np.sum(hit_introns_bases>100)<=1):  # gene structure
                if (int(tmp[1])-(int(tmp[3]) - int(tmp[2]))<500 or de_f<0.05):  # flanking
                    if (de_f<0.05):
                        special_retain.add(tmp[-6].split('|')[-1])
                    entries.append(tmp)
                    load+=1
        print(load)
print('nb entries:', len(entries))
contigs_names = set()
for entry in entries:
    contigs_names.add(entry[0].split('_')[1])
print('nb contigs:', len(contigs_names))
print('\nspecial:', special_retain)

# log exon/intron mapping status
exonintron = {}
for e in entries:
    key = e[0]+':'+e[-6].split('|')[0]
    exons = [float(_) for _ in e[-3].split(',') if _!='']
    introns = [float(_) for _ in e[-1].split(',') if _!='']
    exonintron[key] = [exons, introns]  # exon, intron

fs = ['/path/to/assemblies_contigs/'+f for f in os.listdir('/path/to/assemblies_contigs/')]
contigs = {}
for f in fs:
    with open(f) as file:
        for qname, seq, qual in readfq(file):
            contigs[qname] = seq
print(len(contigs))

#load hs38
hs38 = {}
file = opener('/path/to/hs38.fa')
for qname, seq, qual in readfq(file):
    hs38[qname] = seq

# load gene direction
gene_directions = {}
file = opener('/path/to/hs38.gencode.v31.bed.gz')
for line in file:
    line = line.split('\t')
    gene = line[3].split('|')[-1]
    direc = line[5]
    gene_directions[gene] = direc

# fix gene symbol
genename_map = {
    'AC004086.1':'RPL6',
    'DNAJC25-GNG10':'GNG10',
    'RPL36A-HNRNPH2':'RPL36A'
}
for gene in list(gene_directions.keys()):
    if gene in genename_map:
        gene_directions[genename_map[gene]] = gene_directions[gene]




# find PPG
yes_write = True
# sancheck_list = ['INS_HG00268#VRZF01001217.1_chr12-_125316601-125316601_2811756-2814895']
sancheck_list = []  # debug
SVname_request = ''  # debug

events = {}

def worker(d):
    entries, idx = d  # entries is data block, idx is thread ID
    if species=='human':file_dump = open('/path/to/output/PPGout_part{0}.tsv'.format(idx), 'w')
        
    for entry in entries:
        SVname = entry[0]
        SVname = SVname.replace('chrUn_', 'chrUn-')    
        TYPE, contig_name, chrom_s, coor_ref, coor_contig = SVname.split('_')
        e = '_'.join([TYPE, chrom_s[:-1], entry[-6].split('|')[-1]])
        if len(sancheck_list)>0 and entry[0] not in sancheck_list:continue
        parentgene = entry[-6].split('|')[-1]
        
        # correct a few gene symbols
        if parentgene in genename_map:
            tmp = entry[-6].split('|')
            entry[-6] = '|'.join([tmp[0], tmp[1], genename_map[tmp[-1]]])
            parentgene = genename_map[parentgene]
        parentgene_ID = entry[-6].split('|')[0]

        if TYPE=='INS':
            contig_start, contig_end = coor_contig.split('-')
            contig_start = int(contig_start)
            contig_end = int(contig_end)
            SVlength = contig_end-contig_start
            SVstrand = chrom_s[-1]
            aln_strand = entry[4]
            contig_seq = contigs[contig_name].upper()
            aln_start = int(entry[2])
            aln_end = int(entry[3])

            aln_strand_ori = aln_strand
            SV = contig_seq[contig_start:contig_end]

            if SVstrand=='-':
                contig_seq = contig_seq[::-1].translate(trans_revcomp)
                tmp = contig_start
                contig_start = len(contig_seq)-contig_end
                contig_end = len(contig_seq)-tmp
                tmp = aln_start
                aln_start = SVlength-aln_end
                aln_end = SVlength-tmp

                # correct aln strand
                if aln_strand=='-':aln_strand = '+'
                else: aln_strand = '-'

        else:
            contig_start, contig_end = coor_ref.split('-')
            contig_start = int(contig_start)
            contig_end = int(contig_end)
            aln_start = int(entry[2])
            aln_end = int(entry[3])
            SVstrand = chrom_s[-1]
            aln_strand = entry[4]

            contig_seq = hs38[chrom_s[:-1].replace('-', '_')].upper()
            SV = contig_seq[contig_start:contig_end]
        
        ############### TSD ################
        which = 0
        alignment = 'left'
        FLKleft = contig_seq[contig_start:contig_start+aln_start]#SV[:aln_start]
        FLKright = contig_seq[contig_end:contig_end+50]
        TSD1, TSD2, TSD_offset1, TSD_offset2, TSD_error, TSD_code = find_TSD_r3e(FLKleft, FLKright)

        FLKleft_ = contig_seq[contig_start-50:contig_start][::-1]
        FLKright_ = contig_seq[contig_start+aln_end:contig_end][::-1]#SV[aln_end:][::-1]
        TSD1_, TSD2_, TSD_offset1_, TSD_offset2_, TSD_error_, TSD_code_ = find_TSD_r3e(FLKleft_, FLKright_)

        if (len(sancheck_list)>0 or SVname_request!=''):
            print('left TSD code:', TSD_code, TSD1, TSD2)
            print('right TSD code:', TSD_code_, TSD1, TSD2)
        if min(len(TSD1_), len(TSD2_))>min(len(TSD1), len(TSD2)):
            TSD1 = TSD1_
            TSD2 = TSD2_
            TSD_offset1 = TSD_offset1_
            TSD_offset2 = TSD_offset2_
            TSD_error = TSD_error
            alignment = 'right'
            which = 1
        else:
            alignment = 'left'
        
        ############### polyA ############
        need_rev = False
        if alignment=='left':
            if (aln_strand=='+' and gene_directions[parentgene]=='-') or\
               (aln_strand=='-' and gene_directions[parentgene]=='+'):
                need_rev = True
                polyA = contig_seq[contig_start+len(TSD1):contig_start+aln_start]  # on the left

            else:
                polyA = contig_seq[contig_start+aln_end:contig_end]  # on the right
        else:
            if (aln_strand=='+' and gene_directions[parentgene]=='-') or\
               (aln_strand=='-' and gene_directions[parentgene]=='+'):
                need_rev = True
                polyA = contig_seq[contig_start:contig_start+aln_start]  # polya on the left
            else:
                polyA = contig_seq[contig_start+aln_end:contig_end-len(TSD2)+1]  # on the right
    
        cnt, tail = is_start_with_polyA_allowmismatch(polyA, 'A')
        cnt_, tail_ = is_start_with_polyA_allowmismatch(polyA, 'T')
        A_or_T = 'A'
        if cnt_>cnt:
            cnt = cnt_
            tail = tail_
            A_or_T = 'T'

        # concat and correct base
        the_retrocopy = ''
        if alignment=='left':
            the_retrocopy = contig_seq[contig_start-20:contig_start].lower() +\
                            contig_seq[contig_start:contig_end] + TSD2 +\
                            contig_seq[contig_end+len(TSD2):contig_end+len(TSD2)+20].lower()
        else:
            the_retrocopy = contig_seq[contig_start-len(TSD1)-20:contig_start-len(TSD1)].lower() +\
                            TSD1[::-1] + contig_seq[contig_start:contig_end] +\
                            contig_seq[contig_end:contig_end+20].lower()

        # convert to transcription strand
        if need_rev:
            the_retrocopy = the_retrocopy[::-1].translate(trans_revcomp)
            tmp = TSD1
            TSD1 = TSD2[::-1].translate(trans_revcomp)
            TSD2 = tmp[::-1].translate(trans_revcomp)
            polyA = polyA[::-1].translate(trans_revcomp)
        
            
        
        #summarize
        NMi = int([_ for _ in entry if 'NM:i' in _][0][5:])
        cigar = [_ for _ in entry if 'cg:Z' in _][0][5:]
        aln_length = match_length_from_cigar(cigar)
        exon_state = [float(_) for _ in entry[-3].split(',') if _!='']
        i=0
        while True:
            if exon_state[i]==0:
                i+=1
                continue
            if exon_state[i]>0.9:
                lost5=False; break
            else:
                lost5=True; break
        i=-1
        while True:
            if exon_state[i]==0:
                i-=1
                continue
            if exon_state[i]>0.9:
                lost3=False; break
            else:
                lost3=True; break
        
        # write summary
        if yes_write:
            if NMi/aln_length>0.15:  # human samples better do not allow high seq div
                if 'GGO' not in SVname and 'PTR' not in SVname and 'PAB' not in SVname:
#                         print("skipped {0} ({1})bcs seq div".format(entry[0], entry[-6]))
                    continue
            cleavage = the_retrocopy[16:20].lower() + the_retrocopy[20:24]
            cleavage = cleavage[::-1]
            holder = [SVname, samplename_from_SVname(SVname), 
                      TYPE, contig_name, chrom_s[:-1].replace('-', '_'), chrom_s[-1],
                      coor_ref, coor_contig, len(SV), 
                      NMi,
                      aln_start, aln_end,
                      NMi/aln_length,
                      parentgene_ID,  # gene ID
                      parentgene,  # gene name
                      cnt,  # polyA length
                      polyA,
                      min(len(TSD1), len(TSD2)),  # TSD length
                      TSD1,TSD2,
                      cleavage, 
                      the_retrocopy # 5' TSD ~ 3' TSD
                     ]
            file_dump.write('\t'.join([str(_) for _ in holder])+'\n')
    file_dump.close()
    
with Pool(28) as p:
    step = int(len(entries)/nb_cpu)+1
    packs = [[entries[i*step:(i+1)*step], i] for i in range(nb_cpu)]
    print('each thread gets', step)
    p.map(worker, packs)
    

# merge results from threads
os.system('cat /path/to/output/PPGout*part* > /pth/to/output/PPGout.tsv; rm /path/to/output/PPGout*part*')
    

# filter the results
data = {}
file = open('/path/to/output/PPGOut.tsv')
for line in file:
    line = line.strip().split('\t')
    
    #### tail
    tail_base = 'A'
    if longest_conseq_polymer(line[16], 'T')>longest_conseq_polymer(line[16], 'A'):
        tail_base = 'T'
    
    taillength = longest_conseq_polymer(line[16], tail_base)
    if len(line[16])==0 or taillength==0 or line[16].count(tail_base)/len(line[16])>0.7:
        line[15] = str(taillength)
        pureness = -1
        if len(line[16])!=0:
            pureness = line[16].count(tail_base)/len(line[16])
    else:
        line[15] = '-1'
        pureness = -1
        if len(line[16])!=0:
            pureness = line[16].count(tail_base)/len(line[16])

    tail = line[16]
    if len(tail)==0:
        polyArate = 0
    else:
        polyArate = int(line[15])/len(line[16])
    tail_cnt, tail_seq = is_start_with_polyA_allowmismatch(tail)
            
    ######### TSD and truncation
    TSDlength = len(line[18])
    lost5 = line[20]=='True'
    lost3 = line[21]=='True'
    

    ##### (init)
    flkl = int(line[8]) - (int(line[11])-int(line[10]))
    if line[0] not in data:
        data[line[0]] = {'polyArate':0, 'TSDlength':0, 'pureness':-2,
                            'flk':flkl,
                            'exon':sum(exonintron[line[0].replace('chrUn-', 'chrUn_')+':'+line[13]][0]),
                            'intron':sum(exonintron[line[0].replace('chrUn-', 'chrUn_')+':'+line[13]][1]),
                            '5trunc':True, '3trunc':True, 'entry':'', 'gene':line[14],
                            'tail':tail_seq}
    flkl_old = data[line[0]]['flk']
    if flkl_old<100 and flkl>=100:  # don't ever use a worse case
        continue

    ###### (try to update)
    score = 0
    if data[line[0]]['polyArate'] < polyArate: score+=2
    if data[line[0]]['TSDlength'] <TSDlength: score+=2
    
    if data[line[0]]['polyArate']<0.7 and polyArate<0.7:
        exons = sum(exonintron[line[0].replace('chrUn-', 'chrUn_')+':'+line[13]][0])
        introns = sum(exonintron[line[0].replace('chrUn-', 'chrUn_')+':'+line[13]][1])
#         if introns>data[line[0]]['intron']:continue
        if exons+1<data[line[0]]['exon']:score-=2
        if exons-data[line[0]]['exon']>=0.5 and flkl<=flkl_old:
            score+=2
    
    if score>1:#or (flkl<100 and flkl_old>100):
        data[line[0]]['polyArate'] = polyArate
        data[line[0]]['TSDlength'] = TSDlength
        data[line[0]]['pureness'] = pureness
        data[line[0]]['flk'] = flkl
        data[line[0]]['exon'] =exons
        data[line[0]]['intron'] =introns
        data[line[0]]['entry'] = line
file.close()

file_out = open('/path/to/output/PPGout_filtered.tsv', 'w')

for key in data:
    if data[key]['entry']=='':
        continue
        
    # tail length correction:
    tail = data[key]['entry'][16]
    if len(tail)>0:
        cnt, seq = is_start_with_polyA_allowmismatch(tail)
        data[key]['entry'][15] = str(cnt)

    # a manual filter used in manuscript
    if species=='human' and data[key]['entry'][14]=='AK2':continue
    if species=='human' and data[key]['entry'][14]=='RPL6':continue
    if species=='human' and data[key]['entry'][14]=='RPSA':continue
    if species=='human' and data[key]['entry'][14]=='FMN1':continue
    if species=='human' and data[key]['entry'][14]=='TERT':continue

    #final flk check
    if data[key]['flk']>100:continue

    # final polyA/TSD sancheck
    if min(len(data[key]['entry'][18]), len(data[key]['entry'][19]))<6:  # short TSD
        continue
    if int(data[key]['entry'][15])<6:  # polyA
        cleavage = data[key]['entry'][-2]
        cleavage_left = cleavage[:4]
        if len(data[key]['entry'][18])<10:   # require TSD
#                 print('drop because tail(a):', data[key]['entry'][14])
            continue
        if cleavage.count('A')<3 and cleavage.count('T')<3:   # require cleavage if no tail
#             if cleavage.count('C')+cleavage.count('G')>=3:
#                 print('drop because tail(a):', data[key]['entry'][14])
            continue  

    tail = data[key]['entry'][16]
    if len(tail)>0:
        tail_cnt, tail_seq = is_start_with_polyA_allowmismatch(tail)
    data[key]['entry'][16] = tail_seq

    file_out.write('\t'.join(data[key]['entry'])+'\n')
file_out.close()
