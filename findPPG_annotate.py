import os, sys, gzip, sys, re, pickle, argparse
from readfq import readfq
from fileutil import opener
from xsimple import *
import numpy as np
from multiprocessing import Pool
from collections import Counter

base_for = "ACGT"
base_rev = "TGCA"
convert = str.maketrans(base_for, base_rev)
import editdistance

    
def loadGenesNotTrusted(f_plain, f_conditioned):
    genes_not_trusted = []
    genes_not_trusted_conditioned = {}
    with open(f_plain) as file:
        for line in file:
            if line[0]=='!': continue
            genes_not_trusted.append(line.strip())
    print('will not trust {0} genes'.format(len(genes_not_trusted)))
    if f_conditioned==None: return f_plain
    with open(f_conditioned) as file:
        for line in file:
            if line[0]=='!': continue
            line = line.strip().split('\t')
            if line[1] not in genes_not_trusted_conditioned: genes_not_trusted_conditioned[line[1]] = []
            genes_not_trusted_conditioned[line[1]].append(line[0])
    return genes_not_trusted, genes_not_trusted_conditioned

def filename2shortname(filename, mode):
    """Convert filename to a shorter name. User please configure this if needed."""
    return filename
    
def loadPPG(fs, genes_not_trusted, genes_not_trusted_conditioned, trusted_SVnames, mode='.',
            threshold_least_exons=3, threshold_most_introns=2,
            threshold_exon_cov=0.2, threshold_intron_cov=0.2,
            is_verbose=False, is_debug=False):
    PPGs_SVpov = {}
    PPGs_SVpov_features = ['aln_l', 'aln_qs', 'aln_qe',
                           'aln_strand', 'aln_chrom', 'aln_rs', 'aln_re' ,
                           'NM_i', 'de_f', 'cigar', 'cs',
                           'targetname_long', 'targetname_id', 'targetname_symbol',
                           'PPG_s', 'PPG_e',
                           'cov_bases', 'cov_relative',
                           'cov_bases_introns', 'cov_relative_introns']
    for f in fs:
        shortname = filename2shortname(f.split('/')[-1], mode)
        with open(f) as file:
            print('parsing '+f)
            for line in file:
                line = line.strip().split('\t')
                if trusted_SVnames!=None:
                    if line[0] not in trusted_SVnames: 
                        if is_verbose: print('skipped {0} because SV not in trusted SVnames'.format(line[0]))
                        continue  # CIC case
                parentgene = line[-6].split('|')[-1].split(',')[0]
                if parentgene in genes_not_trusted:
                    if is_verbose: print('skipped {0} because parent gene in genes_not_trusted'.format(line[0]))
                    continue
                if genes_not_trusted_conditioned!=None:
                    if (shortname in genes_not_trusted_conditioned) and (parentgene in genes_not_trusted_conditioned[shortname]):
                        if is_verbose: print('skipped {0} because parent gene in genes_not_trusted_conditioned'.format(line[0]))
                        continue
                longname = shortname+'@'+line[0]
                if 'chrUn' in longname:
                    longname = longname.replace('chrUn_', 'chrUn-')
                if longname not in PPGs_SVpov: 
                    PPGs_SVpov[longname] = {feature:[] for feature in PPGs_SVpov_features}
                # parse info from SV's point of view
                if 'protein_coding' not in line[-6]:
                    print('alignment probably too short at:', longname, '\t', len(line))
                    continue
                PPGs_SVpov[longname]['aln_l'].append(int(line[1]))  # length of splic alignment
                PPGs_SVpov[longname]['aln_qs'].append(int(line[2]))  # query start of splice alignment
                PPGs_SVpov[longname]['aln_qe'].append(int(line[3]))  # query end of splice alignment
                PPGs_SVpov[longname]['aln_strand'].append(line[4])
                PPGs_SVpov[longname]['aln_chrom'].append(line[5])
                PPGs_SVpov[longname]['aln_rs'].append(int(line[7]))  # ref start of the splice alignment
                PPGs_SVpov[longname]['aln_re'].append(int(line[8]))  # ref end of the splice alignment
                PPGs_SVpov[longname]['NM_i'].append(int(line[12][5:]))
                PPGs_SVpov[longname]['de_f'].append(float([_ for _ in line if _.startswith('de:f')][0][5:]))
                PPGs_SVpov[longname]['cigar'].append([_ for _ in line if _.startswith('cg:Z')][0][5:])
                PPGs_SVpov[longname]['cs'].append([_ for _ in line if _.startswith('cs:Z')][0][5:])
                PPGs_SVpov[longname]['targetname_long'].append(line[-6])
                PPGs_SVpov[longname]['targetname_id'].append(line[-6].split('|')[0])
                PPGs_SVpov[longname]['targetname_symbol'].append(line[-6].split('|')[-1].split(',')[0])  #.split('-')[0])  # , case for merged gene reference
                PPGs_SVpov[longname]['PPG_s'].append(int(line[-5].split(',')[0]))
                PPGs_SVpov[longname]['PPG_e'].append(int(line[-5].split(',')[1]))
                PPGs_SVpov[longname]['cov_bases'].append(np.array([float(_) for _ in line[-4].split(',')]))
                PPGs_SVpov[longname]['cov_relative'].append(np.array([float(_) for _ in line[-3].split(',')]))
                PPGs_SVpov[longname]['cov_bases_introns'].append(np.array([float(_) for _ in line[-2].split(',')]))
                PPGs_SVpov[longname]['cov_relative_introns'].append(np.array([float(_) for _ in line[-1].split(',')]))
    print('documented {0} fragments from {1} files'.format(len(PPGs_SVpov), len(fs)))           
    for SVname in list(PPGs_SVpov.keys()):
        if len(PPGs_SVpov[SVname]['aln_l'])==0: 
            del PPGs_SVpov[SVname]
            continue
    print('removed empty entries, left with {0} fragments'.format(len(PPGs_SVpov)))

    # half-legacy: define the best match for each fragment
    sancheck = []
    for longname in PPGs_SVpov:
        tmp_i, tmp_score = 0, 0
        content = PPGs_SVpov[longname]
        for i in range(len(content['aln_l'])):
            cigar = re.findall('[0-9]+[MNID]', content['cigar'][i])
            blocks = []
            for c in cigar:
                if int(c[:-1])<100: 
                    if len(blocks)>0: blocks[-1]+=int(c[:-1])
                    else: blocks.append(int(c[:-1]))
                else: blocks.append(int(c[:-1]))
            score = np.sum(blocks)

            if score>tmp_score:
                tmp_score = score
                tmp_i = i
        PPGs_SVpov[longname]['_best'] = [tmp_i, tmp_score]

    # confident SVs
    san_genes_del, san_genes_ins = [], []
    confident_SVnames = []
    for SVname in PPGs_SVpov:
        content = PPGs_SVpov[SVname]
        flag_at_least_one_alignment_passed_exon, flag_at_least_one_alignment_passed_intron = False, False
        flag_at_least_one_alignment_passed = False
        for idx_best in range(len(content['aln_qs'])):  # check everything, not just the presumed best one
            gene = content['targetname_symbol'][idx_best]
            if gene in genes_not_trusted: 
                continue
            if 'DEL' in SVname: san_genes_del.append(gene)
            if 'INS' in SVname: san_genes_ins.append(gene)
            cov = content['cov_bases'][idx_best]
            cov_relative = content['cov_relative'][idx_best]
            cov_introns = content['cov_bases_introns'][idx_best]
            cov_relative_introns = content['cov_relative_introns'][idx_best]
            if threshold_exon_cov<=1:  # is a percentage
                if np.sum(cov_relative>=threshold_exon_cov)>=threshold_least_exons: 
                    flag_at_least_one_alignment_passed_exon = True
                    cnt_exon = np.sum(cov_relative>=threshold_exon_cov)
            else:  # an absolute value is given; will assume a very forgiving percentage requirement
                if np.sum(cov>=threshold_exon_cov)>=threshold_least_exons and np.sum(cov_relative>=0.2)>=threshold_least_exons:
                    flag_at_least_one_alignment_passed_exon = True
                    cnt_exon = np.sum(cov>=threshold_exon_cov)
            if threshold_intron_cov<=1:  # is a percentage
                nb_intron_covered = np.sum(cov_relative_introns>=threshold_intron_cov)
                if nb_intron_covered<=threshold_most_introns:
                    flag_at_least_one_alignment_passed_intron = True
                    cnt_intron = nb_intron_covered
            else:  # an absolute value is given; will assume a very forgiving percentage requirement
                nb_intron_covered = np.sum(cov_introns>=threshold_intron_cov)
                if nb_intron_covered<=threshold_most_introns and (np.sum(cov_relative_introns>=0.1)<=threshold_most_introns):
                    flag_at_least_one_alignment_passed_intron = True
                    cnt_intron = nb_intron_covered
            if flag_at_least_one_alignment_passed_exon and flag_at_least_one_alignment_passed_intron:
                if cnt_exon==2 and cnt_intron>0:continue  # do not allow 2-exon retrocopies to retain introns
                flag_at_least_one_alignment_passed = True
                break
        if flag_at_least_one_alignment_passed:
            confident_SVnames.append(SVname)
    confident_SVnames = sorted(confident_SVnames)
    nb_del = len([_ for _ in confident_SVnames if 'DEL' in _])
    nb_ins = len([_ for _ in confident_SVnames if 'INS' in _])
    san_genes_del = sorted(list(set(san_genes_del)))
    san_genes_ins = sorted(list(set(san_genes_ins)))
    #print('selected {0} SVs ({1} dels, {2} inss)'.format(len(confident_SVnames), nb_del, nb_ins))
    #print('deleterious:', len(san_genes_del))
    #print(', '.join(san_genes_del))
    #print('insertive:', len(san_genes_ins))
    #print(', '.join(san_genes_ins))
    #print('\n-------------------------------------------------------------------\n')
    return PPGs_SVpov, confident_SVnames

def simpleSW(s1, s2, threshold=None, t_match=3, t_mismatch=-2, t_indel=-2):
    if threshold==None:
        threshold=(min(len(s1), len(s2)))*2
    mat = np.zeros([len(s1)+1, len(s2)+1])
    best_i, best_j, best_score = 0, 0, -100
    for i_ in range(len(s1)):  # rows
        for j_ in range(len(s2)):  # columns
            i = i_+1
            j = j_+1
            upper = mat[i-1, j]+t_indel
            left = mat[i, j-1]+t_indel
            if s1[i_]==s2[j_]:
                leftupper = mat[i-1, j-1] + t_match
            else:
                leftupper = mat[i-1, j-1] +t_mismatch
            score = max([upper, left, leftupper])
            mat[i, j] = score
            if score>=best_score:
                best_score = score
                best_i = i
                best_j = j
    if best_score<threshold: 
        return None, None
    alignment = [[best_i-1, best_j-1]]
    i = best_i
    j = best_j
    while True:
        tmp = [mat[i-1, j], mat[i, j-1]]
        if max(tmp)<=mat[i-1, j-1]: which = 2
        else:
            which = np.argmax(tmp)
        if which==0: 
            alignment.append([i-1-1, j-1])
            i -=1
        elif which==1:
            alignment.append([i-1, j-1-1])
            j -=1
        elif which==2:
            alignment.append([i-1-1, j-1-1])
            i -=1
            j -=1
        else:print('??')
        if i<=1 or j<=1:break
    return [[_[0] for _ in alignment][::-1], [_[1] for _ in alignment][::-1]], best_score

def count_polyA(seq):
    seq = seq.upper()
    count = max(seq.count('T'), seq.count('A'))
    return count

def match_polyA(seq):
    indices = []
    for k in range(8, 100):
        flag = False
        for i in range(0, len(seq)-k+1):
            s  = seq[i:i+k]
            rate = count_polyA(s)/len(s)
            if rate>0.9:
                indices.append([i, k, rate, s])
                flag = True
        if not flag: break
    return indices
            

def find_polyA(left, right, idx_left, idx_right, k_TSD, threshold_polyA=0.7, regardlessofTSD=False):
    if regardlessofTSD:
        indices_left = match_polyA(left)
        l, r = None, None
        if len(indices_left)>0:
            indices_left.sort(key=lambda x:(x[0], x[1], x[2]))
            l = indices_left[0]
        indices_right = match_polyA(right)
        if len(indices_right)>0:
            indices_right.sort(key=lambda x:(len(right)-x[0], x[1], x[2]))
            r = indices_right[0]
        if len(indices_left)==0 and len(indices_right)==0:
            return None
        return [l, r] 
    remnants = [left[idx_left:], right[:idx_right]]
    if len(remnants[0])==0:
        Arate_left = 0
    else:
        Arate_left = count_polyA(remnants[0])/len(remnants[0])
    if len(remnants[1])==0:
        Arate_right = 0
    else:
        Arate_right = count_polyA(remnants[1])/len(remnants[1])
    if max(Arate_left, Arate_right)<threshold_polyA: 
        return None
    if Arate_left>Arate_right:
        which = 1  # right side is 5'
        tail = remnants[0]
        if len(tail)<5: return None  # tail too short
        cleavage = right[idx_right-10:idx_right].lower()+right[idx_right:idx_right+k_TSD].upper() + right[idx_right+k_TSD:idx_right+k_TSD+10].lower()  # because part of the TTTT/AA might be inside the TSD
    else: 
        which = 0  # left side is 5'
        tail = remnants[1]
        if len(tail)<5: return None  # tail too short
        cleavage = left[idx_left-10:idx_left].lower() + left[idx_left:idx_left+k_TSD].upper() + left[idx_left+k_TSD:idx_left+k_TSD+10].lower()
    return tail, cleavage, which

def find_L1hallmarks(left, right, forceSW=False, threshold_polyA=0.7, is_debug=False):
    if len(left)!=len(right):
        raise ValueError
    which = None
    if is_debug:
        forceSW = True

    ret = []
    max_k = 5
    # find perfect matches
    for k in range(5, 20):
        bag1 = set([left[i:i+k] for i in range(len(left)-k+1)])
        bag2 = set([right[i:i+k] for i in range(len(right)-k+1)])
        ovlp = bag1.intersection(bag2)
        if len(ovlp)==0:
            break
        for mer in ovlp:
            idx_left = max([_.span()[1] for _ in re.finditer(mer, left)])
            idx_right = min([_.span()[0] for _ in re.finditer(mer, right)])
            a = find_polyA(left, right, idx_left, idx_right, len(mer))
            if a==None: continue
            tail, cleavage, which = a
            ret.append([len(mer), 0, [idx_left, idx_right], [mer, mer], tail, cleavage, which])  # 0 for "non-SW match"
            max_k = k
    if not forceSW:
        if k>=8 and not forceSW:
            if len(ret)==0:
                return ret, -1
            else:
                ret.sort(reverse=True)
                return ret[:3], 0
    if is_debug:
        print('ret of kmer matching:', ret)
    ret_perfect = ret
    ret = []
    # resort to SW
    if forceSW and max_k>=8:  # if TSD is already significant, only accept if SW gives much better result (i.e. merge two close perfect hits)
        k_start = max_k+5
    else:  # if nothing good was found, be more permissive
        k_start = max_k
    for k in range(k_start, 20):  # require 
        if not is_debug:
            if k<10: nb_mismatch = 1 
            else: nb_mismatch = 2
        else:
            nb_mismatch=5
        bag1 = list(set([left[i:i+k] for i in range(len(left)-k+1)]))
        bag2 = list(set([right[i:i+k] for i in range(len(right)-k+1)]))
        flag_termi = True
        for i in range(len(bag1)):
            best_in_this_iter = None
            max_score = 0
            for j in range(len(bag2)):
                if bag1[i]==bag2[j]: continue
                alignment, SWscore = simpleSW(bag1[i], bag2[j], threshold=-100)
                if SWscore==None: continue
                if SWscore<(3*k-5*nb_mismatch): continue  
                if SWscore<max_score: continue
                idx_left = max([_.span()[1] for _ in re.finditer(bag1[i], left)])
                idx_right = min([_.span()[0] for _ in re.finditer(bag2[j], right)])
                a = find_polyA(left, right, idx_left, idx_right, k)
                if a==None: continue
                tail, cleavage, which = a
                max_score = SWscore
                best_in_this_iter = [j, idx_left, idx_right, tail, cleavage, which]
            if best_in_this_iter==None: continue  # no match
            j, idx_left, idx_right, remnants, cleavage, which = best_in_this_iter
            ret.append([k, 1, [idx_left, idx_right], [bag1[i], bag2[j]], tail, cleavage, which])  # 1 for "SW match"
            max_k = k
            flag_termi = False
        if flag_termi: break
    if is_debug:
        print('ret of SW:', ret)
    if len(ret)==0 and len(ret_perfect)==0:
        return [], -1
    else:
        ret_perfect.sort(reverse=True)
        ret.sort(reverse=True)
        return ret[:3]+ret_perfect[:3], 0
        
def worker(bundle):
    arg, code = bundle
    if code!=2:
        is_unit_test = False
        worker_start, worker_end, confident_SVs = arg
        if code: debug = True
        else: debug=False
    else:
        is_unit_test = True
        SVnames = arg
        debug = False

    alldata = {}
    no = 0
    for worker_idx, SVname in enumerate(sorted(confident_SVs)):
        if not debug:
            if worker_idx<worker_start: continue
            if worker_idx>=worker_end: 
                break
        if is_unit_test:
            if SVname not in SVnames:continue
        alldata[SVname] = {}
        # prepare the sequence
        content = PPG_100flk[SVname]
        SVseq = SVseqs[SVname.split('@')[1]]
        left = SVflks[SVname.split('@')[1]]['l']
        right = SVflks[SVname.split('@')[1]]['r']
        SVseq_ = left+SVseq+right

        ##################################
        # V2
        # - search for TSDs first, score them and choose the longest yet closest possible pair
        # - expect polyA to exist between a TSD and the retrocopy on one end
        ##################################
        alnindices = [i for i in range(len(content['aln_qs']))]
        rets = []
        for i_aln in alnindices:  # iterate over all possible alignments
            strand = content['aln_strand'][i_aln]
            parentgene = content['targetname_symbol'][i_aln]
            parentgeneID = content['targetname_id'][i_aln]
            parentgene_ = parentgene + '|' + SVname.split('@')[1].split('_')[0] + '|' + SVname.split('@')[1].split('_')[2][:-1] +'|'+SVname.split('@')[0]  # gene|type|chrom|sample
            qs = content['aln_qs'][i_aln] + 100  # adjust for the flakning
            qe = content['aln_qe'][i_aln] + 100 # adjust for the flakning
            ql = content['aln_l'][i_aln]
            strand_SV = content['aln_strand'][i_aln]
            strand_gene = transcriptDirection[parentgeneID]
            SVseq = SVseq_[qs:qe]  # retrocopy
            SVflk_left = SVseq_[qs-100:qs]  # retrocopy flk left
            SVflk_right = SVseq_[qe:qe+100]  # retrocopy flk right

            # find TSD
            ret, status = find_L1hallmarks(SVflk_left, SVflk_right, forceSW=True, is_debug=debug)
            if status == -1:
                continue
            else:
                rets.extend([_+[i_aln, ] for _ in ret])
        if len(rets)==0:
            alldata[SVname]['has_TSD'] = False
            no+=1
            idx_best = PPG_100flk[SVname]['_best'][0]
            for key in PPG_100flk[SVname]:
                if key=='_best': continue
                alldata[SVname][key] = PPG_100flk[SVname][key][idx_best]
            tmp = find_polyA(SVflk_left, SVflk_right, None, None, None, regardlessofTSD=True)
            if tmp==None:
                alldata[SVname]['has_polyA'] = False
            else:
                l, r = tmp
                parentgene = PPG_100flk[SVname]['targetname_symbol'][idx_best]
                if transcriptDirection[parentgeneID]==PPG_100flk[SVname]['aln_strand'][idx_best]:  # which side *should* have a tail
                    if r==None:
                        alldata[SVname]['has_polyA'] = False
                        continue
                    idx_tail_start, tail_length, A_rate, tail = r
                    alt = l
                    which = 0
                else:
                    if l==None:
                        alldata[SVname]['has_polyA'] = False
                        continue
                    idx_tail_start, tail_length, A_rate, tail = l
                    alt = r
                    which = 1
                alldata[SVname]['has_polyA'] = True
                alldata[SVname]['idxTailstart'] = idx_tail_start
                alldata[SVname]['tail'] = tail
                alldata[SVname]['tail_alt'] = alt
                alldata[SVname]['which'] = which
        else:
#             rets.sort(key=lambda x: (x[2][0], 100-x[2][1], x[0]), reverse=True)  # sort annotations by TSD location
            rets.sort(key=lambda x:(len(x[3][0]), len(x[4])), reverse=True)  # sort by TSD length and polyA length
            if is_unit_test:  # print all hits
                for tmp in rets:
                    k, alntype, [idx_left, idx_right], [mer1, mer2], tail, cleavage, which, idx_best = tmp
                    if which==0:
                        tmp[-3] = tmp[-3][::-1]
                    print('is_unit_test debug:', tmp)
            k, alntype, [idx_left, idx_right], [mer1, mer2], tail, cleavage, which, idx_best = rets[0]
            for key in PPG_100flk[SVname]:
                if key=='_best': continue
                alldata[SVname][key] = PPG_100flk[SVname][key][idx_best]            
            alldata[SVname]['has_polyA'] = True
            alldata[SVname]['has_TSD'] = True
            alldata[SVname]['idxTSDstart_left'] = idx_left
            alldata[SVname]['idxTSDstart_right'] = idx_right
            alldata[SVname]['TSD'] = [mer1, mer2]
            alldata[SVname]['tail'] = tail
            alldata[SVname]['TSDalntype'] = alntype
            alldata[SVname]['TSDlength'] = k
            alldata[SVname]['which'] = which
            alldata[SVname]['cleavage'] = cleavage
            # log subseqs
            tmpseq = SVflks[SVname.split('@')[1]]['l'].lower() + SVseqs[SVname.split('@')[1]].upper() + SVflks[SVname.split('@')[1]]['r'].lower()
            idx_TSDstart_left = alldata[SVname]['aln_qs']+100-100+idx_left
            idx_TSDstart_right = alldata[SVname]['aln_qe']+100+idx_right
            segs = [tmpseq[:idx_TSDstart_left-k],  # left non-retro stuff
                    tmpseq[idx_TSDstart_left-k:idx_TSDstart_left],  # left TSD
                    tmpseq[idx_TSDstart_left:alldata[SVname]['aln_qs']+100],  # anything between left TSD and start
                    tmpseq[alldata[SVname]['aln_qs']+100:alldata[SVname]['aln_qe']+100],  # retrocopy
                    tmpseq[alldata[SVname]['aln_qe']+100:idx_TSDstart_right],  # anything between end and right TSD
                    tmpseq[idx_TSDstart_right:idx_TSDstart_right+k], # right TSD
                    tmpseq[idx_TSDstart_right+k:]  # right non-retro stuff
                   ]
            alldata[SVname]['segments'] = segs
    # print(no)
    return alldata

            
if __name__=='__main__':
    if len(sys.argv)==1:
        print('annotatePPG.py : annotates PPGpaf generated by findPPG.py')
        exit(0)
    print(os.getcwd())
    parser = argparse.ArgumentParser()
    # basics
    parser.add_argument('-t', action='store', help='number of threads', type=int, default=4)
    parser.add_argument('-f', action='store', help='directory of SV fastas', type=str, default=False)
    parser.add_argument('-l', action='store', help='directory of SV flanking fastas', type=str, default=False)
    parser.add_argument('-p', action='store', help='directory of PPGpafs', type=str, default=False)
    parser.add_argument('-s', action='store', help='a single PPGpaf file to be parsed', default=False)
    parser.add_argument('-o', help='output name', type=str, default=False)
    parser.add_argument('-v', action='store_true', help='is verbose', default=False)
    parser.add_argument('-d', action='store_true', help='is debug', default=False)
    # threhsolds
    parser.add_argument('--trustedSVs', dest='trustedSVs', help='a list of SVnames that are not CIC', type=str, action='store', default=None)
    parser.add_argument('--nb_exons', dest='nb_exons', help='needs to have at least $ exons', type=int, action='store', default=3)
    parser.add_argument('--cov_exons', dest='cov_exons', help='exon coverage threshold', type=float, action='store', default=0.2)
    parser.add_argument('--nb_introns_proportional', dest='nb_introns_proportional', type=float, action='store', default=False)
    parser.add_argument('--nb_introns', dest='nb_introns', help='needs to have at most $ introns', type=int, action='store', default=1)
    parser.add_argument('--cov_introns', dest='cov_introns', help='intron coverage threshold', type=float, action='store', default=0.2)
    
    args, unkown_args = parser.parse_known_args(sys.argv[1:])
    dir_fasta = args.f
    dir_flk = args.l
    dir_PPGpafs = args.p
    dir_one_PPGpaf = args.s
    nb_cpu = args.t
    outputname = args.o
    nb_exons = args.nb_exons
    nb_introns = args.nb_introns
    cov_exons = args.cov_exons
    cov_introns = args.cov_introns
    trustedSVs = args.trustedSVs
    is_verbose = args.v
    is_debug = args.d
    assert (dir_fasta and dir_flk and (dir_PPGpafs or dir_one_PPGpaf) and nb_cpu)
    if trustedSVs==None and is_verbose:
        print("[warning] not specifying list of trusted SVnames. Are you sure there's no CIC case here?")
    if not dir_fasta.endswith('/'): dir_fasta+='/'
    if not dir_flk.endswith('/'): dir_flk+='/'
    if not dir_PPGpafs.endswith('/'): dir_PPGpafs+='/'
    from math import floor

    # load confident SVnames (aka remove CIC)
    if trustedSVs!=None:
        trusted_SVnames = []
        with open(trustedSVs) as file:
            for line in file:
                trusted_SVnames.append(line.strip())
    else:
        trusted_SVnames = None

    # load transcript direction
    transcriptDirection = {}
    file = opener('hs38.gencode.v31.bed.gz')
    for line in file:
        line = line.split('\t')
        geneID = line[3].split('|')[0]
        assert (geneID not in transcriptDirection)
        transcriptDirection[geneID] = line[5]

    genes_not_trusted, genes_not_trusted_conditioned = loadGenesNotTrusted('genes_not_trusted.txt', 'genes_not_trusted.with_condition.txt')  # incoporate manual inspection

    if dir_PPGpafs:
        fs = [dir_PPGpafs+_ for _ in os.listdir(dir_PPGpafs) if _.endswith('PPGpaf')]
    else:
        fs = [dir_one_PPGpaf]
    if len(fs)==0:
        print('no file is selected')
        if dir_PPGpafs:
            print('looking at:', dir_PPGpafs)
        else:
            print('looking at file:', dir_one_PPGpaf)
    print('selecting {0} files'.format(len(fs)))
    PPG_100flk, confident_SVs = loadPPG(fs, genes_not_trusted, genes_not_trusted_conditioned, trusted_SVnames, 
                                        threshold_least_exons=nb_exons, threshold_most_introns=nb_introns,
                                        threshold_exon_cov=cov_exons, threshold_intron_cov=cov_introns,  # CIC cases are filtered here
                                        is_verbose=is_verbose, is_debug=is_debug)

    SVseqs = {}
    fs = [dir_fasta+_ for _ in os.listdir(dir_fasta) if _.endswith('.fa')]
    qnames = [_.split('@')[1] for _ in PPG_100flk.keys()]
    for f in fs:
        with open(f) as file:
            for qname, seq, qual in readfq(file):
                if qname in qnames:
                    SVseqs[qname] = seq
    print('loaded {0} SV seqs from {1} files'.format(len(SVseqs), len(fs)))

    SVflks = {}
    fs = [dir_flk+_ for _ in os.listdir(dir_flk) if _.endswith('.fa')]
    qnames = [_.split('@')[1] for _ in PPG_100flk.keys()]
    for f in fs:
        with open(f) as file:
            for qname, seq, qual in readfq(file):
                if qname[:-2] in qnames:
                    if qname[:-2] not in SVflks:
                        SVflks[qname[:-2]] = {'r':None, 'l':None}
                    SVflks[qname[:-2]][qname[-1]] = seq
    print('loaded {0} SV flankings from {1} files'.format(len(SVflks), len(fs)))

    from math import ceil
    print('using {0} threads'.format(nb_cpu))
    step = ceil(len(confident_SVs)/nb_cpu)
    worker_intervals = [[[i*step, (i+1)*step, confident_SVs], is_debug] for i in range(nb_cpu)]
    mydata = {}
    with Pool(nb_cpu) as p:
        rets = p.map(worker, worker_intervals)
    for ret in rets:
        for key in ret:
            mydata[key] = ret[key]
    print('{0} entries total'.format(len(mydata)))
    if '/' in outputname: 
        outputdir = '/'.join(outputname.split('/')[:-1])+'/'
    else:
        outputdir = './'
    with open(outputdir+'dump.pkl', 'wb') as file_out:
        pickle.dump(mydata, file_out)

    header = ['SVname(debug)', 'sampleName', 'type', 'contigName', 'chrom', 'strand','refPos', 'tigPos',
            'SVlength', 'NM_i', 'splicealn_qs', 'splicealn_qe', 'divergence',
            'parentgene_ID', 'parentgene_symbol',
            'whichSideIs5',
            'hasPolyA', 'tail',
            'hasTSD', 'TSD1_left', 'TSD2_right', 
            'lost_5exon', 'lost_3exon', 'lost_5content', 'lost_3content', 'lost_middlecontent',
            'seg_left_flk', 'seg_left_TSD', 'seg_left_inbetween', 'seg_retrocopy', 'seg_right_inbetween', 'seg_right_TSD', 'seg_right_flk']
    with open(outputname, 'w') as file_out:
        file_out.write('\t'.join(header)+'\n')
        for SVname in confident_SVs:
            c = mydata[SVname]
            samplename, tmp = SVname.split('@')
            SVtype, tigname, chrom_,refpos, tigpos = tmp.split('_') 
            strand = chrom_[-1]
            chrom = chrom_[:-1].replace('-', '_')

            # truncation status
            lost_5, lost_3, trunc_5, trunc_3, mid_lost = [], [], [], [], []
            covs = PPG_100flk[SVname]['cov_relative']
            for idx_cov, cov in enumerate(covs):
                left_lost = int(cov[0]==0)
                right_lost = int(cov[-1]==0)
                left_truncation = int((cov[0]<0.9) and (cov[0]>0))
                right_truncation = int((cov[-1]<0.9) and (cov[-1]>0))
                for i1 in range(len(cov)):
                    if cov[i1]!=0: break
                for i2 in range(len(cov)-1, 0, -1):
                    if cov[i2]!=0: break
                cov_mid = cov[i1:i2+1]
                mid_lost.append(int((np.sum(cov_mid<1)>0)))
                if transcriptDirection[PPG_100flk[SVname]['targetname_id'][idx_cov]]=='-':
                    lost_5.append(right_lost)
                    lost_3.append(left_lost)
                    trunc_5.append(right_truncation)
                    trunc_3.append(left_truncation)
                else:
                    lost_5.append(left_lost)
                    lost_3.append(right_lost)
                    trunc_5.append(left_truncation)
                    trunc_3.append(right_truncation)
            if (0 in lost_5): lost_5 = 'FALSE'
            else: lost_5 = 'TRUE'
            if (0 in lost_3): lost_3 = 'FALSE'
            else: lost_3 = 'TRUE'
            if (0 in trunc_5): trunc_5 = 'FALSE'
            else: trunc_5 = 'TRUE'
            if (0 in trunc_3): trunc_3 = 'FALSE'
            else: trunc_3 = 'TRUE'
            if (0 in mid_lost): mid_lost = 'FALSE'
            else: mid_lost='TRUE'
            
            # divergence
            NM_i = c['NM_i']
            cigar = re.findall('[0-9]+M', c['cigar'])
            M = sum([int(_[:-1]) for _ in cigar])
            cigar = re.findall('[0-9]+[ID]', c['cigar'])
            other = sum([int(_[:-1]) for _ in cigar])
            div = (NM_i-other)/M
            
            # formatting
            if not c['has_polyA']:
                if transcriptDirection[c['targetname_id']]==c['aln_strand']:
                    which = 0
                else:
                    which = 1
                has_polyA = 'FALSE'
                tail = '-'
            else:
                which = c['which']
                has_polyA = 'TRUE'
                tail = c['tail']
            if not c['has_TSD']:
                has_TSD = 'FALSE'
                TSD1 = TSD2 = '-'
            else:
                has_TSD='TRUE'
                TSD1, TSD2 = c['TSD']
            if 'segments' not in c:
                segments = ['-']*7
            else:
                segments = c['segments']
                
            line = [SVname, samplename, SVtype, tigname, chrom, strand, refpos, tigpos,
                    c['aln_l'], NM_i, c['aln_qs'], c['aln_qe'], div, 
                    c['targetname_id'], c['targetname_symbol'],
                    which, 
                    has_polyA, tail,
                    has_TSD, TSD1, TSD2,
                    lost_5, lost_3, trunc_5, trunc_3, mid_lost, 
                    *segments
                ]
            file_out.write('\t'.join([str(_) for _ in line])+'\n')
            
