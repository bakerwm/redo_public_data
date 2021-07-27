#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract records from GTF file, save as BED6 format
filter:
  - feature: gene|CDS|exon|...
  - gene_biotype: protein_coding|...
  
Example:
$ python filt_gtf.py -i ensembl.GRCh37.75.gtf.gz -o to_PolyASite_full -b pas.2.0.hg19.bed -f gene -x --size 5000 --flank 6000
$ python filt_gtf.py -i ensembl.GRCh37.75.gtf.gz -o to_PolyA_DB_full -b human.PAS.bed -f gene -x --size 5000 --flank 6000
  
"""


import os
import io
import pathlib
import argparse
import pybedtools
from xopen import xopen


def parse_desc(s, key='gene_name'):
    """
    Description in GTF:
    gene:
    gene_id "ENSG00000227232.4"; transcript_id "ENSG00000227232.4"; 
    gene_type "pseudogene"; gene_status "KNOWN"; gene_name "WASH7P"; 
    transcript_type "pseudogene"; transcript_status "KNOWN"; 
    transcript_name "WASH7P"; level 2; havana_gene "OTTHUMG00000000958.1";
    """
    d = {}
    if isinstance(s, str):
        for p in s.strip().split(';'):
            if " " in p:
                a, b = p.split()[:2]
                b = b.replace('"', '')
                d.update({a:b})
    # fix gene_type (GENCODE), gene_biotype (ENSEMBL)
    if key == 'gene_type' and 'gene_biotype' in s:
        key = 'gene_biotype'
    elif key == 'gene_biotype' and 'gene_type' in s:
        key = 'gene_type'
    return d.get(key, None)
        

def read_gtf(fh, **kwargs):
    """
    Arguments
    ---------
    feature : str
        gene|exon|CDS|transcript
    gene_biotype : str
        TEC|protein_coding
    """
    feature = kwargs.get('feature', 'gene')
    gene_biotype = kwargs.get('gene_biotype', True)
    for l in fh:
        p = l.strip().split('\t')
        # filtering
        if len(p) < 9:
            continue
        if isinstance(feature, str) and not p[2] == feature:
            continue
        gb = parse_desc(p[8], 'gene_biotype')
        if isinstance(gene_biotype, str) and not gb == gene_biotype:
            continue
        yield p


def gtf_to_bed(gtf, bed, **kwargs):
    print('[1/5] convert gtf to bed')
    with xopen(gtf) as r, open(bed, 'wt') as w:
        for g in read_gtf(r, **kwargs):
            g_chr = g[0] if g[0].startswith('chr') else 'chr'+g[0]
            g_start = str(int(g[3])-1)
            g_name = parse_desc(g[8], 'gene_name')
            gb = parse_desc(g[8], 'gene_biotype')
            bed6 = [g_chr, g_start, g[4], g_name, '254', g[6]]
            if kwargs.get('add_type', False):
                bed6.append(gb)
            bed6 = list(map(str, bed6))
            w.write('\t'.join(bed6)+'\n')


def filt_by_flank(a, b, flank=6000):
    """
    Filt the BED file by size
    1. flanking BED (6kb) not including genes
    """
    print('[2/5] filter by flanking={}'.format(flank))
    # sort
    ba = pybedtools.BedTool(a).sort()
#     ba = pybedtools.BedTool(a)#.filter(lambda x: len(x) > size)
#     ba = ba.sort() # sort
    # filter by flank region
    with open(b, 'wt') as w:
        pre_gene = [None, None] #fwd, rev
        for i in ba:
            if i.strand == '+':
                if pre_gene[0] is None:
                    pre_gene[0] = i
                else:
                    if pre_gene[0].chrom == i.chrom and i.start - pre_gene[0].end <= flank:
                        continue
                    else:
                        pre_gene[0] = i
            elif i.strand == '-':
                if pre_gene[1] is None:
                    pre_gene[1] = i
                else:
                    if pre_gene[1].chrom == i.chrom and i.start - pre_gene[1].end <= flank:
                        continue
                    else:
                        pre_gene[1] = i
            bed6 = list(map(str, i[:7])) # bed6+1
            w.write('\t'.join(bed6)+'\n')
            

def filt_by_pas(a, b, out):
    """
    Filter by another BED (pas)
    extract the distal pas site: tss-pas
    """
    print('[3/5] filt by PAS: {}'.format(b))
    ab = pybedtools.BedTool(a).intersect(b, s=True, wo=True)
    ab = ab.sort() #
    # by pas at 3' end
    with open(out, 'wt') as w:
        pas_pos = None
        pre_bed = None
        for i in ab:
            i_pos = [int(i[8]), int(i[9])] # depends on a (bed6+x)
            if pre_bed is None: # init
                pas_pos = i_pos
                pre_bed = i
            else:
                if i[3] == pre_bed[3]: # update pas_pos
                    if i[5] == '-' and pas_pos[0] < i_pos[0]:
                        pas_pos = i_pos
                    elif i[5] == '+' and pas_pos[1] > i_pos[1]:
                        pas_pos = i_pos
                    else:
                        pass
                else: # save record
                    bed6 = list(map(str, i[:7]))
                    if i[5] == '-':
                        bed6[1] = str(pas_pos[0])
                    else:
                        bed6[2] = str(pas_pos[1])
                    if int(bed6[1]) >= int(bed6[2]):
                        continue # skip
                    # save to bed6
                    w.write('\t'.join(bed6)+'\n')
                    pas_pos = i_pos
                    pre_bed = i
        # last record
        bed6 = list(map(str, i[:6]))
        if i[5] == '-':
            bed6[1] = str(pas_pos[0])
        else:
            bed6[2] = str(pas_pos[1])
        if int(bed6[1]) < int(bed6[2]):
            w.write('\t'.join(bed6)+'\n')

    
def extract_bed(a, b, key='protein_coding', invert_match=False):
    """
    extract bed by column-7: protein_coding
    
    invert_match : bool
        select non-matching lines
    """
    print('[4/5] extract pc')
    with open(a) as r, open(b, 'wt') as w:
        for line in r:
            if line.strip().endswith(key):
                if not invert_match:
                    w.write(line)
            else:
                if invert_match:
                    w.write(line)
        
        
def filt_by_size(a, b, size=5000):
    print('[5/5] filter by size={}'.format(size))
    bed = pybedtools.BedTool(a).filter(lambda x: len(x) > size)
    bed.saveas(b)
    

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='The gene GTF file')
    parser.add_argument('-o', required=True, help='The outdir')
    parser.add_argument('-b', required=True, help='The PAS file in BED format')
    parser.add_argument('-f', '--feature', default='gene', help='feature, default [gene]')
    parser.add_argument('-g', '--gene-biotype', dest='gene_biotype', default=True, help='The gene_biotype for the record')
    parser.add_argument('-x', '--add-type', dest='add_type', action='store_true', help='Add gene_biotype to column-7 in BED')
    parser.add_argument('--size', type=int, default=0, help='The minimum size of BED feature, default [0]')
    parser.add_argument('--flank', type=int, default=0, help='The minimum distance to neighbour gene, default [0]')
    return parser


def main():
    args = vars(get_args().parse_args())
    gtf = args.pop('i', None)
    outdir = args.pop('o', None)
    bed_filt = args.pop('b', None)
    size = args.pop('size', 0)
    flank = args.pop('flank', 0)
    # prepare outdir
    if not isinstance(outdir, str):
        outdir = pathlib.Path.cwd().name
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # prepare files
    prefix = os.path.basename(gtf)
    prefix = prefix.replace('.gtf.gz', '') if prefix.endswith('.gz') else prefix.replace('.gtf', '')
    bed_gene = os.path.join(outdir, prefix+'.gene.bed')
    bed_flank = os.path.join(outdir, prefix+'.gene.f{:d}k.bed'.format(int(flank/1000)))
    bed_pas = bed_flank.replace('.bed', '.pas.bed')
    bed_pc = bed_pas.replace('.gene.', '.pc.') # protein_coding
    bed_size = bed_pc.replace('.bed', '.{}k.bed'.format(int(size/1000)))
    # step1. extract genes from gtf
    gtf_to_bed(gtf, bed_gene, feature=args['feature'], gene_biotype=args['gene_biotype'], add_type=args['add_type'])
    # step2. filt by size
    filt_by_flank(bed_gene, bed_flank, flank=flank)
    # step3. filt by another BED file
    filt_by_pas(bed_flank, bed_filt, bed_pas)
    # step5. extract protein_coding genes
    extract_bed(bed_pas, bed_pc)
    # step4. filt by size
    filt_by_size(bed_pc, bed_size)
    print('Finish, save to file: {}'.format(bed_pas))

if __name__ == '__main__':
    main()
    
# EOF