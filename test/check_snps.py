# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 13:58:54 2018

@author: nilesh-tawari
"""

"""
Script for sanity check snps against the latest (2.0) release of clinvar vcf files 
"""
import os
import allel
import pandas as pd



# STEP1: read clinvar vcf downloaded from: 
#ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20171231.vcf.gz

clinvar_vcf = 'clinvar_20171231.vcf.gz'
df = allel.vcf_to_dataframe(clinvar_vcf, fields=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'ALLELEID', 'CLNHGVS', 'CLNREVSTAT',  'CLNSIG', 'CLNVC', 'ORIGIN', 'RS', 'SSR'] )

df = df.loc[df['CLNVC'] == 'single_nucleotide_variant'] 
df = df.loc[df['ALT_1'] != 'TGAG']  # wrongly annotated entry in clinvar as snv
df = df.loc[df['ALT_1'] != ''] 
df['IDC'] = df.CHROM.astype(str) + ':' + df.POS.astype(str) + '-' + df.REF.astype(str) + '/' +df.ALT_1.astype(str)

# STEP2: read clinvar singles text file downloaded from https://github.com/macarthur-lab/clinvar/blob/master/output/b37/single/clinvar_alleles.single.b37.tsv.gz
clinvar_si = os.path.join('..', 'output', 'b37', 'single', 'clinvar_alleles.single.b37.tsv.gz')
df_c = pd.read_csv(clinvar_si, sep='\t', comment = '#', chunksize=1000, \
                low_memory=False, iterator = True, compression='gzip')
df_c = pd.concat(list(df_c), ignore_index=True)
# take only SNVs
df_c = df_c.loc[df_c['ref'].isin(['A', 'T', 'G', 'C'])]
df_c = df_c.loc[df_c['alt'].isin(['A', 'T', 'G', 'C'])]
#df_ch = df_c.head(1000)

# STEP3: read clinvar multi text file downloaded from https://github.com/macarthur-lab/clinvar/blob/master/output/b37/multi/clinvar_alleles.multi.b37.tsv.gz
clinvar_m = os.path.join('..', 'output', 'b37', 'multi', 'clinvar_alleles.multi.b37.tsv.gz')
df_m = pd.read_csv(clinvar_m, sep='\t', comment = '#', chunksize=1000, \
                low_memory=False, iterator = True, compression='gzip')
df_m = pd.concat(list(df_m), ignore_index=True)
df_m = df_m.loc[df_m['ref'].isin(['A', 'T', 'G', 'C'])]
df_m = df_m.loc[df_m['alt'].isin(['A', 'T', 'G', 'C'])]

# STEP4: merge clinvar text files
df_mc = pd.concat([df_c, df_m], ignore_index=True)
df_mc['IDC'] = df_mc.chrom.astype(str) + ':' + df_mc.pos.astype(str) + '-' + df_mc.ref.astype(str) + '/' +df_mc.alt.astype(str)

# find uniq variantsin clinvar not present in MC text files 
df_cl_uniq = df.loc[~df['IDC'].isin(list(df_mc.IDC))]
df_cl_uniq.reset_index(drop=True, inplace=True)
# get associated refseq ids for these variants
df_cl_uniq['REFIDNC'] = df_cl_uniq.CLNHGVS.str.split('.', expand=True)[0]
not_included = set(df_cl_uniq['REFIDNC'])

df_cl_uniq.to_csv('snps_not_in_tsv_files.txt', sep='\t', index=False)

# conclusion: based on manual search of these variants seems to be new additions to clinvar 
