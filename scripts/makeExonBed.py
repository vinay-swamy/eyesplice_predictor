
# coding: utf-8

# In[32]:


from pybedtools import BedTool
import pandas as pd
import sys
import os
import numpy as np
import random as rd
args=sys.argv
wd=args[1]
ref_gtf_file=args[2]
sample_file=args[3]
tx_quant_file=args[4]
subtissue=args[5]
wsize=int(args[6])
out_bed_file=args[7]
#love stealing my own code 
def read_GTF(file):
    def file_len(fname):
        with open(fname) as f:
            os=0
            for i, l in enumerate(f):
                if l[0]=='#':
                    os+=1
        return ((i + 1 -os , os))

    def parse_attributes(st):
        att_list=st.replace('"','').replace('=', ' ').split(';')[:-1]
        out_dict=dict()
        for att in att_list:
            att=att.lstrip(' ').split(' ')
            out_dict[att[0]]=att[1]
        return(out_dict)
    file_length, os=file_len(file)
    outdf=[None]*file_length
    with open(file) as gtf:
        gtf_header=['seqid','source','type','start','end','score','strand','frame']
        for k in range(os):
            gtf.readline().strip('\n')
        for i in range(file_length):
            line= gtf.readline().strip('\n')
            line_list=line.split('\t')
            line_dict=dict()
            [line_dict.update({value:line_list[j]}) for j,value in enumerate(gtf_header)]
            line_dict.update(parse_attributes(line_list[-1]))
            outdf[i]=line_dict
    df=[k for k in outdf if k is not None]
    df=pd.DataFrame(outdf)
    col=df.columns.tolist()
    df['start']=pd.to_numeric(df['start'])
    df['end']=pd.to_numeric(df['end'])
    col=[cn for cn in col if cn not in gtf_header]
    df=df[gtf_header+col]
    return(df)

def bool_to_split(dfs):
    rd.seed(8774)
    l=[l for l in range(dfs)]
    idx=set(rd.sample(l, int(dfs/2)))
    idx_bool=[True if i in idx else False for i in l]
    return(idx_bool)
'''
Overall process 
    -first, find all exons that have an alternative 3'(end longer) or alternative 5'(start_longer)
    -create a window where the junciton between the long and short forms is the middle of the window
    -

'''


os.chdir(wd)
#os.chdir('/Users/vinayswamy/NIH/eyesplice_predictor')
# wsize=40
# ref_gtf_file='ref/gencodeAno_comp.gtf'
# sample_file='sampleTableESP.tsv'
# tx_quant_file='data/GC_pc_tx_quant.tsv.gz'
# subtissue='RPE_Fetal.Tissue'
# out_bed_file='testing/new_script_old_txquant.bed'

sample_table=pd.read_csv(sample_file, sep='\t', names=['sample', 'run', 'paired', 'tissue', 'subtissue','origin']).query('subtissue == @subtissue')
tx_quant=pd.read_csv(tx_quant_file, sep='\t').loc[:,['transcript_id']+list(sample_table['sample'])]

ref_gtf=read_GTF(ref_gtf_file)
#ref_gtf=ref_gtf.head(1000)
SL_score=111
EL_score=222
R_SL_score=333
R_EL_score=444
I_score=555
E_score=666


# In[33]:


tx_quant.shape
ref_gtf.shape


# In[34]:


'''
A different bed file will be created for each subtissue type, were are only going to build the window bed file from transcripts expressed in the target tissue

'''
keep=(tx_quant
      .iloc[:,1:]
      .sum(axis=1) > len(sample_table.index)  
    )
tx_quant=tx_quant[keep]
ref_gtf=ref_gtf[ref_gtf['transcript_id'].isin(tx_quant["transcript_id"])]


# In[35]:



'''
find all exons that that start on the same coordinate, but end on different coordinates. Do the same but end same and start changes. 
Right now for simplicity sake, Im only selecting alt spliced exons that come in 2 versions
out format for these 2 is BED

'''

end_longer_all=(ref_gtf
                .query('type == "exon"')
                .groupby(['seqid', 'strand', 'start'], as_index=False)
                .agg({ 'end':['min', 'max', 'nunique']})
                .reset_index(drop=True)
               )
end_longer_all.columns=['seqid', 'strand', 'start', 'min_end', 'max_end', 'count']
end_longer_all=(end_longer_all
                .assign(short_length= lambda x: x['min_end'] -x['start'], 
                        long_length= lambda x: x['max_end'] - x['min_end'])
                .rename(columns={'nunique':'count'})
                .query('count > 1')
               )
end_longer= (end_longer_all
             .query('count == 2 & short_length >= @wsize & long_length >=@wsize')
             .assign(wstart= lambda x: x['min_end']-wsize, 
                     wend= lambda x: x['min_end']+wsize,
                     name = lambda x: ['EL_' + str(i) for i in range(len(x.index))])
             .reset_index(drop=True)
             .assign(score=EL_score)
             .loc[:,['seqid', 'wstart', 'wend', 'name', 'score', 'strand']]
            )

start_longer_all= (ref_gtf
                   .query('type == "exon"')
                   .groupby(['seqid','strand', 'end'], as_index=False)
                   .agg({'start':['min', 'max', 'nunique']})
                   .reset_index(drop=True)
                  )
start_longer_all.columns=['seqid', 'strand', 'end', 'min_start', 'max_start','count' ]
start_longer_all=(start_longer_all
                  .assign(short_length= lambda x: x['end'] -x['max_start'], 
                          long_length= lambda x: x['max_start']-x['min_start'])
                  .query('count > 1')
                 )
start_longer= (start_longer_all
               .query("count == 2 & short_length>= @wsize & long_length >= @wsize")
               .assign(wstart= lambda x: x['max_start']- wsize, 
                       wend= lambda x: x['max_start'] + wsize , 
                       name = lambda x: ['SL_' + str(i) for i in range(len(x.index))])
               .reset_index(drop=True)
               .assign(score=SL_score)
               .loc[:,['seqid', 'wstart', 'wend', 'name', 'score', 'strand']]
              )
alt_spliced_windows=pd.concat([end_longer, start_longer]).reset_index(drop=True)


# In[6]:


#********BUG anti_joins are not in pandas (SAD) and the method I'm using can sometimes coerce int64 to floats, which makes bedtools wig out **********************# 
'''
remove any exons got used above using for start/end longer (hence the above problem with antijoins)
antijoin again with bedtools to be real sure 

every exon in below table doesn't get longer/shorter and doesn't overlap with any regions in table above

'''


exon_no_change= (ref_gtf #first, anti_join out known exon locations.
                   .query('type == "exon"').loc[:,['seqid', 'strand', 'start', 'end']]
                   .merge(end_longer_all.loc[:,['seqid','strand', 'start']],
                          how='outer',indicator=True)
                   .query('_merge ==  "left_only"')
                   .drop(columns=['_merge'])
                   .assign(start= lambda x: x['start'].astype(np.int64),
                            end= lambda x: x['end'].astype(np.int64))
                   .merge(start_longer_all.loc[:,['seqid','strand','end']], 
                          how='outer', indicator=True)
                   .query('_merge == "left_only"')
                   .drop(columns=['_merge'])
                   .assign(start= lambda x: x['start'].astype(np.int64), 
                           end= lambda x: x['end'].astype(np.int64),
                           name= 'nc', score=I_score)
                   .loc[:,['seqid', 'start', 'end', 'name', 'score', 'strand']]
                   # I'm concerned that 2 exons might overlap but not share the same start/end so
                   # so also going to subtract out any previous windows.
                   # make sure to keep only exons that are long enough for our window size.
                   .pipe(BedTool.from_dataframe)
                   .intersect(alt_spliced_windows.pipe(BedTool.from_dataframe),
                            s=True, loj=True)
                   .to_dataframe()
                   .rename(columns={'chrom':'seqid'})
                   .assign(length= lambda x: x['end']-x['start'])
                   .query('length >=@wsize & thickEnd == -1 ' )
                   .reset_index(drop=True)
                   .loc[:,['seqid', 'start', 'end', 'name', 'score', 'strand']]
                )


# In[7]:


'''
split the no change exons into 2 tables, one to use to make the no-change exon junction class, and one to use to make exon/intron classes

'''
dfs=len(exon_no_change.index)
spl=np.array(bool_to_split(dfs))
exons_for_junc=exon_no_change.iloc[spl].reset_index(drop=True)
#junction set first 
spl_junc=np.array(bool_to_split(len(exons_for_junc.index)))
junc_sl=(exons_for_junc
         .iloc[spl_junc]
         .reset_index(drop=True)
         .assign(wstart =lambda x: x['start'] -wsize, 
                wend=lambda x: x['start']+wsize, 
                name=lambda x: ['ref_SL_'+str(i) for i in range(len(x.index))],
                score=R_SL_score)
         .loc[:,['seqid', 'wstart', 'wend', 'name', 'score', 'strand']]
        )

junc_el=(exons_for_junc
         .iloc[~spl_junc]
         .reset_index(drop=True)
         .assign(wstart=lambda x: x['end'] -wsize, 
                 wend=lambda x: x['end']+wsize, 
                 name =lambda x: ['ref_EL_' + str(i) for i in range(len(x.index))],
                 score=R_EL_score)
         .loc[:,['seqid', 'wstart', 'wend', 'name', 'score', 'strand']]
        )


ref_altSplice_windows=pd.concat([junc_el, junc_sl]).reset_index(drop=True)

#exon_no_change for sure does not overlap at all with with alt_spliced_windows, but it may overlap with ref_alt_windows, so going to anti_join that bish again.
exons_for_exon=(exon_no_change
                .iloc[~spl]
                .reset_index(drop=True)
                .pipe(BedTool.from_dataframe)
                .intersect(ref_altSplice_windows.pipe(BedTool.from_dataframe),
                          s=True, loj=True)
                .to_dataframe()
                .rename(columns={'chrom':'seqid'})
                .query('thickEnd == -1 ')
                .loc[:,['seqid', 'start', 'end', 'name', 'score', 'strand']]
               )
exons_for_exon.head()


# In[8]:


'''
Now need to get regions of intronic coverage
subtract out all known exons from all transcripts, but pad the exons on wither end, such that pad > wsize. This way we are for sure that the intronic region over lap with none of our previous windows
    - cant overlap with alt_splcied exons, bc window < long exon end. but fo the ref ones, the pad takes care of it.
'''
pad=wsize+5
ref_txs=(ref_gtf
         .query('type == "transcript"')
         .loc[:,['seqid', 'start','end','transcript_id','score','strand']]
         .assign(score=1000)
         .reset_index(drop=True)
        )

ref_exons=(ref_gtf
           .query('type == "exon"')
           .loc[:,['seqid', 'start','end','transcript_id','score','strand']]
           .assign(score=1000, 
                   start= lambda x: x['start'] -pad, 
                   end = lambda x: x['end'] +pad)
           .reset_index(drop=True)
          )

introns=(ref_txs
         .pipe(BedTool.from_dataframe)
         .subtract(ref_exons.pipe(BedTool.from_dataframe), 
                   s=True)
         .to_dataframe()
         .assign(length=lambda x: x['end']-x['start'])
         .query('length >= @wsize*4')# lets make sure were picking from bigger introns
         .rename(columns={'chrom':'seqid'})
         .reset_index(drop=True)
        )


# In[9]:


def make_random_window(ser, ws=wsize *2):
    start=ser[1]
    max_wend=ser[2]-ws
    wstart=rd.randint(start, max_wend)
    wend=wstart+ws
    return((wstart, wend))
intron_windows=(introns
                .apply(make_random_window,axis=1, result_type='expand')
                .rename(columns={0:'wstart', 1:'wend'})
                .merge(introns,how='left', left_index=True, right_index=True)
                .loc[:,['seqid', 'wstart', 'wend', 'name', 'score', 'strand']]
                .assign(name = lambda x: ['intron_' + str(i) for i in range(len(x.index))],
                        score=I_score)
               )


# In[10]:


#exon for exon does not overlap with either ref_alt or alt, so don't need to worry about that
#also cannot overlap with intron version based on above logic, so don't worry about that either
exon_windows= (exons_for_exon.reset_index(drop=True)
                .assign(name='x', 
                        score=E_score,
                        length=lambda x: x['end'] - x['start'])
                .query('length >= @wsize*2 +1 ') # remove potential windows with
              )
exon_windows= (exon_windows
                .apply(make_random_window, axis=1, result_type='expand')
                .rename(columns={0:'wstart', 1:'wend'})
                .merge(exon_windows, left_index=True, right_index=True)
                .loc[:,['seqid', 'wstart', 'wend', 'name', 'score', 'strand']]
                .assign(name=lambda x: ['exon_'+str(i) for i in range(len(x.index))])
                )


# In[11]:


min_samps=min(alt_spliced_windows.shape[0], ref_altSplice_windows.shape[0], 
              intron_windows.shape[0], exon_windows.shape[0])
complete_bed=(pd.concat([alt_spliced_windows, 
                         ref_altSplice_windows.sample(min_samps), 
                         intron_windows.sample(min_samps), 
                         exon_windows.sample(min_samps)
                        ])
              .reset_index(drop=True)
             )
bad=complete_bed.assign(length= lambda x:x['wend']-x['wstart']).query('length<@wsize*2').shape
if bad[0] >0:
    print('Warning: {} rows are less than minimum window size')


# In[13]:


complete_bed.to_csv(out_bed_file, header=False, index=False, sep='\t')

