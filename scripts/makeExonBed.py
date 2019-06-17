import pybedtools
import pandas as pd
import sys
import os
args=sys.argv
#wd=args[1]
#ref_gtf_file=args[2]
#wsize=args[3]  
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
        gtf_header=['seqid','source','type','start','stop','score','strand','frame']
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
    col=[cn for cn in col if cn not in gtf_header]
    df=df[gtf_header+col]
    return(df)

#not going to worray about account for tpm exp, but just gonna have a exp threshold in prep step
os.chdir('/Users/vinayswamy/NIH/eyesplice_predictor')
ref_gtf_file='ref/gencodeAno_comp.gtf'
ref_gtf=read_GTF(ref_gtf_file)
end_longer_all=ref_gtf[ ref_gtf['type'] == 'exon'].group


