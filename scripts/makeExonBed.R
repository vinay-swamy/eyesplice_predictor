library(tidyverse)
setwd('/Volumes/data/eyesplice_predictor/testing/')

# 
# fromGTF_A3ss <- read_tsv('rmats_out/RPE_Fetal.Tissue/fromGTF.A3SS.txt') %>% select(-contains('flanking'), -ID) %>% distinct()  
#     
# fromGTF_A5ss <- read_tsv('rmats_out/RPE_Fetal.Tissue/fromGTF.A5SS.txt') %>% select(-contains('flanking'), -ID) %>% distinct() 
# 
# #we'll try the window size at 40
# ws=30
# short_to_long <- rbind(filter(fromGTF_A5ss, strand == '+'),
#                        filter(fromGTF_A3ss, strand == '-')  
#                         ) %>% mutate(us_length=shortEE-shortES, ds_length=longExonEnd -shortEE) %>% 
#     filter(us_length>=ws, ds_length>=ws)
# long_to_short <- rbind(filter(fromGTF_A5ss, strand == '-'),
#                        filter(fromGTF_A3ss, strand == '+')  
#                         ) %>% mutate(us_length=shortEE-shortES, ds_length= shortES-longExonStart_0base) %>% 
#     filter(us_length>=ws, ds_length>=ws)
# nrow(short_to_long)+nrow(long_to_short)
# #load in stringtie gtf > remove all all above exons > search for exons getting longer shorter > remove those detected in rMATs
# rmats_gtf <- rtracklayer::readGFF('../ref/tissue_gtfs/RPE_Fetal.Tissue_st.gtf')

ref_gtf <- rtracklayer::readGFF('ref/gencodeAno_comp.gtf')
sample_table <- read_tsv('sampleTableESP.tsv', col_names = c('sample', 'run', 'paired','tissue', 'subtissue', 'origin'))
qfiles <- list.files('quant_files', 'quant.sf', full.names = T, recursive = T)
names <- str_split(qfiles, '/') %>% sapply(function(x)x[2])
txi <- tximport::tximport(files = qfiles, type='salmon', txOut = T, countsFromAbundance = 'lengthScaledTPM')
colnames(txi$counts) <- names
#lets try min 5 counts per sample
counts <- txi$counts %>% as.data.frame() %>%  mutate(ID=rownames(.)) %>% 
    select(ID, sample_table %>% filter(tissue=='RPE_Fetal.Tissue') %>% pull(sample)) %>% {filter(., rowSums(.[,-1]) > ncol(.[,-1])*5 )}
#variables for data collection - windowsize,number of isoforms for a given tx
wsize=40 #lets fix window size, and only keep exons that have 2 isoforms.

end_longer <- filter(ref_gtf, type == 'exon',transcript_id %in% counts$ID) %>% group_by(seqid, strand, start ) %>% 
    summarise(n=length(unique(end)),max_end=max(end), min_end=min(end)) %>% 
    mutate(short_length=min_end-start, long_length=max_end-min_end) %>% 
    filter(n==2, short_length>=wsize, long_length>=wsize) %>% ungroup %>% 
    mutate(wstart=min_end-wsize, wend=min_end+wsize, name=paste0('EL_', 1:nrow(.)), score=1000)
start_longer <- filter(ref_gtf, type == 'exon', transcript_id %in% counts$ID) %>% group_by(seqid, strand, end ) %>% 
    summarise(n=length(unique(start)),max_start=max(start), min_start=min(start)) %>%
    mutate(short_length=end-max_start, long_length=max_start-min_start) %>% 
    filter(n==2,short_length>=wsize, long_length>=wsize ) %>% ungroup %>% 
    mutate(wstart=max_start-wsize, wend=max_start+wsize, name=paste0('SL_', 1:nrow(.)), score=1000 ) 
write_tsv(end_longer, 'end_longer_full.tsv')
write_tsv(start_longer,'start_longer_full.tsv')

end_longer %>% select(seqid, wstart, wend, name, score, strand) %>%  write_tsv('end_longer.bed', col_names = F)
start_longer %>% select(seqid, wstart, wend, name, score, strand) %>%  write_tsv('start_longer.bed', col_names = F)
