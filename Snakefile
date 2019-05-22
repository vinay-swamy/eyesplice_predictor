def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'files':info[1].split(','),'paired':True if info[2]=='y' else False, 'tissue':info[3],'subtissue':info[4]}
    return(res)
def salmon_input(id,sample_dict,fql):
    paired=sample_dict[id]['paired']
    id= fql + 'fastq_files/' + id
    if paired:
        return('-1 {s}_1.fastq.gz -2 {s}_2.fastq.gz'.format(s=id))
    else:
        return('-r {}.fastq.gz'.format(id))
#software versioning
salmon_version=config['salmon_version']
stringtie_version=config['stringtie_version']
STAR_version=config['STAR_version']
rmats_version=config['rmats_verson']
R_version=config['R_version']
TransDecoder_version=config['TransDecoder_version']
samtools_version=config['samtools_version']
gffcompare_version=config['gffcompare_version']
hmmer_version=config['hmmer_version']
crossmap_version=config['crossmap_version']
deeptools_version=config['deeptools_version']
mosdepth_version=config['mosdepth_version']
bedtools_version=config['bedtools_version']
working_dir=config['working_dir']


sample_file=config['sampleFile']
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
sample_names=sample_dict.keys()
gtf=config['gtf']
genome=config['genome']
tx_fasta=config['tx_fasta']
fql=config['fastq_path']
bam_path=config['bam_path']
subtissues=['RPE_Fetal_PE', 'synth']
rmats_events=['SE','RI','MXE','A5SS','A3SS']
rule all:
    input: expand('quant_files/{sampleID}/quant.sf', sampleID=sample_names), expand('st_out/{sampleID}.gtf', sampleID=sample_names),\
    expand('rmats_out/{tissue}/{event}.MATS.JC.txt', event=rmats_events, tissue=subtissues)

rule build_STAR_index:
    input:genome, gtf
    output:directory('ref/STARindex')
    shell:
        '''
        module load {STAR_version}
        mkdir -p ref/STARindex
        STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} --sjdbOverhang 100
        '''
rule align_STAR:
    input: fastqs=lambda wildcards: [fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.id),fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.id)] if sample_dict[wildcards.id]['paired'] else fql+'fastq_files/{}.fastq.gz'.format(wildcards.id),
        index='ref/STARindex'
    output: temp(bam_path+'STARbams/{id}/raw.Aligned.out.bam'), bam_path + 'STARbams/{id}/raw.Log.final.out'
    shell:
        '''
        id={wildcards.id}
        mkdir -p {bam_path}/STARbams/$id
        module load {STAR_version}
        STAR --runThreadN 8 --genomeDir {input.index} --outSAMstrandField intronMotif  --readFilesIn {input.fastqs} \
        --readFilesCommand gunzip -c --outFileNamePrefix {bam_path}/STARbams/$id/raw. --outSAMtype BAM Unsorted
        '''

rule sort_bams:
    input:bam_path+'STARbams/{id}/raw.Aligned.out.bam'
    output:bam_path+'STARbams/{id}/Sorted.out.bam'
    shell:
        '''
        module load {samtools_version}
        samtools sort -o {output[0]} --threads 7 {input[0]}
        '''
rule run_stringtie:
    input: bam_path+'STARbams/{sample}/Sorted.out.bam'
    output:'st_out/{sample}.gtf'
    shell:
        '''
        module load {stringtie_version}
        stringtie {input[0]} -o {output[0]} -p 8 -G ref/gencodeAno_bsc.gtf
        '''
rule preprMats_running:
    input: expand(bam_path+'STARbams/{id}/Sorted.out.bam',id=sample_names)
    params: bam_dir=bam_path + 'STARbams/'
    output:expand('ref/rmats_locs/{tissue}.rmats.txt',tissue=subtissues)
    shell:
        #include trailing / for bam_dir
        '''
        bam_path={bam_path}/STARbams/
        suff=/Sorted.out.bam
        grep synth {sample_file} | cut -f1 | while read p; do echo "$bam_path/$p/$suff"; done  > ref/rmats_locs/synth.rmats.txt
        grep -v synth {sample_file} | cut -f1 | while read p; do echo "$bam_path/$p/$suff"; done  > ref/rmats_locs/RPE_Fetal_PE.rmats.txt
        '''

rule runrMATS:
    input: 'ref/rmats_locs/{tissue}.rmats.txt','ref/STARindex',gtf
    output:expand('rmats_out/{{tissue}}/{event}.MATS.JC.txt', event=rmats_events)
    # might have to change read length to some sort of function
    shell:
        '''
        tissue={wildcards.tissue}
        module load {rmats_version}
        rmats --b1 {input[0]} --b2 ref/rmats_locs/synth.rmats.txt  -t paired  \
        --nthread 8  --readLength 130 --gtf {input[2]} --bi {input[1]} --od rmats_out/$tissue
        '''
rule build_salmon_index:
    input:  tx_fasta
    output:'ref/salmonindex_st'
    shell:
        '''
        module load {salmon_version}
        salmon index -t {input} --gencode -i {output} --type quasi --perfectHash -k 31
        '''

rule run_salmon:
    input: fastqs=lambda wildcards: [fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sample_dict[wildcards.sampleID]['paired'] else fql+'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
        index='ref/salmonindex_st'
    params: cmd=lambda wildcards: salmon_input(wildcards.sampleID,sample_dict,fql)
    output: 'quant_files/{sampleID}/quant.sf'
    shell:
        '''
        id={wildcards.sampleID}
        module load {salmon_version}
        salmon quant -p 4 -i {input.index} -l A --gcBias --seqBias  {params.cmd} -o quant_files/$id
        '''
