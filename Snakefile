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
def tissue_to_gtf(tissue, sample_dict):
    res=[]
    for sample in sample_dict.keys():
        if sample_dict[sample]['tissue']==tissue :
            res.append('st_out/{}.gtf'.format(sample))
    return (res)

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
win_size=config['window_size']

sample_file=config['sampleFile']
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
sample_names=sample_dict.keys()
gtf=config['gtf']
genome=config['genome']
tx_fasta=config['tx_fasta']
fql=config['fastq_path']
bam_path=config['bam_path']
subtissues=['RPE_Fetal.Tissue', 'synth']
rmats_events=['SE','RI','MXE','A5SS','A3SS']
rule all:
    input: expand('quant_files/{sampleID}/quant.sf', sampleID=sample_names),\
    'data/all_exons.bed', expand('data/cleaned/{sample}_bp_features.tsv', sample=sample_names)

rule build_STARindex:
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
rule index_bams:
    input: bam_path+'STARbams/{id}/Sorted.out.bam'
    output: bam_path+'STARbams/{id}/Sorted.out.bam.bai'
    shell:
        '''
        module load {samtools_version}
        samtools index -b {input}
        '''

rule calculate_cov:
    input:bam_path+'STARbams/{id}/Sorted.out.bam', bam_path+'STARbams/{id}/Sorted.out.bam.bai'
    output: 'coverage_files/{id}.per-base.bed.gz'
    shell:
        '''
        module load mosdepth
        sample={wildcards.id}
        mosdepth coverage_files/$sample {input[0]}
        '''
rule build_salmon_index:
    input:  tx_fasta
    output: directory('ref/salmonindex_st')
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
rule makeExonBeds:
    input: expand('quant_files/{sampleID}/quant.sf', sampleID=sample_names)
    output: expand('data/{type}_{direction}_full_tab.tsv',type=['grow', 'ref'],direction=['end', 'start']),\
      expand('data/{type}_{direction}_longer.bed',type=['grow', 'ref'], direction=['end', 'start'] )
    shell:
        '''
        module load {R_version}
        Rscript scripts/makeExonBed.R {working_dir} {sample_file} {gtf} {win_size} {output}
        '''
rule mergeBeds:
    input: expand('data/{type}_{direction}_longer.bed',type=['grow','ref'], direction=['end', 'start'] )
    output:'data/all_exons.bed'
    shell:
        '''
        module load {bedtools_version}
        bash scripts/merge_beds_distinct.sh {input} {output}
        '''

rule intersect_coverage:
    input:'data/all_exons.bed', 'coverage_files/{sample}.per-base.bed.gz'
    output:'coverage_files/{sample}_exon_cov.bed'
    shell:
        '''
        module load bedtools
        cut -f1,2,3,4 {input[0]} |
        bedtools intersect -loj -a {input[1]} -b stdin | awk ' $6 != "-1"' - >  {output}
        '''

rule spread_coverage:
    input:'coverage_files/{sample}_exon_cov.bed'
    output: 'data/cleaned/{sample}_bp_features.tsv'
    shell:
        '''
        python3 scripts/makePerBaseFeatureTable.py {working_dir} {input} {output}
        '''
