#!/usr/bin/env nextflow

def usage(){
    log.info """
    Usage:

    For an example human dataset one would run:
    nextflow run main.nf --reads 'fastqs/*_R{1,2}.fastq.gz' --genome GRCh37 --strand reverse

    Required:
    --reads         Path to input data, in quotes
    -profile        Configuration profile
    --fasta         Folder containing reference genome .fasta
    --gtf           Folder containing reference genome .gtf 
    
    additional:
    --strand        reverse, forward or unstranded (default)
    --outdir        Output directory
    --star_index    Pre-computed STAR index file
    --pair          paired or single end reads
    """
}

/*
Step 0: SETUP

Here are some initial helper functions. They are here to parse command line
arguments and setup up our initial channels for:
    -Genome
    -Annotation
    -Paired or single end read
    -STAR index (or create one)
*/

if(params.help){
    usage()
    exit 0
}

clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.tp_clip_r1
three_prime_clip_r2 = params.tp_clip_r2

if(params.pair != "paired" && params.pair != "single"){
    exit 1, "Invalid option for reads: ${params.pair}. It must be 'single' or 'paired'"
}

if(params.reads && params.pair == "paired"){
    log.info """Reading paired end reads"""
    reads = Channel.fromFilePairs(params.reads).ifEmpty{ exit 1, "params.reads was empty - no input files supplied"}
    .into{raw_fastqc; raw_trimgalore}
}

if(params.reads && params.pair == "single"){
    log.info """Reading single end reads"""
    raw_reads_trimgalore = Channel.fromPath(params.reads).ifEmpty { exit 1, "params.readPaths was empty - no input files supplied"}
    .into{raw_fastqc; raw_trimgalore}
}

if(params.gtf){
    Channel.fromPath(params.gtf).ifEmpty{exit 1, "GTF file not found in ${params.gtf}"}
    .into{gen_gtf; idx_gtf}
}else{
    exit 1, "GTF argument not supplied"
}

if(params.fasta){
    Channel.fromPath(params.fasta).ifEmpty{exit 1, "GTF file not found in ${params.fasta}"}.
    set{idx_fasta}
}else{
    exit 1, "FASTA argument not supplied"
}


if(params.star_index){
    star_index = Channel.fromPath(params.star_index).ifEmpty{exit 1, "STAR Index not found in ${params.star_index}. \
    Please select the correct path or leave blank and let the tool create an index using your genome."} 
}else{
    log.info"""Creating index"""
    process make_STAR_index{
        label 'multithreaded'
        publishDir path: { params.save_ref ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.save_ref ? it : null }, mode: 'copy'

        input:
        file fasta from idx_fasta
        file gtf from idx_gtf

        output:
        file "star" into star_index

        script:
        """
        mkdir star
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN ${task.cpus} \\
            --sjdbGTFfile $gtf \\
            --genomeDir star/ \\
            --genomeFastaFiles $fasta
        """
    }
}

/*
Step 1: QUALITY CONTROL

Run FastQC on our files and trim
*/


process fastqc{
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    file r from raw_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $r
    """
}



process trim{
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    file r from raw_trimgalore

    output:
    file "*fq.gz" into trimmed_reads
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports


    script:
    
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
    
    if (params.pair == "single") {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $r
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $r
        """
    }
}


/*
Step 2: ALIGNMENT

Run alignment using STAR
*/

process star{
    label 'multithreaded'
    publishDir "${params.outdir}/STAR", mode: 'copy'

    input:
    file reads from trimmed_reads
    file index from star_index
    file gtf from gen_gtf

    output:
    set file("*Log.final.out"), file ('*.bam') into star_aligned
    file "*.out" into alignment_logs
    file "*SJ.out.tab"
    file "*Log.out" into star_log
    file "${prefix}Aligned.sortedByCoord.out.bam.bai" into bam_index_rseqc, bam_index_genebody
    file "*.raw" into raw_counts

    script:
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/

    """
    STAR --genomeDir $index \\
         --sjdbGTFfile $gtf \\
         --readFilesIn $reads  \\
         --runThreadN ${task.cpus} \\
         --twopassMode Basic \\
         --outWigType bedGraph \\
         --outSAMtype BAM SortedByCoordinate $avail_mem \\
         --readFilesCommand zcat \\
         --runDirPerm All_RWX \\
         --outFileNamePrefix $prefix \\
         --quantMode GeneCounts
        
    samtools index ${prefix}Aligned.sortedByCoord.out.bam
    """
}

/*
Step 3: VARIANT CALLING

Run variant calling pipeline, using Opossum and
Platypus
*/

process opossum{
    publishDir "${params.outdir}/STAR", mode: 'copy'

    input:
    file bamfile from star_aligned

    output:
    file "${bamfile.baseName}_opossum.bam" into platypus_in

    script:
    """
    opossum --BamFile bamfile --OutFile "${bamfile.baseName}_opossum.bam
    """

}

process platypus{
    publishDir "${params.outdir/platypus}", mode: 'copy'

    input:
    file bamfile from platypus_in

    output:
    file "${bamfile.baseName}.vcf" into platypus_vcf

    script:
    """
    platypus callVariants --bamFiles $bamfile --refFile $reference --filterDuplicates 0 --minMapQual 0 \
	--minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o ${bamfile.baseName}.vcf
    """
}

/*
Step 4: HLA TYPING

If human then run HLA type prediction
*/

/*
Step 5: IDENTIFY COMMON VARIANTS

find overlapping variants in the the bed files supplied 
(if there's more than one). Then also remove those variants
which are in the alt reference?

Also annotate the variants using ANNOVAR
*/

process common_variants {
    publishDir "${params.outdir/variants}", mode: 'copy'

    input:
    vcfs from platypus_vcf.collect()
    no_files from platypus_vcf.size()

    output:
    "overlapping_variants.vcf" into variants

    script:
    """
    bcftools isec -n ~$no_files $vcfs -O v -o overlapping_variants.vcf
    """
}

/*
process annovar {
    publishDir "${params.outdir/variants}", mode: 'copy'

    input:
    file variant from variants

    script:
    """
    """
}
*/