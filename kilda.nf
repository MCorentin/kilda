#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Importing processes from the modules:
include { ExtractFastaFromBed }                 from './modules/kiv2_create_kmers_DB.nf'
include { CountKmersRegion }                    from './modules/kiv2_create_kmers_DB.nf'
include { FilterOnOccurence }                   from './modules/kiv2_create_kmers_DB.nf'
include { CountKmersOutsideRegion }             from './modules/kiv2_create_kmers_DB.nf'
include { FilterKmersOccuringOutsideRegion }    from './modules/kiv2_create_kmers_DB.nf'
include { RemoveCommonKmers }                   from './modules/kiv2_create_kmers_DB.nf'
include { OutputFasta }                         from './modules/kiv2_create_kmers_DB.nf'

include { CreateFastaKmers }   from './modules/kiv2_counts.nf'
include { CountKmers }         from './modules/kiv2_counts.nf'
include { DumpKmers }          from './modules/kiv2_counts.nf'
include { CreateSampleMap }    from './modules/kiv2_counts.nf'
include { BamToFastq}          from './modules/kiv2_counts.nf'
include { kilda }              from './modules/kiv2_counts.nf'

// Checking mandatory options:
if(params.kilda_dir == "")    error("ERROR: The path to KILDA directory must be provided (see config: 'kilda_dir')")
if(params.kmer_size < 1)      error("ERROR: The kmer size must be > 0 (see config: 'kmer_size')")

// Filling default working folder if not set:
if(params.wdir == "") {
    println("params.wdir is not set, defaulting to './'")
    wdir = Channel.fromPath("./")
} else {
    wdir = Channel.fromPath(params.wdir)
}

// Initialising kmer DB options
outdir        = params.kmer_DB.outdir       ? Channel.fromPath(params.kmer_DB.outdir)       : Channel.fromPath("./kmers_DB/")
genome_fasta  = params.kmer_DB.genome_fasta ? Channel.fromPath(params.kmer_DB.genome_fasta) : Channel.of()
genome_fai    = params.kmer_DB.genome_fai   ? Channel.fromPath(params.kmer_DB.genome_fai)   : Channel.of()
kiv2_bed      = params.kmer_DB.kiv2_bed     ? Channel.fromPath(params.kmer_DB.kiv2_bed)     : Channel.of()
norm_bed      = params.kmer_DB.norm_bed     ? Channel.fromPath(params.kmer_DB.norm_bed)     : Channel.of()

// Initialising count options:
rsids         = params.rsids_list           ? Channel.fromPath(params.rsids_list)           : Channel.of()
samplesheet   = params.count.samplesheet    ? Channel.fromPath(params.count.samplesheet)    : Channel.of()
kiv2_kmers    = params.count.kiv2_kmers     ? Channel.fromPath(params.count.kiv2_kmers)     : Channel.of()
norm_kmers    = params.count.norm_kmers     ? Channel.fromPath(params.count.norm_kmers)     : Channel.of()

kilda_py      = params.tools.kilda          ? Channel.fromPath(params.tools.kilda)          : Channel.of()


// Workflow to get the kmers unique to the KIV2 regions:
workflow prepare_kmers_kiv2 {
    take:
        genome_fasta
        genome_fai
        kiv2_bed

    main:
        kiv2_fasta = ExtractFastaFromBed(genome_fasta, genome_fai, kiv2_bed)
        kiv2_kmers = CountKmersRegion(kiv2_fasta)
        kiv2_kmers_filt = FilterOnOccurence(kiv2_kmers, 6)
        kiv2_outside_kmers = CountKmersOutsideRegion(genome_fasta, genome_fai, kiv2_bed)
        FilterKmersOccuringOutsideRegion(kiv2_kmers_filt, kiv2_outside_kmers)

    emit:
        FilterKmersOccuringOutsideRegion.out
}

// Workflow to get the kmers unique to the normalisation region(s) (recommended: the LPA gene (without the KIV2 region)):
workflow prepare_kmers_norm {
    take:
        genome_fasta
        genome_fai
        norm_bed

    main:
        norm_fasta = ExtractFastaFromBed(genome_fasta, genome_fai, norm_bed)
        norm_kmers = CountKmersRegion(norm_fasta)
        norm_kmers_filt = FilterOnOccurence(norm_kmers, 1)
        norm_outside_kmers = CountKmersOutsideRegion(genome_fasta, genome_fai, norm_bed)
        FilterKmersOccuringOutsideRegion(norm_kmers_filt, norm_outside_kmers)

    emit:
        FilterKmersOccuringOutsideRegion.out
}

// Workflow to build the kmer database:
workflow prepare_kmers_DB {
    take:
        genome_fasta
        genome_fai
        kiv2_bed
        norm_bed

    main:
        kiv2_kmers_list = prepare_kmers_kiv2(genome_fasta, genome_fai, kiv2_bed)
        norm_kmers_list = prepare_kmers_norm(genome_fasta, genome_fai, norm_bed)

        unique_kmers_lists = RemoveCommonKmers(kiv2_kmers_list, norm_kmers_list)
        OutputFasta(unique_kmers_lists)

    emit:
        RemoveCommonKmers.out
}


// Workflow to count the KIV2 repetitions:
workflow kiv2_counts {
    take:
        kiv2_kmers
        norm_kmers
        rsids
        samplesheet
        kilda_py

    main:
        norm_kiv2_fasta = CreateFastaKmers(kiv2_kmers, norm_kmers, rsids)

        // Read each sample as a tuple (sampleID, list of files) from the samplesheet.
        // Branch to accept a mix of FASTQs and BAMs:
        bam_or_fastq_ch = samplesheet.splitCsv(sep: '\t').
            branch { v ->
                        bams: v[1].endsWith("bam")
                            return [ v[0], file(v[1]), file(v[1]+".bai") ]
                        fastqs: true
                    }

        bams_ch = BamToFastq(kiv2_bed.first(),
                            norm_bed.first(),
                            genome_fasta.first(),
                            genome_fai.first(),
                            bam_or_fastq_ch.bams).fastq.map{ [ it[0], [file(it[1]), file(it[2])] ] }

        fastqs_ch = bam_or_fastq_ch.fastqs.map{ [ it[0], it[1].split(" ").collect(fastq -> file(fastq)) ] }

        input_ch = bams_ch.mix(fastqs_ch)

        jelly_kmers_ch = CountKmers(input_ch, norm_kiv2_fasta.first())
        counts_ch = DumpKmers(jelly_kmers_ch).collect()
        counts_list_ch = CreateSampleMap(counts_ch)
        
        kilda(kiv2_kmers, norm_kmers, kilda_py, counts_list_ch)
}


// Main workflow:
// Builds the kmer database if params.kmer_DB.build_DB is set to true
// Launch the KIV2 counts if params.count.count_kiv2 is set to true
workflow {

    // Check for input and builds the kmer DB:
    if(params.kmer_DB.build_DB) {
        if(!params.kmer_DB.outdir)    println("kmer_DB.outdir not set, defaulting to './kmers_DB'")
        if(!genome_fasta)    error("ERROR: A reference genome is needed to build the Kmer database (see config: 'genome_fasta')")
        if(!genome_fai)      error("ERROR: The reference genome needs to be indexed with samtools (see config: 'genome_fai')")
        if(!norm_bed)        error("ERROR: A bed delimiting the normalisation region(s) is needed to build the Kmer database  (see config: 'norm_bed')")
        if(!kiv2_bed)        error("ERROR: A bed delimiting the KIV2 region is needed to build the Kmer database  (see config: 'kiv2_bed')")

        // We fill the values with the built database:
        kiv2_kmers = "${outdir}/KIV2_kmers_6copies_specific.tsv"
        norm_kmers = "${outdir}/Norm_kmers_1copies_specific.tsv"

        prepare_kmers_DB(genome_fasta, genome_fai, kiv2_bed, norm_bed)
    }

    // Check for input and counts the KIV2 repeats:
    if(params.count.count_kiv2) {
        if(!kiv2_kmers)      error("ERROR: The directory to the list of KIV2 kmers must be provided to count the KIV2 (see config: 'kiv2_kmers')")
        if(!norm_kmers)      error("ERROR: The directory to the list of Normalisation kmers must be provided to count the KIV2 (see config: 'norm_kmers')")
        if(!samplesheet)     error("ERROR: A samplesheet must be provided to count the KIV2 (see config: 'samplesheet')")
        if(!kilda_py)        error("ERROR: The path to kilda.py must be provided (see tools.kilda)")

        kiv2_counts(kiv2_kmers, norm_kmers, rsids, samplesheet, kilda_py)
    }

}
