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

include { CreateFastaKmers }                    from './modules/kiv2_counts.nf'
include { CountKmers }                          from './modules/kiv2_counts.nf'
include { DumpKmers }                           from './modules/kiv2_counts.nf'
include { CreateSampleMap }                     from './modules/kiv2_counts.nf'
include { kilda }                               from './modules/kiv2_counts.nf'

// Initialising the options with default values:
kilda_dir              = "${projectDir}" // projectDir: directory where the main script is located
// Being replaced by the .conf file:
params.wdir            = "${launchDir}"  // launchDir:  directory where the main script was launched 
params.kmer_size       = 31

params.outdir          = "./kmers_DB"
params.genome_fasta    = ""
params.genome_fai      = ""
params.kiv2_bed        = ""
params.norm_bed        = ""

params.samplesheet     = ""
params.kiv2_kmers      = ""
params.norm_kmers      = ""
params.rsids_list      = ""

// Checking the params values:
if(params.kmer_size < 1)                                            error("\nERROR: The kmer size must be > 0 (see config: 'kmer_size')")
// If we make a new kmer DB, then the options to indicate the list of kmers must not be set
// (since we are going to use the one produced from the "kmer_DB" step!)
if(params.build_DB && 
    params.count_kiv2 && 
        (params.kiv2_kmers != "" || params.norm_kmers != ""))       error("\nERROR: conflict in the config file: if 'build_DB' is true, 'kiv2_kmers' and 'norm_kmers' must be empty") 

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

        RemoveCommonKmers(kiv2_kmers_list, norm_kmers_list)
        //OutputFasta(unique_kmers_lists)

    emit:
        kiv2_list = RemoveCommonKmers.out.kiv2_list
        norm_list = RemoveCommonKmers.out.norm_list
}


// Workflow to count the KIV2 repetitions:
workflow kiv2_counts {
    take:
        kiv2_kmers
        norm_kmers
        rsids
        samplesheet

    main:
        // Read each sample as a tuple (sampleID, list of files) from the samplesheet.
        fastqs_ch = samplesheet.splitCsv(sep: '\t').map{ [ it[0], it[1].split(" ").collect(fastq -> file(fastq)) ] }

        norm_kiv2_fasta = CreateFastaKmers(kiv2_kmers, norm_kmers, rsids)
        jelly_kmers_ch = CountKmers(fastqs_ch, norm_kiv2_fasta.first())
        counts_ch = DumpKmers(jelly_kmers_ch).collect()
        counts_list_ch = CreateSampleMap(counts_ch)
        
        kilda(kiv2_kmers, norm_kmers, counts_list_ch)
}


// Main workflow:
// Builds the kmer database if params.build_DB is set to true
// Launch the KIV2 counts if params.count_kiv2 is set to true
workflow {
    
    // Check for input and builds the kmer DB:
    if(params.build_DB) {
        if(params.genome_fasta == "")    error("\nERROR: A reference genome is needed to build the Kmer database (see config: 'genome_fasta')")
        if(params.genome_fai   == "")    error("\nERROR: The reference genome needs to be indexed with samtools (see config: 'genome_fai')")
        if(params.norm_bed     == "")    error("\nERROR: A bed delimiting the normalisation region(s) is needed to build the Kmer database  (see config: 'norm_bed')")
        if(params.kiv2_bed     == "")    error("\nERROR: A bed delimiting the KIV2 region is needed to build the Kmer database  (see config: 'kiv2_bed')")

        println("Working directory: '"+params.wdir+"'")

        prepare_kmers_DB(Channel.fromPath(params.genome_fasta), 
                         Channel.fromPath(params.genome_fai), 
                         Channel.fromPath(params.kiv2_bed), 
                         Channel.fromPath(params.norm_bed))

        kiv2_kmers = prepare_kmers_DB.out.kiv2_list
        norm_kmers = prepare_kmers_DB.out.norm_list
    }

    // Check for input and counts the KIV2 repeats:
    if(params.count_kiv2) {
        if(params.samplesheet == "")    error("\nERROR: A samplesheet must be provided to count the KIV2 (see config: 'samplesheet')")
        if(params.rsids_list  == "")    prinln("rsids file not provided, KILDA will not check for variants in the kmers.")

        if(params.build_DB == false) {
            if(params.norm_kmers == "" || params.kiv2_kmers == "") { 
                error("\nERROR: The path to the lists of kiv2 and normalisation kmers must be provided to count the KIV2 (see config: 'norm_kmers' and 'kiv2_kmers')")
            } else {
                norm_kmers = Channel.fromPath(params.norm_kmers)
                kiv2_kmers = Channel.fromPath(params.kiv2_kmers)
            }
        }

        kiv2_counts(kiv2_kmers,
                    norm_kmers,
                    Channel.fromPath(params.rsids_list),
                    Channel.fromPath(params.samplesheet))
    }
}
