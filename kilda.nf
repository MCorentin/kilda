#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Importing processes from the modules:
include { ExtractFastaFromBed }                  from './modules/kiv2_create_kmers_DB.nf'
include { CountKmersRegion }                     from './modules/kiv2_create_kmers_DB.nf'
include { FilterOnOccurence }                    from './modules/kiv2_create_kmers_DB.nf'
include { CountKmersOutsideRegion }              from './modules/kiv2_create_kmers_DB.nf'
include { FilterKmersOccuringOutsideRegion }     from './modules/kiv2_create_kmers_DB.nf'
include { RemoveCommonKmers }                    from './modules/kiv2_create_kmers_DB.nf'
include { OutputFasta }                          from './modules/kiv2_create_kmers_DB.nf'

include { CreateFastaKmers }     from './modules/kiv2_counts.nf'
include { CountKmers }           from './modules/kiv2_counts.nf'
include { DumpKmers }            from './modules/kiv2_counts.nf'
include { CreateSampleMap }      from './modules/kiv2_counts.nf'


// Workflow to get the kmers unique to the KIV2 regions:
workflow prepare_kmers_kiv2 {
    main:
        kiv2_fasta = ExtractFastaFromBed(file(params.kmer_DB.genome_fasta), 
                                         file(params.kmer_DB.genome_fai), 
                                         file(params.kmer_DB.kiv2_bed))

        kiv2_kmers = CountKmersRegion(kiv2_fasta)
        kiv2_kmers_filt = FilterOnOccurence(kiv2_kmers, 6)
    
        kiv2_outside_kmers = CountKmersOutsideRegion(file(params.kmer_DB.genome_fasta), 
                                                     file(params.kmer_DB.genome_fai), 
                                                     file(params.kmer_DB.kiv2_bed))
        FilterKmersOccuringOutsideRegion(kiv2_kmers_filt, kiv2_outside_kmers)

    emit:
        FilterKmersOccuringOutsideRegion.out
}

// Workflow to get the kmers unique to the normalisation region(s) (recommended: the LPA gene (without the KIV2 region)):
workflow prepare_kmers_norm {
    main:
        norm_fasta = ExtractFastaFromBed(file(params.kmer_DB.genome_fasta), 
                                         file(params.kmer_DB.genome_fai), 
                                         file(params.kmer_DB.norm_bed))

        norm_kmers = CountKmersRegion(norm_fasta)
        norm_kmers_filt = FilterOnOccurence(norm_kmers, 1)
        
        norm_outside_kmers = CountKmersOutsideRegion(file(params.kmer_DB.genome_fasta), 
                                                     file(params.kmer_DB.genome_fai), 
                                                     file(params.kmer_DB.norm_bed))
        FilterKmersOccuringOutsideRegion(norm_kmers_filt, norm_outside_kmers)
        
    emit:
        FilterKmersOccuringOutsideRegion.out
}

// Workflow to build the kmer database:
workflow prepare_kmers_DB {
   kiv2_kmers_list = prepare_kmers_kiv2()
   norm_kmers_list = prepare_kmers_norm()

   unique_kmers_lists = RemoveCommonKmers(kiv2_kmers_list, norm_kmers_list)
   OutputFasta(unique_kmers_lists)

   emit:
        RemoveCommonKmers.out
}


// Workflow to count the KIV2 repetitions:
workflow kiv2_counts {
    // Optional inputs: we use an empty list "[]" when the file is not provided ([] is a valid path in nextflow):
    rsids_ch = Channel.of(params.rsids_list).map{ f ->        if(f) { return file(f) } else { return [] } }
    fasta_ch = Channel.of(params.kmer_DB.genome_fasta).map{ f ->    if(f) { return file(f) } else { return [] } }
    fai_ch = Channel.of(params.kmer_DB.genome_fai).map{ f ->        if(f) { return file(f) } else { return [] } }
    kiv2_bed_ch = Channel.of(params.kmer_DB.kiv2_bed).map{ f ->     if(f) { return file(f) } else { return [] } }
    norm_bed_ch = Channel.of(params.kmer_DB.norm_bed).map{ f ->     if(f) { return file(f) } else { return [] } }

    norm_kiv2_fasta = CreateFastaKmers(file(params.count.kiv2_kmers),
                                       file(params.count.norm_kmers),
                                       rsids_ch)

    // Read each sample as a tuple (sampleID, list of files) from the samplesheet.
    // Branch to accept a mix of FASTQs and BAMs:
    bam_or_fastq_ch = Channel.fromPath("${params.count.samplesheet}").splitCsv(sep: '\t').
     branch { v ->
                bams: v[1].endsWith("bam")
                    return [ v[0], file(v[1]), file(v[1]+".bai") ]
                fastqs: true
            }

    bams_ch = bam_to_fastq(kiv2_bed_ch.first(),
                           norm_bed_ch.first(),
                           fasta_ch.first(),
                           fai_ch.first(),
                           bam_or_fastq_ch.bams).fastq.map{ [ it[0], [file(it[1]), file(it[2])] ] }
 
    fastqs_ch = bam_or_fastq_ch.fastqs.map{ [ it[0], it[1].split(" ").collect(fastq -> file(fastq)) ] }
  
    input_ch = bams_ch.mix(fastqs_ch)

    jelly_kmers_ch = CountKmers(input_ch, norm_kiv2_fasta.first())
    counts_ch = DumpKmers(jelly_kmers_ch).collect()
    counts_list_ch = CreateSampleMap(counts_ch)
    
    kilda(file(params.count.kiv2_kmers),
          file(params.count.norm_kmers),
          file(params.tools.kilda),
          counts_list_ch)
}


// Main workflow:
// Builds the kmer database if the lists of kmers are not specified in the config file
// Launch the KIV2 counts if a samplesheet is provided in the config file
workflow {
    // Checking mandatory options:
    if(params.kilda_dir == "")          error: "ERROR: The path to KILDA installation directory must be provided (see config: 'kilda_dir')"
    if(params.wdir == "")               error: "ERROR: The path to the working directory must be provided (see config: 'wdir')"
    if(params.kmer_size < 1)      error: "ERROR: The kmer size must be > 0 (see config: 'kmer_size')"

    if(params.kmer_DB.build_DB) {
        if(params.kmer_DB.outdir == "")   error: "ERROR: You need to provide a kmer DB directory (see config: 'outdir')"
        if(params.kmer_DB.genome_fasta == "")     error: "ERROR: A reference genome is needed to build the Kmer database (see config: 'genome_fasta')"
        if(params.kmer_DB.genome_fai == "")       error: "ERROR: The reference genome needs to be indexed with samtools (see config: 'genome_fai')"
        if(params.kmer_DB.norm_bed == "")         error: "ERROR: A bed delimiting the normalisation region(s) is needed to build the Kmer database  (see config: 'norm_bed')"
        if(params.kmer_DB.kiv2_bed == "")         error: "ERROR: A bed delimiting the KIV2 region is needed to build the Kmer database  (see config: 'kiv2_bed')"

        // We fill the values with the built database:
        params.count.kiv2_kmers     = "${params.kmer_DB.outdir}/KIV2_kmers_6copies_specific.tsv"
        params.count.norm_kmers     = "${params.kmer_DB.outdir}/Norm_kmers_1copies_specific.tsv"

        prepare_kmers_DB()
        // Add some QCs ?
    }

    // If a samplesheet is provided, we count the kmers for the samples:
    if(params.count.count_kiv2) {
        if(params.count.kiv2_kmers == "")   error: "ERROR: The directory to the list of KIV2 kmers must be provided to count the KIV2 (see config: 'kiv2_kmers')"
        if(params.count.norm_kmers == "")   error: "ERROR: The directory to the list of Normalisation kmers must be provided to count the KIV2 (see config: 'norm_kmers')"
        if(params.count.samplesheet == "")  error: "ERROR: A samplesheet must be provided to count the KIV2 (see config: 'samplesheet')"

        kiv2_counts()
    }
}