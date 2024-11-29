
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Importing processes from the modules:
import { ExtractFastaFromBed }                  from './modules/kiv2_create_kmers_DB.nf'
import { CountKmersRegion }                     from './modules/kiv2_create_kmers_DB.nf'
import { CountKmersOutsideRegion }              from './modules/kiv2_create_kmers_DB.nf'
import { FilterKmersOccuringOutsideRegion }     from './modules/kiv2_create_kmers_DB.nf'
import { RemoveCommonKmers }                    from './modules/kiv2_create_kmers_DB.nf'
import { OutputFasta }                          from './modules/kiv2_create_kmers_DB.nf'

import { CreateFastaKmers }     from './modules/kiv2_counts.nf'
import { CountKmers }           from './modules/kiv2_counts.nf'
import { DumpKmers }            from './modules/kiv2_counts.nf'
import { CreateSampleMap }      from './modules/kiv2_counts.nf'


// Workflow to get the kmers unique to the KIV2 regions:
workflow prepare_kmers_kiv2 {
    main:
        kiv2_fasta = ExtractFastaFromBed(file(params.input.genome_fasta), 
                                         file(params.input.genome_fai), 
                                         file(params.input.kiv2_bed))

        kiv2_kmers = CountKmersRegion(kiv2_fasta)
        kiv2_kmers_filt = FilterOnOccurence(kiv2_kmers, 6)
    
        kiv2_outside_kmers = CountKmersOutsideRegion(file(params.input.genome_fasta), 
                                                     file(params.input.genome_fai), 
                                                     file(params.input.kiv2_bed))
        FilterKmersOccuringOutsideRegion(kiv2_kmers_filt, kiv2_outside_kmers)

    emit:
        FilterKmersOccuringOutsideRegion.out
}

// Workflow to get the kmers unique to the normalisation region(s) (recommended: the LPA gene (without the KIV2 region)):
workflow prepare_kmers_norm {
    main:
        // norm_bed = Channel.fromPath("${params.input.norm_bed}")
        norm_fasta = ExtractFastaFromBed(file(params.input.genome_fasta), 
                                         file(params.input.genome_fai), 
                                         file(params.input.norm_bed))

        norm_kmers = CountKmersRegion(norm_fasta)
        norm_kmers_filt = FilterOnOccurence(norm_kmers, 1)
        
        // norm_bed = Channel.fromPath("${params.input.norm_bed}")
        norm_outside_kmers = CountKmersOutsideRegion(file(params.input.genome_fasta), 
                                                     file(params.input.genome_fai), 
                                                     file(params.input.norm_bed))
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
    rsids_ch = Channel.of(params.input.rsids_list).map{ f -> if(f) { return file(f) } else { return [] } }
    fasta_ch = Channel.of(params.input.genome_fasta).map{ f -> if(f) { return file(f) } else { return [] } }
    fai_ch = Channel.of(params.input.genome_fai).map{ f -> if(f) { return file(f) } else { return [] } }
    kiv2_bed_ch = Channel.of(params.input.kiv2_bed).map{ f -> if(f) { return file(f) } else { return [] } }
    norm_bed_ch = Channel.of(params.input.norm_bed).map{ f -> if(f) { return file(f) } else { return [] } }

    norm_kiv2_fasta = CreateFastaKmers(file(params.input.kiv2_kmers),
                                       file(params.input.norm_kmers),
                                       rsids_ch)

    // Read each sample as a tuple (sampleID, list of files) from the samplesheet.
    // Branch to accept a mix of FASTQs and BAMs:
    bam_or_fastq_ch = Channel.fromPath("${params.input.samplesheet}").splitCsv(sep: '\t').
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
    
    kilda(file(params.input.kiv2_kmers),
          file(params.input.norm_kmers),
          file(params.tools.kilda),
          counts_list_ch)
}


// Main workflow:
// Builds the kmer database if the lists of kmers are not specified in the config file
// Launch the KIV2 counts if a samplesheet is provided in the config file
workflow {
    build_DB = false

    // Checking mandatory options:
    if(params.kilda_dir == "")          error: "ERROR: The path to KILDA installation directory must be provided (see config: 'kilda_dir')"
    if(params.wdir == "")               error: "ERROR: The path to the working directory must be provided (see config: 'wdir')"
    if(params.input.kmer_size < 1)      error: "ERROR: The kmer size must be > 0 (see config: 'kmer_size')"

    // If theses files are missing we are building the kmer DB from scratch:
    if(params.input.kmer_DB_dir == "")  build_DB = true
    if(params.input.kiv2_kmers == "")   build_DB = true
    if(params.input.norm_kmers == "")   build_DB = true

    if(build_DB == true) {
        if(params.input.genome_fasta == "")     error: "ERROR: A reference genome is needed to build the Kmer database (see config: 'genome_fasta')"
        if(params.input.genome_fai == "")       error: "ERROR: The reference genome needs to be indexed with samtools (see config: 'genome_fai')"
        if(params.input.norm_bed == "")         error: "ERROR: A bed delimiting the normalisation region(s) is needed to build the Kmer database  (see config: 'norm_bed')"
        if(params.input.kiv2_bed == "")         error: "ERROR: A bed delimiting the KIV2 region is needed to build the Kmer database  (see config: 'kiv2_bed')"

        // We fill the values with the built database:
        params.input.kmer_DB_dir = "${kilda_dir}/data/kmer_DB/"
        params.input.kiv2_kmers  = "${kmer_DB_dir}/KIV2_kmers_6copies_specific.tsv"
        params.input.norm_kmers  = "${kmer_DB_dir}/Norm_kmers_1copies_specific.tsv"

        prepare_kmers_DB()
        // Add some QCs ?
    }

    // If a samplesheet is provided, we count the kmers for the samples:
    if(params.input.samplesheet) {
        if(params.input.kiv2_kmers == "")   error: "ERROR: The directory to the list of KIV2 kmers must be provided to count the KIV2 (see config: 'kiv2_kmers')"
        if(params.input.norm_kmers == "")   error: "ERROR: The directory to the list of Normalisation kmers must be provided to count the KIV2 (see config: 'norm_kmers')"
        
        kiv2_counts()
    }
}