#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process extract_fasta_from_bed {
    input:
        path(fasta)
        path(fai)
        path(region_bed)
        
    output:
        path(region_fasta)
        
    script:
        region_fasta = "${region_bed.baseName}.fasta"
        
        """
        set -eo pipefail
        
        awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' ${region_bed} > region.regions
        ${params.tools.samtools} faidx ${fasta} -r region.regions > ${region_fasta}
        """
}

process count_kmers_region {
    input:
        path(region_fasta)
        
    output:
        path(kmers)
    
    script:
        kmers = "${region_fasta.baseName}_kmers"
        
        """
        set -eo pipefail
        ${params.tools.jellyfish} count -C -t ${task.cpus} -m ${params.input.kmer_size} -s 10G -o ${kmers} ${region_fasta}
        """
}


process filter_on_occurence {
    input:
        path(kmers)
        val(occurence)
        
    output:
        path(kmers_filtered_dump)
    
    script:
        kmers_filtered_dump = "${kmers}_${occurence}copies.dump"
        
        """
        set -eo pipefail
        ${params.tools.jellyfish} dump -c -t -L ${occurence} -U ${occurence} -o ${kmers_filtered_dump} ${kmers}
        """
}


process count_kmers_outside_region {
    input:
        path(fasta)
        path(fai)
        path(region_bed)
        
    output:
        path(kmers_outside_region_count)
    
    script:
        kmers_outside_region_count = "genome_without_region.kmers"
        
        """
        set -eo pipefail
        
        awk -F '\t' '{ print \$1"\\t0\\t"\$2 }' '${fai}' > genome.bed
        
        ${params.tools.bedtools} intersect -v -a genome.bed -b ${region_bed} | awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' > genome_without_region.regions
        ${params.tools.samtools} faidx ${fasta} -r genome_without_region.regions > genome_without_region.fasta
        
        ${params.tools.jellyfish} count -C -t ${task.cpus} -m ${params.input.kmer_size} -s 10G -o ${kmers_outside_region_count} genome_without_region.fasta
        """
}


process filter_kmers_occuring_outside_region {
    input:
        path(kmers_filtered_dump)
        path(kmers_outside_region_count)
        
    output:
        path(kmers_specific)
    
    script:
        kmers_filtered_fasta = "${kmers_filtered_dump.baseName}.fasta"
        kmers_specific = "${kmers_filtered_dump.baseName}_specific.fasta"
        
        """
        set -eo pipefail
        
        cut -f 1 ${kmers_filtered_dump} | awk '{print ">kmer_"NR"\\n"\$1 }' > ${kmers_filtered_fasta}
        ${params.tools.jellyfish} query -s ${kmers_filtered_fasta} -o kmers_outside_occurences.counts ${kmers_outside_region_count}
        
        awk '\$2 == 0 { print \$1 }' kmers_outside_occurences.counts > specific_kmers.list
        cut -f 1 specific_kmers.list > ${kmers_specific}
        """
}


process remove_common_kmers {
    publishDir "${params.kmer_DB_outdir}/", mode: 'copy'
    
    input:
        path(kiv2_kmers_list)
        path(norm_kmers_list)
        
    output:
        tuple(path(kiv2_unique_kmers_list), path(norm_unique_kmers_list))
        
    script:
        kiv2_unique_kmers_list = "${kiv2_kmers_list.baseName}.tsv"
        norm_unique_kmers_list = "${norm_kmers_list.baseName}.tsv"
    
        """
        set -eo pipefail
        
        grep -vf ${kiv2_kmers_list} ${norm_kmers_list} > ${norm_unique_kmers_list}
        grep -vf ${norm_kmers_list} ${kiv2_kmers_list} > ${kiv2_unique_kmers_list}
        """
}


process output_fasta {
    publishDir "${params.kmer_DB_outdir}/", mode: 'copy'
    
    input:
        tuple(path(kiv2_unique_kmers_list), path(norm_unique_kmers_list))
    
    output:
        tuple(path(kiv2_unique_kmers_fasta), path(norm_unique_kmers_fasta))

        
    script:
        kiv2_unique_kmers_fasta = "${kiv2_unique_kmers_list.baseName}.fasta"
        norm_unique_kmers_fasta = "${norm_unique_kmers_list.baseName}.fasta"
        
        """
        set -eo pipefail
        
        awk '{print ">kmer_KIV2_"NR"\\n"\$1 }' ${kiv2_unique_kmers_list} > ${kiv2_unique_kmers_fasta}
        awk '{print ">kmer_NORM_"NR"\\n"\$1 }' ${norm_unique_kmers_list} > ${norm_unique_kmers_fasta}
        """
}


workflow prepare_kmers_kiv2 {
    main:
        KIV2_bed = Channel.fromPath("${params.input.kiv2_bed}")
        
        KIV2_fasta = extract_fasta_from_bed(file(params.input.genome_fasta),file(params.input.genome_fai),KIV2_bed)
        kiv2_kmers = count_kmers_region(KIV2_fasta)
        
        kiv2_kmers_filt = filter_on_occurence(kiv2_kmers, 6)
        
        kiv2_bed = Channel.fromPath("${params.input.kiv2_bed}")
        kiv2_outside_kmers = count_kmers_outside_region(file(params.input.genome_fasta), file(params.input.genome_fai), kiv2_bed)
        
        filter_kmers_occuring_outside_region(kiv2_kmers_filt, kiv2_outside_kmers)

    emit:
        filter_kmers_occuring_outside_region.out
}

// Normalisation region (recommended: the exons of LDLR, APOB and PCSK9): 
workflow prepare_kmers_norm {
    main:
        norm_bed = Channel.fromPath("${params.input.norm_bed}")
        
            norm_fasta = extract_fasta_from_bed(file(params.input.genome_fasta),file(params.input.genome_fai),norm_bed)
            norm_kmers = count_kmers_region(norm_fasta)
        
        norm_kmers_filt = filter_on_occurence(norm_kmers, 1)
        
        norm_bed = Channel.fromPath("${params.input.norm_bed}")
        norm_outside_kmers = count_kmers_outside_region(file(params.input.genome_fasta),file(params.input.genome_fai),norm_bed)
        
        filter_kmers_occuring_outside_region(norm_kmers_filt, norm_outside_kmers)
        
    emit:
        filter_kmers_occuring_outside_region.out
}


workflow {
   kiv2_kmers_list = prepare_kmers_kiv2()
   norm_kmers_list = prepare_kmers_norm()

   unique_kmers_lists = remove_common_kmers(kiv2_kmers_list, norm_kmers_list)
   output_fasta(unique_kmers_lists)
}
