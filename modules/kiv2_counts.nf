#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CreateFastaKmers {
    publishDir "${params.wdir}/", mode: 'copy'

    input:
       path(kiv2_kmers)
       path(norm_kmers)
       path(rsids_list)
    
    output:
        path(norm_kiv2_fasta)
    
    script:
        norm_kiv2_fasta = "NORM_KIV2_kmers.fasta"

        // rsids is optional, if the file was not provided in the config file, do not add the rsid kmers to the FASTA:
        add_rsids = ""
        if(rsids_list) add_rsids = "awk '{ print \">\"\$1\"_ref\\n\"\$2\"\\n>\"\$1\"_alt\\n\"\$3 }' ${rsids_list} >> ${norm_kiv2_fasta};"

        """
        set -eo pipefail

        awk '{ print ">KIV2_"NR"\\n"\$1 }' ${kiv2_kmers} > ${norm_kiv2_fasta};
        awk '{ print ">NORM_"NR"\\n"\$1 }' ${norm_kmers} >> ${norm_kiv2_fasta};

        ${add_rsids}
        """
}


// Counting all the kmers in the fastq files
process CountKmers {
    input:
        tuple val(sampleID), path(fastqs)
        val(norm_kiv2_fasta)
    
    output:
        path(jellyfish_kmers)
    
    script:
        jellyfish_kmers = "${sampleID}.kmers"

        """
        set -eo pipefail
        ${params.tools.jellyfish} count -t ${task.cpus} -m ${params.kmer_size} -s 100M -C -o ${jellyfish_kmers} --if=${norm_kiv2_fasta} <(zcat -f ${fastqs.join(" ")});
        """
}


// Querying the relevant kmers (KIV2 + normalisation regions)
process DumpKmers {
    publishDir "${params.wdir}/counts/", mode: 'copy'

    input:
        path(jellyfish_kmers)
        
    output:
        path(counts)
    
    script:
        counts = "${jellyfish_kmers}".replaceAll(/.kmers/, ".counts")
        
        """
        set -eo pipefail
        ${params.tools.jellyfish} dump ${jellyfish_kmers} -c -t > ${counts}
        """
}
        

// Generate the list of counts, to be used as input by the "kilda" process
process CreateSampleMap {
    publishDir "${params.wdir}/", mode: 'copy'

    input:
        path(counts)

    output:
        path(counts_list), emit: counts_list

    script:	
        counts_list = "counts.list"
        list_content = ""
        
        for(String count : counts) {
            // "split("/").last()" as a basename, then removing the ".counts" to get the sampleID (cf previous processes for naming conventions)
            sampleID = count.split("/").last().split(".counts").first()    
            list_content = list_content + sampleID + "\t${params.wdir}/counts/" + count + "\n"
        }
        // Remove the last trailing '\n' (empty final line)
        list_content = list_content.trim()

        """
        set -eo pipefail
        echo "${list_content}" > ${counts_list}
        """
}

// Counting the number of kiv2, from all the counts files
process kilda {
    publishDir "${params.wdir}/", mode: 'copy'
    
    input:
       path(kiv2_kmers)
       path(norm_kmers)
       path(counts_list)
    
    output:
        path(kiv2_outdir)
    
    script:
        kiv2_outdir = "kilda_kiv2_CNs/"
        
        rsid_param = ""
        if("${params.rsids_list}" != "") {
            rsid_param = " -r ${params.rsids_list}"
        }

        """
        set -eo pipefail

        kilda.py \
        -c ${counts_list} -o ${kiv2_outdir} -v -p \
        -k ${kiv2_kmers} -l ${norm_kmers} \
        ${rsid_param}
        """
}


process BamToFastq {
    tag "${sample}"
    afterScript "rm -rf TMP"

    input:
        path(kiv2bed)
        path(normbed)
        path(fasta)
        path(fai)
        tuple val(sample), path(bam), path(bai)
    
    output:
        tuple val(sample), path("${sample}.R1.fq.gz"), path("${sample}.R2.fq.gz"), emit: fastq
    
    script:
        // Here, a reference FASTA is needed to convert BAM/CRAM to FASTQ:
        if(!fasta || !fai) error("A reference FASTA (+fai) should be provided in the configuration file if BAM/CRAM files are present in the samplesheet.")

        cmd_cat_bed = ""

        // If beds not provided: we extract all the reads:
        cmd_collate = "${params.tools.samtools} collate -f --threads ${Math.max(1,(task.cpus as int) -2)} -O -u --no-PG --reference ${fasta} ${bam} TMP/tmp.collate"
        cmd_fastq = "${params.tools.samtools} fastq -n --threads 1 -1 TMP/${sample}.R1.fq.gz -2 TMP/${sample}.R2.fq.gz -s /dev/null -0 /dev/null"
        cmd_b2f = "${cmd_collate} | ${cmd_fastq}"

        // If beds provided: we extract only the corresponding regions:
        if(kiv2bed && normbed) {
            cmd_cat_bed = "cat ${kiv2bed} ${normbed} > TMP/select.bed"

            cmd_view = "${params.tools.samtools} view -F '3844' --uncompressed -O BAM  -M -L TMP/select.bed --threads 1 --reference ${fasta} ${bam}"
            cmd_collate = "${params.tools.samtools} collate -f --threads ${Math.max(1,(task.cpus as int) -2)} -O -u --no-PG --reference ${fasta} \"-\" TMP/tmp.collate"
            cmd_b2f = "${cmd_view} | ${cmd_collate} | ${cmd_fastq}"
        } 

        """
        set -eo pipefail

        mkdir -p TMP

        ${cmd_cat_bed}
        
        ${cmd_b2f}

        mv -v TMP/${sample}.R1.fq.gz ./
        mv -v TMP/${sample}.R2.fq.gz ./
        """
}

