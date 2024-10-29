#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process create_fasta_kmers {
    publishDir "${params.wdir}/", mode: 'copy'
    input:
       path(kiv2_kmers)
       path(norm_kmers)
       path(rsids_list)
    output:
        path(norm_kiv2_fasta)
    
    script:
        norm_kiv2_fasta = "NORM_KIV2_kmers.fasta"
    
        """
        set -eo pipefail
        
        awk '{ print ">KIV2_"NR"\\n"\$1 }' ${kiv2_kmers} > ${norm_kiv2_fasta};
        awk '{ print ">NORM_"NR"\\n"\$1 }' ${norm_kmers} >> ${norm_kiv2_fasta};
        # Adding the ref (2nd column) and alt (3rd column) kmers for each rsid (1st column) in "rsids_list"
        # with "rsid_ref" and "rsid_alt" as headers:
        awk '{ print ">"\$1"_ref\\n"\$2"\\n>"\$1"_alt\\n"\$3 }' ${rsids_list} >> ${norm_kiv2_fasta};
        """
}


// Counting all the kmers in the fastq files
process count_kmers {
    input:
        tuple val(sampleID), path(fastqs)
        val(norm_kiv2_fasta)
    
    output:
        path(jellyfish_kmers)
    
    script:
        jellyfish_kmers = "${sampleID}.kmers"

        """
        set -eo pipefail
        ${params.tools.jellyfish} count -t ${task.cpus} -m ${params.input.kmer_size} -s 100M -C -o ${jellyfish_kmers} --if=${norm_kiv2_fasta} <(zcat -f ${fastqs.join(" ")});
        """
}


// Querying the relevant kmers (KIV2 + normalisation regions)
process dump_kmers {
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
       path(kilda_py)
       path(counts_list)
    output:
        path(kiv2_outdir)
    
    script:
        kiv2_outdir = "kilda_kiv2_CNs/"
        
        rsid_param = ""
        if("${params.input.rsids_list}" != "") {
            rsid_param = " -r ${params.input.rsids_list}"
        }

        """
        set -eo pipefail

        ${params.tools.python} ${kilda_py} \
        -c ${counts_list} -o ${kiv2_outdir} -v -p \
        -k ${kiv2_kmers} -l ${norm_kmers} \
        ${rsid_param}
        """
}


process bam_to_fastq {
tag "${sample}"
afterScript "rm -rf TMP"
input:
  path(kiv2bed)
  path(normbed)
  path(fasta)
  path(fai)
  tuple val(sample),path(bam), path(bai)
output:
  tuple val(sample),path("${sample}.R1.fq.gz"),path("${sample}.R2.fq.gz"),emit:fastq
script:
"""
set -eo pipefail
mkdir -p TMP

cat "${kiv2bed}" "${normbed}" > TMP/select.bed

${params.tools.samtools} view -F '3844' --uncompressed -O BAM  -M -L TMP/select.bed --threads 1 --reference "${fasta}" "${bam}"  |\\
${params.tools.samtools} collate -f --threads ${Math.max(1,(task.cpus as int) -2)} -O -u --no-PG --reference "${fasta}" "-" TMP/tmp.collate |\\
${params.tools.samtools} fastq -n --threads 1 -1 TMP/${sample}.R1.fq.gz -2 TMP/${sample}.R2.fq.gz -s /dev/null -0 /dev/null
mv -v TMP/${sample}.R1.fq.gz ./
mv -v TMP/${sample}.R2.fq.gz ./
"""
}

workflow {
    norm_kiv2_fasta = create_fasta_kmers(file(params.input.kiv2_kmers), file(params.input.norm_kmers), file(params.input.rsids_list) )

    // Read each sample as a tuple (sampleID, list of fastqs) from the samplesheet:
    //fastqs_ch = Channel.fromPath("${params.input.samplesheet}").splitCsv(sep: '\t').map{ [it[0], it[1].split(" ")] }

    bam_or_fastq_ch = Channel.fromPath("${params.input.samplesheet}").splitCsv(sep: '\t').
     branch{v->
            bam: v[1].endsWith(".bam")
            fastq: true
            }

    bam2fq_ch = bam_to_fastq(
		file(params.input.kiv2_bed),
		file(params.input.norm_bed),
		file(params.input.genome_fasta),
		file(params.input.genome_fai),
		bam_or_fastq_ch.bam).fastq.
                     map{[it[0],[it[1],it[2]]]}
 
    fastqs_ch = bam_or_fastq_ch.fastq.map{ [it[0], it[1].split(" ").collect(fastq -> file(fastq)) ] }.mix(bam2fq_ch)
  
    jelly_kmers_ch = count_kmers(fastqs_ch, norm_kiv2_fasta)
    counts_ch = dump_kmers(jelly_kmers_ch).collect()
    counts_list_ch = CreateSampleMap(counts_ch)
    
    kilda(file(params.input.kiv2_kmers), file(params.input.norm_kmers) , file(params.tools.kilda), counts_list_ch)
}
