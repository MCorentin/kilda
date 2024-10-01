#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process create_fasta_kmers {
    publishDir "${params.wdir}/", mode: 'copy'
    
    output:
        path(norm_kiv2_fasta)
    
    script:
        norm_kiv2_fasta = "NORM_KIV2_kmers.fasta"
    
        """
        set -eo pipefail
        
        awk '{ print ">KIV2_"NR"\\n"\$1 }' ${params.input.kiv2_kmers} > ${norm_kiv2_fasta};
        awk '{ print ">NORM_"NR"\\n"\$1 }' ${params.input.norm_kmers} >> ${norm_kiv2_fasta};
        # Adding the ref (2nd column) and alt (3rd column) kmers for each rsid (1st column) in "rsids_list"
        # with "rsid_ref" and "rsid_alt" as headers:
        awk '{ print ">"\$1"_ref\\n"\$2"\\n>"\$1"_alt\\n"\$3 }' ${params.input.rsids_list} >> ${norm_kiv2_fasta};
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

        ${params.tools.python} ${params.tools.kilda} \
        -c ${counts_list} -o ${kiv2_outdir} -v -p \
        -k ${params.input.kiv2_kmers} -l ${params.input.norm_kmers} \
        ${rsid_param}
        """
}


workflow {
    norm_kiv2_fasta = create_fasta_kmers()

    // Read each sample as a tuple (sampleID, list of fastqs) from the samplesheet:
    //fastqs_ch = Channel.fromPath("${params.input.samplesheet}").splitCsv(sep: '\t').map{ [it[0], it[1].split(" ")] }

    fastqs_ch = Channel.fromPath("${params.input.samplesheet}").splitCsv(sep: '\t').map{ [it[0], it[1].split(" ").collect(fastq -> file(fastq)) ] }
  
    jelly_kmers_ch = count_kmers(fastqs_ch, norm_kiv2_fasta)
    counts_ch = dump_kmers(jelly_kmers_ch).collect()
    counts_list_ch = CreateSampleMap(counts_ch)
    
    kilda(counts_list_ch)
}
