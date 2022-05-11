#!/usr/bin/env nextflow

process index {
    cache 'lenient'
    label 'index'
    // denote a list of modules with : seperator 
    module '/isg/shared/modulefiles/bwa/0.7.17'
    input:
    path ref from params.reference
    output:
    // error if not a glob pattern
    path 'index.*' into index_ch
    script:
    """
    bwa index -p index ${ref}
    """
    stub:
    """
    touch index.amb
    """
}

// set channel on pooled forward and reverse files
Channel
    .fromFilePairs( params.fq_files, size: 2 ) // size: number of inputs expected in grouping
    .ifEmpty { error "Cannot find any reads matching: ${params.fq_files}" }
    .set { read_pairs_ch }

process fastp {
    cache 'lenient'
    tag "fastp on ${pair_id}" // associates each process execution with a custom label
    module '/isg/shared/modulefiles/fastp/0.23.2' 
    publishDir "${params.outdir}/01-fastp_trimmed", mode: 'copy', pattern: "${pair_id}.html" // publish only html reports

    input:
    // takes shared string from paired reads [pair_id_1, [/path/to/file/,/path/to/file],...]
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    tuple val(pair_id), path("trim_${pair_id}_{1,2}.fq.gz") into fastp_ch, fastp_ch1 // place output into a channel for tunneling to next process

    script:
    """
    fastp -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} \
    -o trim_${reads[0]} -O trim_${reads[1]} \
    --correction -h ${pair_id}.html &> ${pair_id}.log
    """
    // replaces the actual process when the -stub-run or -stub command line option is invoked 
    stub:
    """
    touch trim_Golden{001..012}_{1,2}.fq.gz
    touch Golden{001..012}_{1,2}.{html,log}
    """
}

/*
read_pairs_ch
    .groupTuple(by:0)
    .map {prefix, reads -> tuple(prefix, reads.sort{it.name} ) }
    .view {"read pairs $it"}
Channel
    .from(1..12)
    .map {fq -> tuple("Goldenchr", path("/some/path/foo.${chr}.indels.vcf"), path("/other/path/foo.snvs.${chr}.vcf")}
    .view{ "numbers $it"}
*/
Channel
    .fromPath(params.barcodes).set {barcodes}

//fastp_ch1.view {"fastp Channel $it"}
// takes items from channel and puts them into a sorted list
//barcodes.toSortedList().view {" barcodes $it"}

// split reads into respective individual files
process demultiplex {
    cache 'lenient'
    tag "demultiplexing on pool: ${pair_id}"
    module '/isg/shared/modulefiles/stacks/2.53'
    //publishDir "${params.outdir}/01-process", mode: 'copy', patter: '*.log'

    input:
    tuple val(pair_id), path(trimmed_reads) from fastp_ch
    path barcode from barcodes.collect() // place all symlinks barcodes in directory for grabbing, outputs all items in the channel

    output:
    // Missing value declared as output parameter: ind_id
    tuple val('Golden????'), path('Golden????.{1,2}.fq.gz') into demult_output_ch
    //path "*.log" into demultLog_ch

    script:
    """
    pool="\$(echo ${pair_id} | sed 's/Golden//')"
    process_radtags -i gzfastq -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} \
    -b barcodes_"\${pool}".txt -o ./ \
    -c -q -r \
    --inline_null \
    --renz_1 sbfI --renz_2 mseI    
    """
    stub:
    """
    touch Golden{1..3}A{01..12}.{1,2}.fq.gz
    """
}

// provide file id's based on barcodes file 

// function to determine if paired end or not from https://github.com/nf-core/demultiplex/blob/dev/workflows/demultiplex.nf
def getFastqPairName(fqfile) {
    // variable sampleName with regex pattern to filter for with .find()
    def sampleName = (fqfile =~ /.*\/(.+).[12]\.fq\.gz/)
    // creates a matcher for the input String
    if (sampleName.find()) {
        //  prints the pattern matched following regex matcher
        return sampleName.group(1)
    }
    return fqfile
}

demult_pruned_ch = demult_output_ch
    //.collect ()
    .flatten ()
    // closure with predicate
    .filter { it =~/.*fq.gz/ }
    // creates a new tuple with matching value
    .map { fastq -> [getFastqPairName(fastq), fastq] }
    .groupTuple ()
    //.view {"filter demultiplex channel $it"}
    //.set { demult_pruned_ch }

process bwa {
    cache 'lenient'
    tag "aligning $pair_id"
    module '/isg/shared/modulefiles/bwa/0.7.17'
    module '/isg/shared/modulefiles/samtools/1.9'

    input:
    // add literal string that is assigned the forward or reverse read and conveys that in the string meaning
        // forward_trim_read=contains("F")
    tuple val(pair_id), path(t_d_reads) from demult_pruned_ch
    path index from index_ch
    
    output:
    path("${pair_id}.bam") into aln_reads_ch, aln_ch
    path ("${pair_id}.bam.bai") into aln_indx_ch
    path("${pair_id}.stat") into stat_write_report_ch

    """
    # change index to explicit literal string
    bwa mem -t ${task.cpus} -R "@RG\\tID:${pair_id}\\tSM:${pair_id}" index ${t_d_reads[0]} ${t_d_reads[1]} | \
        samtools view -@ ${task.cpus} -S -h -u - | \
        samtools sort -@ ${task.cpus} - > ${pair_id}.bam
    samtools index -@ ${task.cpus} ${pair_id}.bam
    samtools stat -@ ${task.cpus} ${pair_id}.bam > ${pair_id}.stat
    """
}

// writes individual reports, need to concatanate reports maybe through collect and bash for loop
// remove bam records with low alignment rates before calling variants
process bwa_stats {
    tag "bwa Stats: ${bam_stat}"
    label 'stats'
    publishDir "${params.outdir}/bwa", mode: 'copy', pattern: '*.tsv'
    module '/isg/shared/modulefiles/samtools/1.9'

    input:
    path bam_stat from stat_write_report_ch.collect()

    output:
    path 'mapRate.tsv' into maprate_ch
    
    // tried placing bash code into bin dir. caused only one individual to be written into mapRate.tsv, even with .collect()
    script:
    """
    bam_files="\$( echo ${bam_stat} )"
    for STAT in \${bam_files[@]}; do
        IND="\$(basename \${STAT} .stats)"
        SEQ="\$(cat \${STAT} | grep -P \$'SN\tsequences:' | awk '{print \$3}')"
        ALIG="\$(cat \${STAT} | grep -P \$'SN\treads mapped:' | awk '{print \$4}')"
        MRATE="\$(echo "scale=4; \${ALIG} / \${SEQ}" | bc)"
        echo -e "\${IND}\t\${MRATE}" >> mapRate.tsv
    done
    exit
    """
    stub:
    """
    touch bamstats.tsv
    """
}

id_aln_ch = aln_reads_ch // convert channel to [matcher, bam]
    .map { bam -> 
        def key = bam.name.toString().tokenize('.').get(0) // name removes path items
        return tuple(key, bam)}

cln_ali_ch = maprate_ch // Channel correction based on stats file from samtools stats 
    .splitCsv(sep:'\t') // converts file into channel formating (each row). else prints file literal string
    .map { // assign variable to column subsets and convert 2nd col. to class float
        def key = it[0].toString().tokenize('.').get(0) // similar to cut -d 
        def mappingrate = it[1].toFloat()
        [ key, mappingrate ]
        }
    .filter ({ key, mappingrate -> mappingrate >= .75}) // retain samples with a mapping greater than 75%
    .join( id_aln_ch ) // outputs [key, stat, bam]
    .map { it[2] } // retain only bam records (3rd column)
    //.into {clean_ali_ch; test_ch}

process index_reference {   
    input:
    path ref from params.reference
    output:
    path '*fai' into ref_index_ch
    script:
    """
    samtools faidx ${ref}
    """
}

intervals_ch = ref_index_ch
    .splitCsv(sep: '\t')
    .map { row ->
        // rows of interval lists
        if (row[0][0] != "@") {
            def interval_start = row[1].toLong()
            def interval_length = row[2].toLong()
            long start
            long end
            int width = 5000000

            if (!params.intervals) {
                // update interval start and length for .fai
                interval_start = 0
                interval_length = row[1].toLong()
            }

            while(interval_start < interval_length) {
                start = interval_start
                // add a slight overlap
                end = interval_start + width + 1000
                interval_start = end - 1000
                if (end > interval_length) {
                    end = interval_length
                    interval_start = end
                }
                // add the interval to the channel
                intervals_ch.bind( "${row[0]}:${start}-${end}" )
            }
        }
    }

process variant_calling {
    //cache 'lenient'
    
    module '/isg/shared/modulefiles/freebayes/1.3.4'
    tag "interval: ${interval}"
    cpus "${task.cpus}" // specify cpus for nextflow to run parallel with each
    input:
    // from https://github.com/brwnj/freebayes-nf/blob/master/main.nf
    path aln from cln_ali_ch.collect()
    path ref from params.reference
    path indx from aln_indx_ch.collect()
    path faidx from ref_index_ch
    each interval from intervals_ch
    
    output:
    path "fb_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz" into (splitvcf_ch, makelist_ch)
    path "fb_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz.csi" into vcfidx_ch

    script:
    """
    # Noah: write text file using bash
    # ls *.bam > bamlist.txt
    freebayes \
        --region ${interval} \
        --fasta-reference ${ref} \
        ${aln.collect { "--bam $it" }.join(" ")} \
        -m 30 -q 20 --min-coverage 100 --skip-coverage 50000 \
        | bgzip -c > fb_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz
    bcftools index fb_${interval.replaceAll(~/\:|\-|\*/, "_")}.vcf.gz
    """
    stub:
    """
    touch fb.vcf.gz
    """
}

// order vcf_ch

process make_vcf_list {
    label 'little_demon'
    input:
    path vcf from makelist_ch.collect()
    output:
    path 'vcfList.txt' into vcfList_ch 
    script:
    template 'makelist.py'
}

process merge_vcfs {
    publishDir "${params.outdir}/fb", mode:'copy', pattern: '*.vcf.gz'

    input:
    path vcf from splitvcf_ch.collect()
    path vcfidx from vcfidx_ch.collect()
    path fasta from params.reference
    path faidx from ref_index_ch
    path vcfList from vcfList_ch

    output:
    path "fb.vcf.gz" into fb_popmap_ch, fb_filtering_ch
    path "fb.vcf.gz.tbi"

    script:
    """
    # vcflib sort 50kb window
    gunzip -cd \$(cat ${vcfList}) | vcffirstheader | bgzip -c > fb_dirty.vcf.gz
    gsort fb_dirty.vcf.gz ${faidx} | vcfuniq \
        | bgzip -c > fb_dirty_sorted.vcf.gz
    bcftools norm -c all -f ${fasta} --multiallelics - --threads ${task.cpus} \
        --output fb.vcf.gz --output-type z \
        fb_dirty_sorted.vcf.gz
    tabix -p vcf fb.vcf.gz
    """
}

process variant_filtering {
    input:
    path vcf from fb_filtering_ch
    path meta from params.meta
    path probSamples from params.probSamples

    output:
    //path 'consensus.recode.vcf' into clean_vcf_ch
    tuple path("consensus.vcf"), path("final_indv.txt") into clean_vcf_ch, vcf2gwas_ch
    path ("final_indv.txt") into vcfIndv_ch

    script:
    template 'fb_F_SITE-PIPELINE.sh'

    stub:
    """
    touch consensus.vcf.gz
    """
}

// split meta into env files for groups w/ data a.k.a genocide
process meta2grps {
    label 'little_demon'

    input:
    path meta from params.meta
    path writegrps_fun from params.writegrps_function
    path vcfIndv from vcfIndv_ch
    output:
    path "*.grp" into grp2baypass_ch, grp2meta2env_ch

    script:
    """
    #!/usr/bin/env Rscript
    source("${writegrps_fun}")
    writegrps("${meta}","${vcfIndv}")
    """
}

process vcf2gwas_ch {
    cache 'lenient'

    input:
    tuple path(vcf), path(indv) from vcf2gwas_ch
    each grp from grp2baypass_ch.flatten()
    //tuple path(vcf), path(indvs) from vcf2baypass_ch
    
    output:
    path ("*.vcf") into vcf2rda_ch
    path ("*.ped") into ped2lea_ch, ped2prune_ch
    path ("*.frq") into frq2bay_ch
    tuple val (grp), path ("*.indv"), path ("*.pop_order") into reorder_bay2env_ch
    path ("*.loci") into baypass_loci_ch
    
    script:
    """
    env="\$(basename ${grp} .grp)"
    awk -F"\t" '{print \$1}' ${grp} > "\${env}_indv.list"
    bgzip ${vcf} 
    bcftools view -I --force-samples -S "\${env}_indv.list" ${vcf}.gz > "\${env}_${vcf}"
    bcftools query -l "\${env}_${vcf}" > "\${env}.indv"
    plink --vcf "\${env}_${vcf}" --recode --allow-extra-chr --out "\${env}"
    cat "\${env}_${vcf}" | perl ${baseDir}/bin/vcf2baypass.pl ${grp} "\${env}.frq"
    """
}

process meta2env {
    cache 'lenient'

    input:
    path (meta2env) from params.meta2env_function
    path (meta) from params.meta
    tuple val (grps), path (indv) , path (grp_order) from reorder_bay2env_ch
    path (grp) from grp2meta2env_ch.collect()

    output:
    //tuple path("*.lfmm"), path("*.grp")
    path ("*.env") into env2lea_ch
    path ("*.bayenv") into env2bay_ch

    script:
    """
    #!/usr/bin/env Rscript
    source("${meta2env}")
    meta2env("${meta}", "${indv}", "${grps}", "${grp_order}")
    """
}


lea_ch = ped2lea_ch
    .mix ( env2lea_ch )
    .map { env -> 
        def key = env.name.toString().tokenize('.').get(0)
        return tuple(key, env)}
    .groupTuple(size:2)

process lea {
    label 'gwas'
    cache 'lenient'
    cpus "${task.cpus}"
    publishDir "${params.outdir}/gwas/output_files", mode: 'copy'

    input:
    // [env, [vcf,env],]
    tuple val(env), path(files) from lea_ch
    path fst_func from params.fst_function
    path lea_func from params.lea_function

    output:
    tuple val (env), path ("*_zscores.txt") into lea_zscores_ch
    stdout result

    script:
    """
    #!/usr/bin/env Rscript
    library("LEA", lib="${params.rlibrariesPath}")
    library("dplyr", lib="${params.rlibrariesPath}")
    source("${lea_func}")
    source('${fst_func}')
    lea("${env}","${files[0]}","${files[1]}")
    """
}

baypass_ch = frq2bay_ch
    .mix(env2bay_ch)
    .map { env -> 
        def key = env.name.toString().tokenize('.').get(0)
        return tuple(key, env)}
    .groupTuple(size:2)

channel.from(1..3).into {baypass_iter_ch; maxiter_ch}

maxiter_ch.max(). into {rep_ch; repgroup_ch}

process baypass {

    input:
    tuple val (env), path (files) from baypass_ch
    each rep from baypass_iter_ch

    output:
    tuple val (env), path ("baypass_core_${env}_${rep}_summary_{betai_reg,pi_xtx}.out") into baypass_rawoutput_ch, bay2spliter_ch
    
    script:
    """
    EPOP="\$(head -n 1 "${files[1]}" | awk '{print NF}')"
    POP="\$(head -n 1 "${files[0]}" | awk '{print NF}')"
    if ["\$EPOP" != "\$POP"]; then 
        echo "different number of populations in environment and genetic dataframes"
    else
        baypass -npop \${EPOP} -gfile "${files[0]}" -efile "${files[1]}" \
	    -nval 100 -thin 25 -burnin 50 -npilot 30 -pilotlength 10 \
	    -outprefix baypass_core_${env}_${rep} \
	    -nthreads ${task.cpus} -seed "\${RANDOM}"
    fi
    """
}

// add variable that is max value from rep_ch * 2
bay2merge_ch = bay2spliter_ch.transpose().groupTuple(size:6)

process merge_baypass {

    input:
    tuple val (env), path (files) from bay2merge_ch
    path (loci) from baypass_loci_ch
    path (merge_bay) from params.merge_baypass_function
    val (maxiter) from rep_ch

    output:
    tuple val (env), path ("${env}_med_baypass.txt") into baypass_med_ch

    script:
    """
    #!/usr/bin/env Rscript
    library("dplyr", lib="${params.rlibrariesPath}")
    library("matrixStats", lib="${params.rlibrariesPath}")
    source("${merge_bay}")
    merge_baypass("${env}","${maxiter}", "${loci}")
    """
}

baypass_med_ch
    .join(lea_zscores_ch)
    .view()


