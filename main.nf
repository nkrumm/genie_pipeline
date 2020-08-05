

// parameters; default is for test data
def sample               = params.sample_id ?: "20-90206-1"
def assay                = params.assay ?: 'test_assay'
def ref_build            = params.genome ?: 'hg38' 

// base paths -- typically defined in config
def ref_base             = params.ref_base
def assay_base           = params.assay_base
def sample_base          = params.sample_base

assert params.assays[assay].compatible_references.contains(ref_build)


// assay-specific files
assay_intervals      = Channel.fromPath("${assay_base}/${params.assays[assay].intervals}").collect()

// reference sequence
ref_fasta            = file("${ref_base}/${params.references[ref_build].fasta}", checkIfExists: true)
ref_fasta_fai        = file("${ref_base}/${params.references[ref_build].fastaFai}", checkIfExists: true)

// BWA index for alignment
bwa_index            = Channel.fromPath("${ref_base}/${params.references[ref_build].bwaIndex}", checkIfExists: true).collect()

// Files for base recalibration
ref_fasta_dict       = file("${ref_base}/${params.references[ref_build].dict}", checkIfExists: true)
dbsnp                = file("${ref_base}/${params.references[ref_build].dbsnp}", checkIfExists: true)
dbsnp_index          = file("${ref_base}/${params.references[ref_build].dbsnpIndex}", checkIfExists: true)
known_indels         = Channel.fromPath("${ref_base}/${params.references[ref_build].knownIndels}", checkIfExists: true).collect()
known_indels_index   = Channel.fromPath("${ref_base}/${params.references[ref_build].knownIndelsIndex}", checkIfExists: true).collect()

// annotation
snpeff_db            = Channel.fromPath("${ref_base}/${params.references[ref_build].snpEffDb}").collect()
snpeff_config        = file("${ref_base}/${params.references[ref_build].snpEffConfig}")
vcfanno_config_files = Channel.fromPath("${assay_base}/${params.assays[assay].annotationConfig}", checkIfExists: true).collect()
vcfanno_ref_files    = Channel.fromPath("${ref_base}/Homo_sapiens/GATK/GRCh38/Annotation/{clinvar,gnomAD,OMIM,CPDX,HGMD}/*").collect()

// xls report generation
xls_config           = file("${assay_base}/${params.assays[assay].xlsConfig}", checkIfExists: true)


// Print inputs/config
println("Sample: " + sample)
println("Assay: " + assay)
println("Reference: " + ref_build)
println("Reference Base: " + ref_base)
println("Assay Base: " + assay_base)
println("Sample Base: " + sample_base)
println('XLS Config: ' + xls_config)


Channel.fromPath("${params.sample_base}/${sample}/exome/libraries/**.fastq.gz")
    .map { fastq -> 
        def (filename, readgroup_id, library_id, _f1, _f2, sample_id, rest) = fastq.toString().tokenize('/').reverse() // tokenize path
        return [readgroup_id, fastq]
    }
    .groupTuple()
    .map { readgroup_id, fastqs ->
        def (filename, _r, library_id, _f1, _f2, sample_id, rest) = fastqs[0].toString().tokenize('/').reverse() // tokenize path
        def (fcid, lane, barcodes) = readgroup_id.tokenize(".")
        def config = [
            sample_id: sample_id,
            library_id: library_id,
            readgroup_id: readgroup_id,
            fcid: fcid, lane: lane, barcodes: barcodes,
        ]
        return tuple(readgroup_id, config, fastqs)
    }
    .set { mapping_source_fastqs_ch }


// alignment of individual read groups (and sort by *queryname*)
process map_reads {
    label 'alignment'
    tag "${sample_id}:${config.readgroup_id}"
    echo true
    input:
        tuple sample_id, config, file(fastqs) from mapping_source_fastqs_ch
        path ref_fasta
        path bwa_index

    output:
        tuple sample_id, file("${config.readgroup_id}.bam") into readgroup_bams_ch

    script:
        // sample_id (SM): YY-XXXXX
        // library_id (LB): {sample_id}-{library_id}
        // platform_unit (PU): {FCID}.{Lane}
        // readgroup_id (ID): {FCID}.{Lane}.{barcode-barcode2}
        sample_id = config.sample_id
        platform_unit = "${config.fcid}.${config.lane}"
        library_id = config.library_id
        readgroup_string = "@RG\\tID:${config.readgroup_id}\\tSM:${sample_id}\\tLB:${config.library_id}\\tPL:Illumina\\tPU:${platform_unit}"
        def (fastq1, fastq2) = fastqs
        """
        bwa mem \
        -R "${readgroup_string}" \
        -K 100000000 \
        -t ${task.cpus}  \
        ${ref_fasta} *.fastq.gz 2> log.txt \
        | samtools sort -n -m4G \
          - -o ${config.readgroup_id}.bam

        """
}

// merge readgroups together, mark duplicates, and output coordinate-sorted bam
process merge_and_markdups {
    label 'gatk_merge_and_markdups'
    tag "${sample_id}"

    input:
        tuple sample_id, file(bams) from readgroup_bams_ch.groupTuple()
    output:
        tuple sample_id, file("${sample_id}.sorted.deduped.bam"), file("${sample_id}.sorted.deduped.bam.bai") into merged_deduped_bams_ch
    script:
        inputs = bams.collect{"-I ${it}"}.join(' ')

        """
        gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        MarkDuplicatesSpark \
        ${inputs} \
        -O ${sample_id}.sorted.deduped.bam \
        -M mark_duplicates.metrics \
        --tmp-dir . \
        --create-output-bam-index true \
        --spark-master local[${task.cpus}]

        # SetNmMdAndUqTags?
        """
}



// recalibrate bases with GATK

process recalibrate_bases {
    label 'gatk_base_recalibrator'
    tag "${sample_id}"

    input:
        tuple sample_id, file(bam), file(bai) from merged_deduped_bams_ch
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path dbsnp
        path dbsnp_index
        path known_indels
        path known_indels_index

    output:
        tuple sample_id, file("${sample_id}.recalibrated.bam"), file("${sample_id}.recalibrated.bai") into recalibrated_bams_ch
   
    script:
        
        // carry over NCGL options
        // --deletions_default_quality 45 \
        // --indels_context_size 3 \
        // --insertions_default_quality 45 \
        // --low_quality_tail 2 \
        // --maximum_cycle_value 500 \
        // --mismatches_context_size 2 \
        // --mismatches_default_quality -1 \

        known_sites_all_indels = known_indels.collect{"--known-sites ${it}"}.join(' ')
        """
        gatk --java-options "-Xmx${task.memory.toGiga()}g" \
            BaseRecalibrator \
            --reference ${ref_fasta} \
            --input ${bam} \
            --known-sites ${dbsnp} \
            ${known_sites_all_indels} \
            --output ${sample_id}.recal_table 2> log.txt

        gatk --java-options "-Xmx${task.memory.toGiga()}g" \
            ApplyBQSR \
            --reference ${ref_fasta} \
            --input ${bam} \
            --bqsr-recal-file ${sample_id}.recal_table \
            --create-output-bam-index \
            --output ${sample_id}.recalibrated.bam 2>> log.txt
        """
}



recalibrated_bams_ch.into { haplotype_caller_input; strelka2_input }

process haplotype_caller {
    label 'gatk_haplotype_caller'
    tag "${sample_id}"

    input:
        tuple sample_id, file(bam), file(bai) from  haplotype_caller_input
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path dbsnp
        path dbsnp_index
        path "intervals/*" from assay_intervals
        // path known_indels
        // path known_indels_index
    output:
        tuple sample_id, file("${sample_id}.gatk.g.vcf.gz") into genotype_gvcf_in

    script:
        """
        gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${bam} \
        --dbsnp ${dbsnp} \
        --intervals intervals/*.bed \
        -O ${sample_id}.gatk.g.vcf.gz \
        -ERC GVCF
        """
}


// TODO: strelka denovo.py (multisample)
process strelka2 {
    label 'strelka2'
    tag "${sample_id}"

    input:
        tuple sample_id, file(bam), file(bai) from strelka2_input
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path "intervals/*" from assay_intervals

    output:
        tuple val(sample_id), file("${sample_id}.strelka.vcf.gz"), file("${sample_id}.strelka.vcf.gz.tbi"), val("strelka") into strelka2_out

    script:
        """
        configureStrelkaGermlineWorkflow.py \
            --bam ${bam} \
            --referenceFasta ${ref_fasta} \
            --exome --callRegions intervals/*.gz \
            --runDir Strelka

        python Strelka/runWorkflow.py -m local -j ${task.cpus}

        mv Strelka/results/variants/genome.*.vcf.gz     ${sample_id}.strelka.g.vcf.gz
        mv Strelka/results/variants/genome.*.vcf.gz.tbi ${sample_id}.strelka.g.vcf.gz.tbi
        mv Strelka/results/variants/variants.vcf.gz     ${sample_id}.strelka.vcf.gz
        mv Strelka/results/variants/variants.vcf.gz.tbi ${sample_id}.strelka.vcf.gz.tbi
        """
}

// TODO: GenomicsDB import step here

process genotype_gvcf {
    label 'gatk_genotype_gvcf'
    tag "${sample_id}"

    input:
        tuple sample_id, file(gvcf) from genotype_gvcf_in
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
        path dbsnp
        path dbsnp_index
        path "intervals/*" from assay_intervals

    output:
        tuple val(sample_id), file("${sample_id}.gatk.vcf.gz"), file("${sample_id}.gatk.vcf.gz.tbi"), val("gatk") into genotype_gvcf_out

    script:
        """
        gatk --java-options -Xmx${task.memory.toGiga()}g \
            IndexFeatureFile \
            -I ${gvcf}
        gatk --java-options -Xmx${task.memory.toGiga()}g \
            GenotypeGVCFs \
            -R ${ref_fasta} \
            --intervals intervals/*.bed \
            --dbsnp ${dbsnp} \
            -V ${gvcf} \
            -O ${sample_id}.gatk.vcf.gz
        gatk --java-options -Xmx${task.memory.toGiga()}g \
            IndexFeatureFile \
            -I ${sample_id}.gatk.vcf.gz
        """
}
normalize_vcf_in = genotype_gvcf_out.mix(strelka2_out)

process normalize_vcf {
    label 'vcf_normalize'
    tag "${sample_id}"
    
    input:
        tuple val(sample_id), file(vcf), file(ix), val(variant_caller) from normalize_vcf_in
        path ref_fasta
        path ref_fasta_fai
        path ref_fasta_dict
    output:
        tuple sample_id, file("${vcf.baseName}.normalized.vcf.gz"), file("${vcf.baseName}.normalized.vcf.gz.tbi"), variant_caller into normalize_vcf_out
    shell:
        """
        bcftools norm ${vcf} --fasta-ref ${ref_fasta} \
            --check-ref=w \
            --multiallelics=-both \
            --output-type=z -o ${vcf.baseName}.normalized.vcf.gz
        bcftools index --tbi ${vcf.baseName}.normalized.vcf.gz
        """
}

//process intersect_and_merge {
    //// intersect GATK + Strelka variants -- take strelka passing which are not in GATK
    //// create separate file or merge into merged file with GATK?
//
//}


process annotation {
    // will use combined docker container with snpeff/snpsift, vcf_anno, +/- slivar
    label 'annotation'
    tag "${sample_id}"

    publishDir "publish/${sample_id}/", mode: 'copy'
    input:
        tuple val(sample_id), file(vcf), file(ix), val(variant_caller) from normalize_vcf_out
        path snpeff_config from snpeff_config
        path "snpeff/*" from snpeff_db
        path vcfanno_config_files
        path vcfanno_ref_files
    output:
        tuple val(sample_id), file("${sample_id}.${variant_caller}.annotated.vcf"), val(variant_caller) into annotated_vcf_out

    shell:
    """
    # basic function annotation (ANN field)
    java -Xmx4g -jar /snpEff/snpEff.jar \
        -c snpeff.config \
        -dataDir data/ \
        -noStats \
        -noDownload \
        -no-intergenic -no-upstream -no-downstream \
        -noMotif -noNextProt \
        GRCh38.86 ${vcf} > ${sample_id}.snpeff.vcf

    cat ${sample_id}.snpeff.vcf \
    | vcfanno_linux64 -lua custom.lua gnomad_annotation.toml /dev/stdin \
    | vcfanno_linux64 -lua custom.lua clinvar_annotation.toml /dev/stdin \
    | vcfanno_linux64 -lua custom.lua omim_annotation.toml /dev/stdin \
    | vcfanno_linux64 -lua custom.lua cpdx_annotation.toml /dev/stdin \
    | vcfanno_linux64 -lua custom.lua hgmd_annotation.toml /dev/stdin \
    > ${sample_id}.${variant_caller}.annotated.vcf
    """
}

process make_xls {
    label 'annotation'
    tag "${sample_id}"

    publishDir "publish/${sample_id}", mode: 'copy'

    input:
        tuple val(sample_id), file(vcf), val(variant_caller) from annotated_vcf_out
        path config from xls_config
    output:
        tuple sample_id, file("${sample_id}.${variant_caller}.report.xls")

    shell:
    """
    vcf2xlsx.py ${vcf} ${config} ${sample_id}.${variant_caller}.report.xlsx
    """
}
//process slivar {
    // https://github.com/brentp/slivar
    // https://brentp.github.io/post/trio-duo-solo/
    // https://brentp.github.io/post/variant-filter/
//}