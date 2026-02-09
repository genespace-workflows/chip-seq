nextflow.enable.dsl=2

params.samplesheet   = params.samplesheet ?: "samplesheet.csv"
params.outdir        = params.outdir ?: "results"

// Trimmomatic
params.trimmomatic_jar = params.trimmomatic_jar ?: "trimmomatic.jar"
params.adapters_fa     = params.adapters_fa ?: "adapters.fa"
params.trim_threads    = params.trim_threads ?: 4

// Bowtie2
params.bt2_index     = params.bt2_index ?: null   // prefix of Bowtie2 index, required
params.align_threads = params.align_threads ?: 8

// PhantomPeakQualTools / SPP
params.spp_script    = params.spp_script ?: "run_spp.R"  // usually provided by PhantomPeakQualTools
params.spp_threads   = params.spp_threads ?: 4

// Peak calling inputs
params.mnase_input_bam = params.mnase_input_bam ?: null  // required
params.bg_mouse_bam    = params.bg_mouse_bam ?: null     // optional, if you really need an extra background

// Blacklist
params.blacklist_bed = params.blacklist_bed ?: "ENCFF547MET.bed" // DAC blacklist (BED)

// Genome size for MACS2 (mm)
params.macs2_gsize   = params.macs2_gsize ?: "mm"

// DESeq2
params.deseq2_r      = params.deseq2_r ?: "${projectDir}/scripts/deseq2.R"

workflow {

    if( !params.bt2_index )
        error "Parameter --bt2_index is required (Bowtie2 index prefix)."

    if( !params.mnase_input_bam )
        error "Parameter --mnase_input_bam is required (MNase-treated input BAM)."

    // Read samplesheet
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample = row.sample as String
            def group  = row.group  as String
            def fastq  = file(row.fastq)
            tuple(sample, group, fastq)
        }
        .set { CH_SAMPLES }

    // 1) Trimming
    TRIM(CH_SAMPLES)
        .set { CH_TRIMMED }  // (sample, group, trimmed_fastq)

    // 2) Alignment -> sorted BAM + BAI
    ALIGN(CH_TRIMMED)
        .set { CH_BAMS }     // (sample, group, bam, bai)

    // 3) QC metrics NSC/RSC (PhantomPeakQualTools)
    QC_SPP(CH_BAMS)
        .set { CH_QC }       // (sample, qc_table, qc_plot, qc_log)

    // 4) Pooled BAM for MACS2 (pooled ChIP BAMs)
    MERGE_BAMS(CH_BAMS.map{ s,g,bam,bai -> tuple("pooled_chip", bam) })
        .set { CH_POOLED }   // (name, pooled_bam, pooled_bai)

    // 5) MACS2 peak calling with single-nucleosome parameters
    //    Uses pooled_chip as -t and MNase input as -c.
    //    If you need extra background mouse BAM, you can incorporate it into control before run.
    MACS2_CALLPEAK(CH_POOLED.map{ name,bam,bai -> bam })
        .set { CH_MACS2 }    // (nuc_peaks_bed, summits_bed, macs2_dir)

    // 6) Filter peaks by blacklist (DAC)
    FILTER_BLACKLIST(CH_MACS2.map{ peaks, summits, dir -> peaks })
        .set { CH_PEAKS_FILTERED } // (filtered_peaks_bed)

    // 7) Count reads in peaks per sample (MAPQ>=10)
    COUNT_IN_PEAKS(CH_BAMS, CH_PEAKS_FILTERED)
        .set { CH_COUNTS }   // (count_matrix_tsv)

    // 8) DESeq2 differential analysis
    DESEQ2(CH_COUNTS)
        .set { CH_DE }       // (deseq2_results_tsv, normalized_counts_tsv)

    // Publish QC summary etc. (optional)
    CH_QC.view()
    CH_DE.view()
}

process TRIM {
    tag "${sample}"
    publishDir "${params.outdir}/01_trimmed", mode: 'copy'

    input:
    tuple val(sample), val(group), path(fastq)

    output:
    tuple val(sample), val(group), path("${sample}.trimmed.fastq.gz")

    cpus params.trim_threads

    script:
    """
    java -jar ${params.trimmomatic_jar} SE -threads ${task.cpus} -phred33 \
      ${fastq} ${sample}.trimmed.fastq.gz \
      ILLUMINACLIP:${params.adapters_fa}:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process ALIGN {
    tag "${sample}"
    publishDir "${params.outdir}/02_align", mode: 'copy'

    input:
    tuple val(sample), val(group), path(fq)

    output:
    tuple val(sample), val(group),
          path("${sample}.sorted.bam"),
          path("${sample}.sorted.bam.bai")

    cpus params.align_threads

    script:
    """
    set -euo pipefail

    bowtie2 -p ${task.cpus} -x ${params.bt2_index} -U ${fq} \
      | samtools view -bS - \
      | samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam

    samtools index ${sample}.sorted.bam
    """
}

process QC_SPP {
    tag "${sample}"
    publishDir "${params.outdir}/03_qc_spp", mode: 'copy'

    input:
    tuple val(sample), val(group), path(bam), path(bai)

    output:
    tuple val(sample),
          path("${sample}.spp.txt"),
          path("${sample}.spp.pdf"),
          path("${sample}.spp.log")

    cpus params.spp_threads

    script:
    """
    set -euo pipefail

    # PhantomPeakQualTools typically runs via run_spp.R
    # Output format may differ depending on your version; adjust flags if needed.
    Rscript ${params.spp_script} \
      -c=${bam} -p=${task.cpus} -savp \
      -out=${sample}.spp.txt \
      > ${sample}.spp.log 2>&1

    # Many versions write PDF as 'phantompeakqualtools.pdf' or similar; standardize if exists
    if ls *.pdf 1>/dev/null 2>&1; then
      pdf=\$(ls *.pdf | head -n 1)
      mv "\$pdf" ${sample}.spp.pdf
    else
      # create placeholder to satisfy output (or remove from outputs)
      echo "No PDF produced by SPP" > ${sample}.spp.pdf
    fi
    """
}

process MERGE_BAMS {
    tag "${name}"
    publishDir "${params.outdir}/04_pooled", mode: 'copy'

    input:
    tuple val(name), path(bams)

    output:
    tuple val(name),
          path("${name}.merged.bam"),
          path("${name}.merged.bam.bai")

    script:
    """
    set -euo pipefail

    samtools merge -f ${name}.merged.bam ${bams}
    samtools sort -o ${name}.merged.bam ${name}.merged.bam
    samtools index ${name}.merged.bam
    """
}

process MACS2_CALLPEAK {
    tag "macs2"
    publishDir "${params.outdir}/05_macs2", mode: 'copy'

    input:
    path(pooled_chip_bam)

    output:
    tuple path("single_nuc_peaks.bed"), path("summits.bed"), path("macs2_out")

    script:
    """
    set -euo pipefail
    mkdir -p macs2_out

    # MACS2 callpeak; single-nucleosome parameters
    macs2 callpeak \
      -t ${pooled_chip_bam} \
      -c ${params.mnase_input_bam} \
      -f BAM -g ${params.macs2_gsize} \
      -n pooled_singleNuc \
      --shift 37 --extsize 73 \
      --outdir macs2_out

    # Get summits
    # Usually: pooled_singleNuc_summits.bed
    cp macs2_out/pooled_singleNuc_summits.bed summits.bed

    # Define single nucleosome peak as 147bp around summit:
    # [summit-73, summit+74) in BED
    awk 'BEGIN{OFS="\\t"}{
      s=$2-73; e=$2+74;
      if(s<0) s=0;
      print $1, s, e, $4, $5, $6
    }' summits.bed > single_nuc_peaks.bed
    """
}

process FILTER_BLACKLIST {
    tag "blacklist"
    publishDir "${params.outdir}/06_blacklist_filtered", mode: 'copy'

    input:
    path(peaks_bed)

    output:
    path("single_nuc_peaks.filtered.bed")

    script:
    """
    set -euo pipefail

    # Exclude peaks overlapping DAC blacklist
    bedtools intersect -v \
      -a ${peaks_bed} \
      -b ${params.blacklist_bed} \
      > single_nuc_peaks.filtered.bed
    """
}

process COUNT_IN_PEAKS {
    tag "count"
    publishDir "${params.outdir}/07_counts", mode: 'copy'

    input:
    tuple val(sample), val(group), path(bam), path(bai) from CH_BAMS
    path(peaks_bed) from CH_PEAKS_FILTERED

    output:
    tuple val(sample), val(group), path("${sample}.counts.tsv")

    script:
    """
    set -euo pipefail

    # filter MAPQ>=10
    samtools view -b -q 10 ${bam} > ${sample}.q10.bam
    samtools index ${sample}.q10.bam

    # counts per peak (one file per sample)
    # output: chr start end name score strand  +  count in last column
    bedtools coverage -counts -a ${peaks_bed} -b ${sample}.q10.bam \
      | awk 'BEGIN{OFS="\\t"}{print $1,$2,$3,NR,$NF}' \
      > ${sample}.counts.tsv

    rm -f ${sample}.q10.bam ${sample}.q10.bam.bai
    """
}

process DESEQ2 {
    tag "deseq2"
    publishDir "${params.outdir}/08_deseq2", mode: 'copy'

    input:
    tuple val(meta), path(files) from COUNT_IN_PEAKS.out.collect()

    output:
    path("deseq2_results.tsv")
    path("normalized_counts.tsv")

    script:
    // meta is a list of [sample, group] pairs; files is list of *.counts.tsv paths in same order
    def rows = (0..<meta.size()).collect { i ->
        def s = meta[i][0]
        def g = meta[i][1]
        def f = files[i]
        "${s}\t${g}\t${f}"
    }.join("\n")

    """
    set -euo pipefail

    echo -e "sample\\tgroup\\tcounts" > samples.tsv
    cat >> samples.tsv << 'EOF'
    ${rows}
    EOF

    Rscript ${params.deseq2_r} \
      --samples samples.tsv \
      --out_results deseq2_results.tsv \
      --out_norm normalized_counts.tsv
    """
}

workflow.onComplete {
    println "Done. Results in: ${params.outdir}"
}
