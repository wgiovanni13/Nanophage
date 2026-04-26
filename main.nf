#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ===========================================
//               PARAMETERS
// ===========================================

params.file_id     = '15trDWYRK9HhbywG2a6NuHeCCcolieNtG'
params.filename    = '1_024_O.fastq.gz'
params.outDir      = 'results'
params.min_len     = 500
params.min_qual    = 10

// ===========================================
//                  WORKFLOW
// ===========================================

workflow {
    reads_ch = DOWNLOAD()
    READ_QC(reads_ch)
    filtered_ch = FILTER_READS(reads_ch)
    assembly_ch = ASSEMBLE(filtered_ch.reads)
    ASSEMBLY_QC(assembly_ch.assembly, filtered_ch.reads)
    GENOMAD(assembly_ch.assembly)
    BLAST_NT(assembly_ch.assembly)
    BACPHLIP(assembly_ch.assembly)
    TRNASCAN(assembly_ch.assembly)
    ABRICATE(assembly_ch.assembly)
    PHAROKKA(assembly_ch.assembly)

    REPORT(
        READ_QC.out.stats,
        FILTER_READS.out.stats,
        ASSEMBLE.out.info,
        ASSEMBLY_QC.out.quast_tsv,
        ASSEMBLY_QC.out.coverage,
        ASSEMBLY_QC.out.depth,
        ASSEMBLY_QC.out.checkv_summary,
        GENOMAD.out.summary,
        GENOMAD.out.genes,
        BLAST_NT.out.hits,
        BACPHLIP.out.results,
        TRNASCAN.out.results,
        ABRICATE.out.summary,
        PHAROKKA.out.plot
    )

}

// ===========================================
//                PROCESSES
// ===========================================

// 1. DOWNLOAD RAW READS FROM GOOGLE DRIVE

process DOWNLOAD {
    publishDir "${params.outDir}/00_raw", mode: 'copy'
    conda "conda-forge::python=3.12 conda-forge::gdown=5.2.0"
    cpus 1
    memory '4.GB'
    time '1h'

    output:
    path "${params.filename}", emit: reads

    script:
    """
    gdown ${params.file_id} -O ${params.filename}
    """
}

// 2. RAW READ QUALITY ASSESSMENT

process READ_QC {
    tag "${reads.simpleName}"
    publishDir "${params.outDir}/01_read_qc", mode: 'copy'
    conda "bioconda::nanoplot=1.43.0 conda-forge::python=3.12 conda-forge::python-kaleido=0.2.1"
    cpus 4
    memory '16.GB'
    time '2h'

    input:
    path reads

    output:
    path "nanoplot/NanoStats.txt", emit: stats
    path "nanoplot/"

    script:
    """
    NanoPlot --fastq ${reads} -o nanoplot -t ${task.cpus} --tsv_stats --no_static
    """
}

// 3. FILTER LOW QUALITY AND SHORT READS

process FILTER_READS {
    tag "${reads.simpleName}"
    publishDir "${params.outDir}/02_filtered", mode: 'copy'
    conda "bioconda::chopper=0.9.0 bioconda::seqkit=2.9.0"
    cpus 4
    memory '8.GB'
    time '1h'

    input:
    path reads

    output:
    path "filtered.fastq.gz",    emit: reads
    path "filtering_report.txt", emit: stats

    script:
    """

    # Pre filter statistics

    seqkit stats -T ${reads} > pre_stats.tsv
    PRE_READS=\$(tail -1 pre_stats.tsv | cut -f4)
    PRE_BASES=\$(tail -1 pre_stats.tsv | cut -f5)

    # Filter: remove reads < ${params.min_len} (500) bp and < Q${params.min_qual} (Q10)

    gunzip -c ${reads} \
        | chopper -q ${params.min_qual} -l ${params.min_len} -t ${task.cpus} \
        | gzip > filtered.fastq.gz

    # Post filter statistics

    seqkit stats -T filtered.fastq.gz > post_stats.tsv
    POST_READS=\$(tail -1 post_stats.tsv | cut -f4)
    POST_BASES=\$(tail -1 post_stats.tsv | cut -f5)

    # Calculate removed reads

    REMOVED=\$(( \$PRE_READS - \$POST_READS ))

    cat <<EOF > filtering_report.txt
    Filtering parameters: min_length=${params.min_len} bp, min_quality=Q${params.min_qual}
    Before: \${PRE_READS} reads, \${PRE_BASES} bases
    After:  \${POST_READS} reads, \${POST_BASES} bases
    Removed: \${REMOVED} reads
    EOF
    """
}

// 4. GENOME ASSEMBLY

process ASSEMBLE {
    tag "flye"
    publishDir "${params.outDir}/03_assembly", mode: 'copy'
    conda "bioconda::flye=2.9.5"
    cpus 8
    memory '16.GB'
    time '2h'

    input:
    path filtered_reads

    output:
    path "flye_out/assembly.fasta",    emit: assembly
    path "flye_out/assembly_info.txt", emit: info

    script:
    """
    flye --nano-raw ${filtered_reads} \
        -o flye_out \
        -t ${task.cpus}
    """
}

// 5. ASSEMBLY QUALITY ASSESSMENT

process ASSEMBLY_QC {
    tag "assembly_qc"
    publishDir "${params.outDir}/04_assembly_qc", mode: 'copy'
    conda "bioconda::quast=5.2.0 bioconda::minimap2=2.28 bioconda::samtools=1.21 bioconda::checkv=1.0.3"
    cpus 4
    memory '16.GB'
    time '4h'

    input:
    path assembly
    path filtered_reads

    output:
    path "quast/report.tsv",           emit: quast_tsv
    path "coverage.tsv",               emit: coverage
    path "depth.tsv",                  emit: depth
    path "checkv/quality_summary.tsv", emit: checkv_summary
    path "quast/",                     emit: quast_dir
    path "checkv/",                    emit: checkv_dir

    script:
    """

    # Assembly statistics, coverage and GC percentage with QUAST, Minimap2 and CheckV

    quast ${assembly} -o quast --min-contig 500 -t ${task.cpus}

    minimap2 -a -x map-ont -t ${task.cpus} ${assembly} ${filtered_reads} \
        | samtools sort -@ ${task.cpus} -o mapped.bam
    samtools index mapped.bam
    samtools coverage mapped.bam > coverage.tsv
    samtools depth -a mapped.bam > depth.tsv

    checkv download_database checkv_db
    checkv end_to_end ${assembly} checkv -t ${task.cpus} -d checkv_db/checkv-db-v1.5
    """
}

// 6. CONTIG CHARACTERIZATION: taxonomy, lifestyle, gene content, tRNAs, AMR and phage functional annotation

// 6.1. Viral identification, taxonomy, gene annotation

process GENOMAD {
    tag "genomad"
    publishDir "${params.outDir}/05_characterization/genomad", mode: 'copy'
    conda "bioconda::genomad=1.8.1"
    cpus 8
    memory '32.GB'
    time '4h'

    input:
    path assembly

    output:
    path "results/*_summary/*_virus_summary.tsv", emit: summary
    path "results/*_summary/*_virus_genes.tsv",   emit: genes
    path "results/",                               emit: full

    script:
    """
    mkdir -p db
    genomad download-database db
    genomad end-to-end --cleanup ${assembly} results db/genomad_db -t ${task.cpus}
    """
}

// 6.2. Closest relative via remote BLASTn

process BLAST_NT {
    tag "blast"
    publishDir "${params.outDir}/05_characterization/blast", mode: 'copy'
    conda "bioconda::blast"
    cpus 1
    memory '4.GB'
    time '2h'

    input:
    path assembly

    output:
    path "blast_hits.tsv", emit: hits

    script:
    """
    blastn -query ${assembly} -db nt -remote \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -max_target_seqs 5 -evalue 1e-10 \
        -out blast_hits.tsv \
        || echo "no_hits" > blast_hits.tsv
    """
}

// 6.3. Lifestyle prediction

process BACPHLIP {
    tag "bacphlip"
    publishDir "${params.outDir}/05_characterization/bacphlip", mode: 'copy'
    conda "bioconda::bacphlip=0.9.6 bioconda::hmmer=3.3.2 conda-forge::python=3.8 conda-forge::scikit-learn=0.23.1 conda-forge::numpy=1.23.5"
    cpus 1
    memory '4.GB'
    time '1h'

    input:
    path assembly

    output:
    path "*.bacphlip", emit: results

    script:
    """
    bacphlip -i ${assembly}
    """
}

// 6.4. tRNA detection

process TRNASCAN {
    tag "trnascan"
    publishDir "${params.outDir}/05_characterization/trnascan", mode: 'copy'
    conda "bioconda::trnascan-se"
    cpus 1
    memory '4.GB'
    time '1h'

    input:
    path assembly

    output:
    path "trna_results.txt", emit: results

    script:
    """
    tRNAscan-SE -B -o trna_results.txt ${assembly} || touch trna_results.txt
    """
}

// 6.5. AMR screening

process ABRICATE {
    tag "abricate"
    publishDir "${params.outDir}/05_characterization/abricate", mode: 'copy'
    conda "bioconda::abricate"
    cpus 1
    memory '4.GB'
    time '1h'

    input:
    path assembly

    output:
    path "amr_summary.tsv", emit: summary

    script:
    """
    abricate --db ncbi ${assembly} > amr_ncbi.tsv || true
    abricate --db card ${assembly} > amr_card.tsv || true
    abricate --db vfdb ${assembly} > amr_vfdb.tsv || true
    abricate --summary amr_ncbi.tsv amr_card.tsv amr_vfdb.tsv > amr_summary.tsv
    """
}

// 6.6. Phage functional annotation and genome map

process PHAROKKA {
    tag "pharokka"
    publishDir "${params.outDir}/05_characterization/pharokka", mode: 'copy'
    conda "bioconda::pharokka"
    cpus 8
    memory '16.GB'
    time '4h'

    input:
    path assembly

    output:
    path "pharokka_out/pharokka.gff",    emit: gff
    path "pharokka_out/pharokka.gbk",    emit: gbk
    path "pharokka_plot/pharokka_plot.png",               emit: plot
    path "pharokka_out/",                emit: full

    script:
    """
    install_databases.py -o pharokka_db
    pharokka.py -i ${assembly} -o pharokka_out -d pharokka_db -t ${task.cpus} --fast -f
    mkdir -p pharokka_plot
    pharokka_plotter.py -i ${assembly} \
        --gff pharokka_out/pharokka.gff \
        --genbank pharokka_out/pharokka.gbk \
        -o pharokka_plot -f \
        || touch pharokka_plot/pharokka_plot.png
    """
}

// 7. GENERATE HTML REPORT WITH EMBEDDED FIGURES FROM PREVIOUS PROCESSES

process REPORT {
    tag "report"
    publishDir "${params.outDir}/06_report", mode: 'copy'
    conda "conda-forge::python=3.12 conda-forge::matplotlib=3.9 conda-forge::pandas"
    cpus 1
    memory '4.GB'
    time '1h'

    input:
    path raw_stats
    path filter_stats
    path assembly_info
    path quast_tsv
    path coverage
    path depth
    path checkv_summary
    path genomad_summary
    path genomad_genes
    path blast_hits
    path bacphlip
    path trna
    path amr
    path pharokka_plot

    output:
    path "nanophage_report.html", emit: report

    script:
    """
    generate_report.py \
        --raw-stats ${raw_stats} \
        --filter-stats ${filter_stats} \
        --assembly-info ${assembly_info} \
        --quast ${quast_tsv} \
        --coverage ${coverage} \
        --depth ${depth} \
        --checkv ${checkv_summary} \
        --genomad-summary ${genomad_summary} \
        --genomad-genes ${genomad_genes} \
        --blast ${blast_hits} \
        --bacphlip ${bacphlip} \
        --trna ${trna} \
        --amr ${amr} \
        --pharokka-plot ${pharokka_plot} \
        -o nanophage_report.html
    """
}
