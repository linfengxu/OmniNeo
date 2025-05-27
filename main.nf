#!/usr/bin/env nextflow

// define help information
def helpMessage() {
    log.info"""
    =============================================================================================
                                OmniNeo PIPELINE v1.0
    =============================================================================================

    Usage:
    nextflow run main.nf [options]

    Main Parameters:
      --samplesheet        Path to sample information table (TSV format, required)
      --outdir             Output directory [default: ${params.outdir ?: 'results'}]
      --hg38db             Reference genome database path [default: ${params.hg38db ?: '/path/to/hg38db'}]
      --tmp_dir            Temporary directory [default: ${params.tmp_dir ?: '/tmp'}]
      --refseq_peptide     RefSeq peptide reference file [default: ${params.refseq_peptide ?: 'NA'}]
      --dbvep_dir          VEP archive file path [default: ${params.dbvep_dir ?: 'NA'}]

    Resource Parameters:
      --max_memory         Maximum memory per job [default: ${params.max_memory ?: '64.GB'}]
      --max_cpus           Maximum CPU cores per job [default: ${params.max_cpus ?: '16'}]
      --max_time           Maximum run time per job [default: ${params.max_time ?: '24.h'}]

    Fastp Parameters:
      --fastp.threads             Threads for fastp [default: ${params.fastp?.threads ?: '30'}]
      --fastp.adapters_auto       Auto-detect adapters [default: ${params.fastp?.adapters_auto ?: 'true'}]
      --fastp.min_length          Minimum read length [default: ${params.fastp?.min_length ?: '50'}]
      --fastp.cut_mean_quality    Quality threshold [default: ${params.fastp?.cut_mean_quality ?: '20'}]

    BWA Parameters:
      --bwa.threads               Threads for BWA [default: ${params.bwa?.threads ?: '16'}]
      --bwa.container             Container for BWA [default: ${params.bwa?.container ?: 'NA'}]

    GATK Parameters:
      --gatk.container                   GATK container path [default: ${params.gatk?.container ?: 'broadinstitute/gatk:latest'}]
      --gatk.known_sites.dbsnp           Known SNP sites file [default: ${params.gatk?.known_sites?.dbsnp ?: 'NA'}]
      --gatk.known_sites.mills_1000g     Mills & 1000G known indels file [default: ${params.gatk?.known_sites?.mills_1000g ?: 'NA'}]
      --gatk.read_groups.id              Read group ID [default: ${params.gatk?.read_groups?.id ?: '4'}]
      --gatk.read_groups.library         Read group library [default: ${params.gatk?.read_groups?.library ?: 'lib1'}]
      --gatk.read_groups.platform        Read group platform [default: ${params.gatk?.read_groups?.platform ?: 'illumina'}]

    STAR-Fusion Parameters:
      --star_fusion.genome_lib_dir       STAR-Fusion genome library directory [default: ${params.star_fusion?.genome_lib_dir ?: 'NA'}]
      --star_fusion.container            STAR-Fusion container [default: ${params.star_fusion?.container ?: 'NA'}]

    Kallisto Parameters:
      --kallisto.transcriptome_fasta     Transcriptome FASTA file [default: ${params.kallisto?.transcriptome_fasta ?: 'NA'}]
      --kallisto.index_file              Kallisto index file [default: ${params.kallisto?.index_file ?: 'NA'}]
      --kallisto.bootstrap_samples       Number of bootstrap samples [default: ${params.kallisto?.bootstrap_samples ?: '100'}]
      --kallisto.container               Kallisto container [default: ${params.kallisto?.container ?: 'NA'}]

    OptiType Parameters:
      --optitype.hla_reference_rna       HLA reference file for RNA [default: ${params.optitype?.hla_reference_rna ?: 'NA'}]
      --optitype.container               OptiType container [default: ${params.optitype?.container ?: 'NA'}]

    Other Tool Parameters:
      --bcftools.container               BCFtools container [default: ${params.bcftools?.container ?: 'NA'}]
      --samtools.container               Samtools container [default: ${params.samtools?.container ?: 'NA'}]
      --picard.container                 Picard container [default: ${params.picard?.container ?: 'NA'}]
      --generate_pep.container           Peptide generation container [default: ${params.generate_pep?.container ?: 'NA'}]
      --generate_pep.threads             Peptide generation threads [default: ${params.generate_pep?.threads ?: '4'}]

    Other Options:
      --help                  Display this help message and exit
      --skip_dna              Skip DNA analysis workflow
      --skip_rna              Skip RNA analysis workflow
      --skip_combine_ms       Skip combining peptides for MS analysis
      --skip_maxquant         Skip MaxQuant analysis
      --skip_netmhcpan        Skip netMHCpan analysis

    Sample Sheet Format (TSV):
    sampleID  DNANormalReads1  DNANormalReads2  DNATumorReads1  DNATumorReads2  RNANormalReads1  RNANormalReads2  RNATumorReads1  RNATumorReads2  MSTumor
    L041      path/to/L041_dna_normal_R1.fq.gz  path/to/L041_dna_normal_R2.fq.gz  path/to/L041_dna_tumor_R1.fq.gz  path/to/L041_dna_tumor_R2.fq.gz  path/to/L041_rna_normal_R1.fq.gz  path/to/L041_rna_normal_R2.fq.gz  path/to/L041_rna_tumor_R1.fq.gz  path/to/L041_rna_tumor_R2.fq.gz  path/to/L041_ms1.raw,path/to/L041_ms2.raw
    
    Note: The MSTumor column contains paths to MS raw files for MaxQuant analysis, comma-separated for multiple files.
    
    Configuration:
    A detailed configuration file is available at 'nextflow.config'. Advanced users can modify parameters directly in this file.
    
    Example:
    nextflow run main.nf --samplesheet samples.tsv --outdir results --hg38db /path/to/hg38 -profile singularity
    
    Using Default Parameters:
    If you run the pipeline without specifying parameters, it will use all default values from nextflow.config:
    nextflow run main.nf -profile singularity
    
    Pipeline Info:
    For more information and bug reports, please visit: https://github.com/linfengxu/OmniNeo
    """
}

// parse command line parameters
// only set necessary control parameters
params.help = false
params.skip_dna = false
params.skip_rna = false
params.skip_combine_ms = false
params.skip_maxquant = false
params.skip_netmhcpan = false

// display help information
if (params.help) {
    helpMessage()
    exit 0
}

// check necessary parameters
if (!params.samplesheet) {
    log.error "Error: Sample sheet file not specified, please use the --samplesheet option"
    exit 1
}

if (!params.hg38db) {
    log.error "Error: Reference genome path not specified, please use the --hg38db option"
    exit 1
}

// print workflow startup information
log.info """
=============================================================================================
                                OmniNeo PIPELINE v1.0
=============================================================================================
🎉 Pipeline Execution Started! 🎉

Sample Sheet: ${params.samplesheet}
Output Directory: ${params.outdir}
Reference Genome: ${params.hg38db}
Maximum Memory: ${params.max_memory}
Maximum CPUs: ${params.max_cpus}
Maximum Run Time: ${params.max_time}
Skip DNA Analysis: ${params.skip_dna}
Skip RNA Analysis: ${params.skip_rna}
Skip MS Database Creation: ${params.skip_combine_ms}
Skip MaxQuant Analysis: ${params.skip_maxquant}
Skip netMHCpan Analysis: ${params.skip_netmhcpan}
"""

// import sub workflows
include { DNA_WORKFLOW } from './workflows/dna_workflow'
include { RNA_WORKFLOW } from './workflows/rna_workflow'
include { DNA_DOWNSTREAM } from './workflows/dna_downstream'
include { RNA_DOWNSTREAM } from './workflows/rna_downstream'
include { COMBINE_MS } from './workflows/Combine_MS'
include { MAXQUANT_ANALYSIS } from './modules/maxquant/maxquant'
include { NETMHCPAN_ANALYSIS } from './workflows/all_netMHCpan'

// define main workflow
workflow {
    // create main input channel
    samples_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true, sep:'\t')
        .map { row -> 
            tuple(
                row.sampleID,
                row.DNANormalReads1,
                row.DNANormalReads2,
                row.DNATumorReads1,
                row.DNATumorReads2,
                row.RNANormalReads1,
                row.RNANormalReads2,
                row.RNATumorReads1,
                row.RNATumorReads2
            )
        }

    // initialize result channels
    dna_results = Channel.empty()
    rna_results = Channel.empty()
    ms_database_results = Channel.empty()

    // 运行 DNA 分析
    if (!params.skip_dna) {
        dna_results = DNA_WORKFLOW(
            samples_ch,
            params.hg38db
        )
        
        // 将DNA工作流的输出传递给下游分析工作流
        dna_downstream_results = DNA_DOWNSTREAM(
            dna_results.annotation_results
        )
        
        // 输出整合结果
        dna_downstream_results.merged_fasta.subscribe { sample_id, merged_fasta ->
            log.info "DNA merged peptides for sample: ${sample_id}, file: ${merged_fasta}"
        }
    }

    // 运行 RNA 分析
    if (!params.skip_rna) {
        // 执行RNA工作流
        rna_results = RNA_WORKFLOW(
            samples_ch,
            params.hg38db
        )
        
        // 检查融合结果是否为空
        if (rna_results.fusion_results == null) {
            log.warn "No fusion results detected from RNA workflow, downstream analysis may be affected"
        }
        
        // 将RNA工作流的输出传递给RNA下游分析工作流
        rna_downstream_results = RNA_DOWNSTREAM(
            rna_results.annotation_results,
            rna_results.kallisto_abundance,
            rna_results.fusion_results
        )
        
        // 输出RNA分析结果日志
        rna_downstream_results.all_short_peptides.subscribe { sample_id, peptide_file ->
            log.info "RNA short peptides for sample: ${sample_id}, file: ${peptide_file}"
        }
        
        rna_downstream_results.all_long_peptides.subscribe { sample_id, peptide_file ->
            log.info "RNA long peptides for sample: ${sample_id}, file: ${peptide_file}"
        }
    }
    
    // 合并DNA和RNA肽段用于质谱分析
    if (!params.skip_combine_ms) {
        // 准备参考数据库路径
        uniprot_db = file(params.combine_ms.uniprot_fasta)
        crap_db = file(params.combine_ms.crap_fasta)
        
        // 整合质谱数据库
        ms_database_results = COMBINE_MS(
            dna_downstream_results.merged_fasta,
            rna_downstream_results.all_short_peptides,
            rna_downstream_results.all_long_peptides,
            uniprot_db,
            crap_db
        )
        
        // 运行MaxQuant分析（如果未跳过）
        if (!params.skip_maxquant) {
            // 从samplesheet中获取MS文件路径
            ms_files_ch = Channel.fromPath(params.samplesheet)
                .splitCsv(header:true, sep:'\t')
                .map { row -> tuple(row.sampleID, row.MSTumor.split(',')) }
            
            // 运行MaxQuant分析
            maxquant_results = MAXQUANT_ANALYSIS(
                ms_database_results.clean_ms_db,
                ms_files_ch
            )
        }
    }
    
    // 运行netMHCpan分析
    if (!params.skip_netmhcpan) {


        // 运行netMHCpan分析工作流
        netmhcpan_results = NETMHCPAN_ANALYSIS(
            dna_downstream_results.merged_fasta,
            rna_downstream_results.all_short_peptides,
            rna_downstream_results.all_long_peptides,
            rna_results.hla_results,
            rna_results.hla2_results
        )
        
        // 输出netMHCpan分析结果
        netmhcpan_results.class1_8mer_results.subscribe { sample_id, results ->
            log.info "HLA-I 8mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class1_9mer_results.subscribe { sample_id, results ->
            log.info "HLA-I 9mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class1_10mer_results.subscribe { sample_id, results ->
            log.info "HLA-I 10mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class1_11mer_results.subscribe { sample_id, results ->
            log.info "HLA-I 11mer predictions completed for sample: ${sample_id}"
        }
        
        // HLA-II predictions (12-25mer)
        netmhcpan_results.class2_12mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 12mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_13mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 13mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_14mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 14mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_15mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 15mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_16mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 16mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_17mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 17mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_18mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 18mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_19mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 19mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_20mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 20mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_21mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 21mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_22mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 22mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_23mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 23mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_24mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 24mer predictions completed for sample: ${sample_id}"
        }
        netmhcpan_results.class2_25mer_results.subscribe { sample_id, results ->
            log.info "HLA-II 25mer predictions completed for sample: ${sample_id}"
        }
        
        netmhcpan_results.netctl_results.subscribe { sample_id, results ->
            log.info "netCTL predictions completed for sample: ${sample_id}"
        }
    }
}

// 工作流完成消息
workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Execution duration: ${workflow.duration}"
    
    // 添加跳过的进程统计
    def skipped = workflow.stats.succeededCount - workflow.stats.cachedCount
    log.info """
    Pipeline statistics:
    ------------------------------------------------------------------------
    Completed processes : ${workflow.stats.succeededCount}
    Cached (skipped)    : ${workflow.stats.cachedCount}
    Actually run       : ${skipped}
    Failed processes   : ${workflow.stats.failedCount}
    ------------------------------------------------------------------------
    """
}


