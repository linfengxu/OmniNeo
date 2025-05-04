#!/usr/bin/env nextflow

// 定义帮助信息
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
    For more information and bug reports, please visit: https://github.com/yourusername/pipeline
    """
}

// 解析命令行参数
// 仅设置必要的控制参数
params.help = false
params.skip_dna = false
params.skip_rna = false
params.skip_combine_ms = false
params.skip_maxquant = false

// 显示帮助信息
if (params.help) {
    helpMessage()
    exit 0
}

// 检查必要参数
if (!params.samplesheet) {
    log.error "Error: Sample sheet file not specified, please use the --samplesheet option"
    exit 1
}

if (!params.hg38db) {
    log.error "Error: Reference genome path not specified, please use the --hg38db option"
    exit 1
}

// 打印工作流启动信息
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
"""

// 导入子工作流
include { DNA_WORKFLOW } from './workflows/dna_workflow'
include { RNA_WORKFLOW } from './workflows/rna_workflow'
include { DNA_DOWNSTREAM } from './workflows/dna_downstream'
include { RNA_DOWNSTREAM } from './workflows/rna_downstream'
include { COMBINE_MS } from './workflows/Combine_MS'
include { MAXQUANT_ANALYSIS } from './modules/maxquant/maxquant'

// 定义主工作流
workflow {
    // 创建主输入通道
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

    // 初始化通道变量
    dna_merged_fasta_ch = Channel.empty()
    rna_short_peptides_ch = Channel.empty()
    rna_long_peptides_ch = Channel.empty()

    // 运行 DNA 分析
    if (!params.skip_dna) {
        dna_results = DNA_WORKFLOW(
            samples_ch,
            params.hg38db
        )
        
        // 将DNA工作流的输出传递给下游分析工作流
        dna_downstream_results = DNA_DOWNSTREAM(
            dna_results.annotation_vcf
        )
        
        // 获取DNA合并的肽段文件
        dna_merged_fasta_ch = dna_downstream_results.merged_fasta
        
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
        
        // 获取RNA短肽段和长肽段
        rna_short_peptides_ch = rna_downstream_results.all_short_peptides
        rna_long_peptides_ch = rna_downstream_results.all_long_peptides
        
        // 输出合并后的肽段结果
        rna_downstream_results.all_short_peptides.subscribe { sample_id, peptide_file ->
            log.info "RNA short peptides for sample: ${sample_id}, file: ${peptide_file}"
        }
        
        rna_downstream_results.all_long_peptides.subscribe { sample_id, peptide_file ->
            log.info "RNA long peptides for sample: ${sample_id}, file: ${peptide_file}"
        }
        
        // 输出融合肽段结果的特定日志
        rna_downstream_results.fusion_short_peptides.subscribe { sample_id, peptide_file ->
            log.info "RNA fusion short peptides for sample: ${sample_id}, file: ${peptide_file}"
        }
        
        rna_downstream_results.fusion_long_peptides.subscribe { sample_id, peptide_file ->
            log.info "RNA fusion long peptides for sample: ${sample_id}, file: ${peptide_file}"
        }
    }
    
    // 合并DNA和RNA肽段用于质谱分析
    if (!params.skip_combine_ms) {
        // 准备参考数据库路径
        uniprot_db = file(params.combine_ms.uniprot_fasta)
        crap_db = file(params.combine_ms.crap_fasta)
        
        // 运行合并MS数据库进程
        ms_results = COMBINE_MS(
            dna_merged_fasta_ch,        // DNA合并的肽段
            rna_short_peptides_ch,      // RNA短肽段
            rna_long_peptides_ch,       // RNA长肽段
            uniprot_db,                 // UniProt参考数据库
            crap_db                     // cRAP污染物数据库
        )
        
        ms_results.combined_ms_db.subscribe { sample_id, ms_db ->
            log.info "Mass spectrometry database created for sample: ${sample_id}"
            log.info "MS database location: ${params.outdir}/${sample_id}/ms_database/${ms_db.getName()}"
        }
        
        ms_results.database_stats.subscribe { stats_file ->
            log.info "Mass spectrometry database statistics available in: ${stats_file}"
        }
        
        // 运行MaxQuant分析（如果未跳过）
        if (!params.skip_maxquant) {
            // 运行MaxQuant分析，使用来自samplesheet的MS文件路径
            maxquant_results = MAXQUANT_ANALYSIS(
                ms_results.combined_ms_db,
                params.samplesheet
            )
            
            // 输出MaxQuant配置文件和结果信息
            maxquant_results.mqpar_files.subscribe { sample_id, mqpar_file ->
                log.info "MaxQuant configuration generated for sample: ${sample_id}"
                log.info "Configuration file: ${params.maxquant.mqpar_dir}/mqpar_${sample_id}.xml"
            }
            
            if (!params.maxquant.skip_run) {
                maxquant_results.results.subscribe { sample_id, results ->
                    log.info "MaxQuant analysis completed for sample: ${sample_id}"
                    log.info "Results available in: ${params.maxquant.results_dir}/${sample_id}"
                }
            }
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


