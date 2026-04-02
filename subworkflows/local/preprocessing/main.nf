include { FALCO as FALCO_RAW            } from '../../../modules/nf-core/falco/main'
include { FALCO as FALCO_AFTER_FASTP    } from '../../../modules/nf-core/falco/main'
include { FALCO as FALCO_AFTER_MERGE    } from '../../../modules/nf-core/falco/main'
include { FASTP as FASTP_TRIM           } from '../../../modules/nf-core/fastp/main'
include { FASTP as FASTP_MERGE          } from '../../../modules/nf-core/fastp/main'
include { FASTK_FASTK                   } from '../../../modules/nf-core/fastk/fastk/main'
include { FASTK_HISTEX                  } from '../../../modules/nf-core/fastk/histex/main'
include { GENESCOPEFK as GENESCOPEFK_P1 } from '../../../modules/nf-core/genescopefk/main'
include { GENESCOPEFK as GENESCOPEFK_P2 } from '../../../modules/nf-core/genescopefk/main'

workflow PREPROCESSING {

    take:

    ch_samplesheet // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    FALCO_RAW (
        ch_samplesheet
    )
    ch_versions = ch_versions.mix( FALCO_RAW.out.versions_falco )

    def fastp_adapter_fasta = params.fastp_adapter_fasta ? file(params.fastp_adapter_fasta, checkIfExists: true) : []

    ch_samplesheet_fastp_trim = ch_samplesheet.map { meta, reads -> [ meta, reads, fastp_adapter_fasta ] }

    FASTP_TRIM (
        ch_samplesheet_fastp_trim,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix( FASTP_TRIM.out.versions_fastp )

    ch_samplesheet_fastp_merge = FASTP_TRIM.out.reads.map { meta, reads -> [ meta, reads, [] ] }

    FASTP_MERGE (
        ch_samplesheet_fastp_merge,
        false,
        false,
        true
    )
    ch_versions = ch_versions.mix( FASTP_MERGE.out.versions_fastp )

    FALCO_AFTER_FASTP (
        FASTP_TRIM.out.reads
    )
    ch_versions = ch_versions.mix( FALCO_AFTER_FASTP.out.versions_falco )

    FALCO_AFTER_MERGE (
        FASTP_MERGE.out.reads_merged
    )
    ch_versions = ch_versions.mix( FALCO_AFTER_MERGE.out.versions_falco )

    FASTK_FASTK (
        FASTP_TRIM.out.reads
    )
    ch_versions = ch_versions.mix( FASTK_FASTK.out.versions_fastk )

    FASTK_HISTEX (
        FASTK_FASTK.out.hist
    )
    ch_versions = ch_versions.mix( FASTK_HISTEX.out.versions_fastk )

    GENESCOPEFK_P1 (
        FASTK_HISTEX.out.hist
    )
    ch_versions = ch_versions.mix( GENESCOPEFK_P1.out.versions_genescopefk )

    GENESCOPEFK_P2 (
        FASTK_HISTEX.out.hist
    )
    ch_versions = ch_versions.mix( GENESCOPEFK_P2.out.versions_genescopefk )

    emit:
    fastp_reads            = FASTP_TRIM.out.reads
    fastp_reads_merged     = FASTP_MERGE.out.reads_merged
    falco_raw_html         = FALCO_RAW.out.html
    falco_after_fastp_html = FALCO_AFTER_FASTP.out.html
    falco_after_merge_html = FALCO_AFTER_MERGE.out.html
    fastk_ktab             = FASTK_FASTK.out.ktab
    fastk_hist             = FASTK_FASTK.out.hist
    histex_txt             = FASTK_HISTEX.out.hist
    genomescope_summary    = GENESCOPEFK_P1.out.summary.mix( GENESCOPEFK_P2.out.summary )
    versions               = ch_versions
}
