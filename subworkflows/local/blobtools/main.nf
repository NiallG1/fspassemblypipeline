include { BLOBTK_PLOT               } from '../../../modules/nf-core/blobtk/plot/main'
include { BLOBTOOLKIT_CREATEBLOBDIR } from '../../../modules/local/blobtoolkit_create/main'
//include { SAMTOOLS_INDEX            } from '../../../modules/local/samtools/index/main'

workflow BLOBTOOLS {

    take:
    ch_samplesheet
    


    main:
    ch_versions = channel.empty()


    //
    // Create Blobtools dataset files
    //
    BLOBTOOLKIT_CREATEBLOBDIR ( windowstats, busco, blastp, config, taxdump )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_CREATEBLOBDIR.out.versions.first() )


    //
    // Update Blobtools dataset files
    //
    BLOBTOOLKIT_UPDATEBLOBDIR ( BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir, syn_tsv, cat_tsv, blastx, blastn, taxdump )
    ch_versions = ch_versions.mix ( BLOBTOOLKIT_UPDATEBLOBDIR.out.versions.first() )


    emit:
    blobdir  = BLOBTOOLKIT_UPDATEBLOBDIR.out.blobdir  // channel: [ val(meta), path(dir) ]
    versions = ch_versions                            // channel: [ versions.yml ]
}
