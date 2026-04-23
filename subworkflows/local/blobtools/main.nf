include { BLOBTK_PLOT               } from '../../../modules/nf-core/blobtk/plot/main'
include { BLOBTOOLKIT_CREATEBLOBDIR } from '../../../modules/local/blobtoolkit_create/main'
//include { SAMTOOLS_INDEX            } from '../../../modules/local/samtools/index/main'

workflow BLOBTOOLS {

    take:
    ch_samplesheet 
   
    main:
    // im copying the tiara samplesheet parsing block as initally i just want to create
    // a blank blobdir with just a fasta file.
    ch_btk = ch_samplesheet
        .map { meta, files ->
            tuple(
                [id: meta.id],      // Create new meta with just id
                files[0]            // First file in the list is the fasta
            )
        }

    //ch_versions = channel.empty()-think this version control is outdated

    //
    // Create Blobtools dataset files
    //
    BLOBTOOLKIT_CREATEBLOBDIR(ch_btk)
    //ch_versions = ch_versions.mix ( BLOBTOOLKIT_CREATEBLOBDIR.out.versions.first() )-think this version control is outdated


    //
    // Update Blobtools dataset files
    //
    //BLOBTOOLKIT_UPDATEBLOBDIR ( BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir, syn_tsv, cat_tsv, blastx, blastn, taxdump )
    //ch_versions = ch_versions.mix ( BLOBTOOLKIT_UPDATEBLOBDIR.out.versions.first() )


    emit:
    blobdir  = BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir
    //versions = ch_versions                            // channel: [ versions.yml ]
}
