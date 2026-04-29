include { BLOBTK_PLOT               } from '../../../modules/nf-core/blobtk/plot/main'
include { BLOBTOOLKIT_CREATEBLOBDIR } from '../../../modules/local/blobtoolkit_create/main'


workflow BLOBTOOLS {

    take:
    ch_samplesheet 
   
    main:

        // Create a channel from the BUSCO file path in config
        ch_busco = Channel.fromPath(params.busco_file)
    
        ch_btk = ch_samplesheet
        .map { meta, files ->
            tuple(
                [id: meta.id],      // Create new meta with just id
                files[0],         // First file in the list is the fasta
                files[1]            // Second file in the list is the bam

            )
        }
        .combine(ch_busco)

    //
    // Create Blobtools dataset files
    //
    BLOBTOOLKIT_CREATEBLOBDIR(ch_btk)
    //ch_versions = ch_versions.mix ( BLOBTOOLKIT_CREATEBLOBDIR.out.versions.first() )-think this version control is outdated



    emit:
    blobdir  = BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir
    //versions = ch_versions                            // channel: [ versions.yml ]
}
