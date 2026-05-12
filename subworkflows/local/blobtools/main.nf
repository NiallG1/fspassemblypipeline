include { BLOBTK_PLOT               } from '../../../modules/nf-core/blobtk/plot/main'
include { BLOBTOOLKIT_CREATEBLOBDIR } from '../../../modules/local/blobtoolkit_create/main'
include { CREATE_PROJECT_YAML       } from '../../../modules/local/createyml/main'


workflow BLOBTOOLS {

    take:
    ch_samplesheet 
   
    main:

        // Create a channel from the BUSCO file path in config
        ch_busco = Channel.fromPath(params.busco_file)
    
        // Create YAML files from samplesheet metadata
        ch_yaml_input = ch_samplesheet
            .map { meta, files ->
                tuple(meta, files[0])  // meta and fasta file
            }
        
        CREATE_PROJECT_YAML(ch_yaml_input)
        
        // Debug: Check what YAML output looks like
        CREATE_PROJECT_YAML.out.yaml.view { meta, yaml -> 
            "YAML created: ${meta.id} -> ${yaml}"
        }

        // Prepare channels for joining
        ch_samplesheet_keyed = ch_samplesheet
            .map { meta, files ->
                tuple(meta.id, meta, files[0], files[1])
            }
            .view { "Samplesheet keyed: ${it}" }
        
        ch_yaml_keyed = CREATE_PROJECT_YAML.out.yaml
            .map { meta, yaml -> 
                tuple(meta.id, yaml)
            }
            .view { "YAML keyed: ${it}" }

        // Join and combine
        ch_btk = ch_samplesheet_keyed
            .join(ch_yaml_keyed)
            .view { "After join: ${it}" }
            .combine(ch_busco)
            .view { "After combine: ${it}" }
            .map { id, meta, fasta, bam, yaml, busco ->
                tuple(meta, fasta, bam, yaml, busco)
            }
            .view { "Final input to BLOBTOOLKIT: ${it}" }

        //
        // Create Blobtools dataset files
        //
        BLOBTOOLKIT_CREATEBLOBDIR(ch_btk)

    emit:
    blobdir  = BLOBTOOLKIT_CREATEBLOBDIR.out.blobdir
}