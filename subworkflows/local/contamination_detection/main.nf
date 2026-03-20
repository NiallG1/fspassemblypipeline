
//include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
//include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include   { TIARA_TIARA        } from '../../../modules/nf-core/tiara/tiara/main'
include   { FCSGX_RUNGX        } from '../../../modules/nf-core/fcsgx/rungx/main' 
include   { CONVERTFCSRPT      } from '../../../modules/local/convertfcsrpt/main'
include   { COMPARISON         } from '../../../modules/local/comparison/main'

workflow CONTAMINATION_DETECTION {

    take:
    // TODO nf-core: edit input (take) channels
    ch_tiara_input // [val(meta), path(assemblies)]
    ch_ramdisk_path
    ch_db_path // channel: val(gxdb)
   

    main:
    TIARA_TIARA(ch_tiara_input)  // call module

    ch_fcs_gx = ch_tiara_input.map {meta, assembly -> [meta, params.taxid , assembly]}

    

    FCSGX_RUNGX(ch_fcs_gx, params.db_path, params.ramdisk_path ?:[])



    CONVERTFCSRPT(FCSGX_RUNGX.out.taxonomy_report)

    COMPARISON(
    TIARA_TIARA.out.classifications,
    CONVERTFCSRPT.out.fcs_report_reformatted
)
    
    emit:
    tiara_classifications  = TIARA_TIARA.out.classifications
    taxonomy_report        = FCSGX_RUNGX.out.taxonomy_report   
    fcsgx_reformatted      = CONVERTFCSRPT.out.fcs_report_reformatted  
    
    // versions
    versions = TIARA_TIARA.out.versions.mix(FCSGX_RUNGX.out.versions) 


   
    
    
}




