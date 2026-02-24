// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

//include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
//include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include   { TIARA_TIARA        } from '../../../modules/nf-core/tiara/tiara/main'
include   { FCSGX_RUNGX        } from '../../../modules/nf-core/fcsgx/rungx/main' 

workflow CONTAMINATION_DETECTION {

    take:
    // TODO nf-core: edit input (take) channels
    ch_tiara_input // [val(meta), path(assemblies)]
    ch_ramdisk_path
    ch_db_path // channel: val(gxdb)
   

    main:
    TIARA_TIARA(ch_tiara_input)  // call module

     ch_fcs_gx = ch_tiara_input.map {meta, assembly -> [meta, params.taxid , assembly]}

    

    // Current (incorrect):
    FCSGX_RUNGX(ch_fcs_gx, params.db_path, params.ramdisk_path ?:[])

    // Fixed - pass the directory containing the database:
   // ch_gxdb_dir = file(params.db_path).parent  // Gets /home/nga10kg/FSP/pipeline/gx_test_db/test-only
   // ch_gxdb_name = file(params.db_path).name   // Gets test-only

    //FCSGX_RUNGX(ch_fcs_gx, ch_gxdb_dir, params.ramdisk_path ?:[])
    
    emit:
    classifications  = TIARA_TIARA.out.classifications  
    versions = TIARA_TIARA.out.versions 
    
    
}




