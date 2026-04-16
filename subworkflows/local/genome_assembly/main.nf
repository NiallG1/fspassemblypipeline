include { SEQKIT_STATS                         } from '../../../modules/nf-core/seqkit/stats/main'
// include { SEQKIT_STATS as SEQKIT_STATS_MERGED  } from '../../../modules/nf-core/seqkit/stats/main' I can't work on this if the preprocessing subworkflow is not updated to ouput merged reads.
include { KMERGENIE                            } from '../../../modules/nf-core/kmergenie/main'
include { GETKMERGENIEK                         } from '../../../modules/local/getkmergeniek/main'
include { FASTK_FASTK                          } from '../../../modules/nf-core/fastk/fastk/main'
include { SPADES                               } from '../../../modules/nf-core/spades/main'
include { MEGAHIT                              } from '../../../modules/nf-core/megahit/main'
include { MINIA                                } from '../../../modules/nf-core/minia/main'
include { ABYSS_ABYSSPE                        } from '../../../modules/nf-core/abyss/abysspe/main'
include { SPARSEASSEMBLER                      } from '../../../modules/local/sparseassembler/main'
include { RENAME_ASSEMBLIES                    } from '../../../modules/local/rename_assemblies/main'
include { BUSCO_BUSCO                          } from '../../../modules/nf-core/busco/busco/main'
include { BUSCO_BUSCO as BUSCO_SPECIFIC        } from '../../../modules/nf-core/busco/busco/main'
include { MERQURYFK_MERQURYFK                  } from '../../../modules/nf-core/merquryfk/merquryfk/main'
include { QUAST                                } from '../../../modules/nf-core/quast/main'

workflow GENOME_ASSEMBLY {

    take:
    ch_fastp_reads // channel: [ val(meta), path(reads) ]

    main:

    SEQKIT_STATS ( ch_fastp_reads )
//    SEQKIT_STATS_MERGED
    FASTK_FASTK  ( ch_fastp_reads )

    ch_input_reads_kmergenie = ch_fastp_reads.map { meta, reads -> [ meta, reads[0], reads[1] ] }
    KMERGENIE    ( ch_fastp_reads )

    GETKMERGENIEK ( KMERGENIE.out.html )

    // Spades needs a tuple with 4 elements as inputs, so we need to map the channel to add empty lists for the other 2 inputs (see PREPROCESSING subworkflow for example)
    // SPADES: [ meta, illumina, pacbio, nanopore ]
    ch_input_reads_spades = ch_fastp_reads.map { meta, reads -> [ meta, reads, [], [] ] }

    SPADES       ( ch_input_reads_spades,
    [],
    []
    )

    // Megahit needs a tuple with 3 elements as input. I can't use ch_fastp_reads directly because R1 and R2 paths there are in a single list element. So I need to map the channel to split R1 and R2 into separate list elements.
    // MEGAHIT: [ meta, reads1, reads2 ]
    ch_input_reads_megahit = ch_fastp_reads.map { meta, reads -> [ meta, reads[0], reads[1] ] }

    MEGAHIT      ( ch_input_reads_megahit )

    MINIA        ( ch_fastp_reads )

    ch_abyss_input = ch_fastp_reads.map { meta, reads -> [ meta, reads, [] ] }
    ABYSS_ABYSSPE ( ch_abyss_input, params.abyss_kmer )

    SPARSEASSEMBLER ( ch_fastp_reads, params.sparseassembler_kmer, params.sparseassembler_genome_size, params.sparseassembler_expected_coverage )

    // input channel for renaming the assemblies. I need to change the meta.id to include the assembler and avoid conflicts in the output names.
    def ch_draft_assemblies_input = SPADES.out.scaffolds.map { meta, scaffolds ->
        // add assembler name to meta.id to ensure unique output names
        def assembler = 'spades'
        def new_meta = meta + [assembly_id: "${meta.id}_${assembler}", assembler: 'spades', id: meta.id]
        return [ new_meta, scaffolds, "${meta.id}_${assembler}.fa" ]
    }
    .mix( MEGAHIT.out.contigs.map { meta, contigs ->
        def assembler = 'megahit'
        def new_meta = meta + [assembly_id: "${meta.id}_${assembler}", assembler: 'megahit', id: meta.id]
        return [ new_meta, contigs, "${meta.id}_${assembler}.fa" ]
    } )
    .mix( MINIA.out.contigs.map { meta, contigs ->
        def assembler = 'minia'
        def new_meta = meta + [assembly_id: "${meta.id}_${assembler}", assembler: 'minia', id: meta.id]
        return [ new_meta, contigs, "${meta.id}_${assembler}.fa" ]
    } )
    .mix( ABYSS_ABYSSPE.out.contigs.map { meta, contigs ->
        def assembler = 'abyss'
        def new_meta = meta + [assembly_id: "${meta.id}_${assembler}", assembler: 'abyss', id: meta.id]
        return [ new_meta, contigs, "${meta.id}_${assembler}.fa" ]
    } )
    .mix( SPARSEASSEMBLER.out.scaffolds
    .concat(SPARSEASSEMBLER.out.contigs)
    .unique { meta, assembly -> meta.id }
    .map { meta, assembly ->
        def assembler = 'sparseassembler'
        def new_meta = meta + [assembly_id: "${meta.id}_${assembler}", assembler: 'sparseassembler', id: meta.id]
        return [ new_meta, assembly, "${meta.id}_${assembler}.fa" ]
    } )

    RENAME_ASSEMBLIES ( ch_draft_assemblies_input )

    BUSCO_BUSCO ( RENAME_ASSEMBLIES.out.renamed_assemblies, params.busco_mode, params.busco_lineage, params.busco_lineages_path ?:[], params.busco_config_file ?:[], params.busco_clean_intermediates )

    // Map samples to their best BUSCO database based on taxonomy
    // Load BUSCO lineages once
    def busco_lineages_file = file(params.lineages_list_file)

    def busco_lineages = busco_lineages_file
        .readLines() as Set

    def ext = "_${params.busco_db_extension}"
    def fallback = params.busco_lineage

    // Add lineage to each sample and split into synchronized channels
    ch_assemblies_with_lineage = RENAME_ASSEMBLIES.out.renamed_assemblies
        .map { meta, assembly ->
            // Find best matching lineage for this sample
            def lineage = ['family', 'order', 'class', 'phylum']
                .findResult { level ->
                    def taxon = meta[level]
                    if (!taxon) {
                        log.info "  ${level}: not set"
                        return null
                    }
                    def candidate = "${taxon.toLowerCase()}${ext}".toString()
                    def found = candidate in busco_lineages
                    found ? candidate : null
                } ?: fallback

            // Return tuple with lineage for splitting
            tuple(meta, assembly, lineage)
        }
        .multiMap { meta, assembly, lineage ->
            assemblies: tuple(meta, assembly)
            lineages: lineage
        }

    // Specific BUSCO with per-sample lineage
    BUSCO_SPECIFIC(
        ch_assemblies_with_lineage.assemblies,
        params.busco_mode,
        ch_assemblies_with_lineage.lineages,
        params.busco_lineages_path ?: [],
        params.busco_config_file ?: [],
        params.busco_clean_intermediates
    )



    // input channel for the first input required by merquryfk: tuple val(meta) , path(fastk_hist), path(fastk_ktab), path(assembly), path(haplotigs)
    // to obtain this:
    // 1. join fastk hist and ktab in a single list and map to meta.id to be use as key to then join with the assemblies
    def ch_combined_fastk = FASTK_FASTK.out.hist.join(FASTK_FASTK.out.ktab, by: 0).map { meta, hist, ktab -> [ meta.id, hist, ktab ] }
    // 2. map renamed assemblies to original meta.id (to be used as key to then join with combined fastk results)
    def ch_draft_assemblies_mapped_to_id = RENAME_ASSEMBLIES.out.renamed_assemblies.map { meta, renamed_assembly -> [ meta.id, meta, renamed_assembly ]}
    // 3. join combined fastk with renamed assemblies using meta.id as key
    def ch_merquryfk_input = ch_combined_fastk.combine( ch_draft_assemblies_mapped_to_id, by: 0 ).map { sample_id, hist, ktab, meta, assembly -> [ meta, hist, ktab, assembly, [] ] }

    MERQURYFK_MERQURYFK ( ch_merquryfk_input, [[],[]], [[],[]] ) // no mathernal and pathernal haplotypes for trio mode

    // input channel for quast: I want to run quast once per sample, so I have to group the different assemblies per sample name
    // ch_draft_assemblies_mapped_to_id is: [ sample_id, meta, assembly ]
    // groupTuple(by: 0) groups by position 0 (sample_id)
    // Result: [ sample_id, [meta1, meta2, meta3], [assembly1, assembly2, assembly3] ]
    def ch_quast_input = ch_draft_assemblies_mapped_to_id.groupTuple( by: 0 ).map { sample_id, metas, assemblies -> [ [id: sample_id], assemblies ]}

    QUAST ( ch_quast_input,[[],[]], [[],[]] ) // no reference fasta or gff for quast

    emit:
    seqkit_stats                       = SEQKIT_STATS.out.stats           // channel: [ val(meta), [ bam ] ]
    fastk_ktab                         = FASTK_FASTK.out.ktab             // channel: [ val(meta), path('*.ktab') ]
    fastk_hist                         = FASTK_FASTK.out.hist             // channel: [ val(meta), path('*.hist') ]
    spades_scaffolds                   = SPADES.out.scaffolds             // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    megahit_contigs                    = MEGAHIT.out.contigs              // channel: [ val(meta), path('*.contigs.fa.gz') ]
    minia_contigs                      = MINIA.out.contigs                // channel: [ val(meta), path('*.contigs.fa') ]
    abyss_scaffolds                    = ABYSS_ABYSSPE.out.scaffolds      // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    sparseassembler_scaffolds          = SPARSEASSEMBLER.out.scaffolds    // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    renamed_assemblies                 = RENAME_ASSEMBLIES.out.renamed_assemblies // channel: [ val(meta), path('*.fa.gz') ]
    busco_batch_summary                = BUSCO_BUSCO.out.batch_summary  // channel: [ val(meta), path('*.busco.batch_summary.txt') ]
    busco_short_summaries_txt          = BUSCO_BUSCO.out.short_summaries_txt  // channel: [ val(meta), path('short_summary.*.txt') ]
    busco_full_table                   = BUSCO_BUSCO.out.full_table  // channel: [ val(meta), path('full_table.*.txt') ]
    busco_batch_summary_specific       = BUSCO_SPECIFIC.out.batch_summary  // channel: [ val(meta), path('*.busco.batch_summary.txt') ]
    busco_short_summaries_txt_specific = BUSCO_SPECIFIC.out.short_summaries_txt  // channel: [ val(meta), path('short_summary.*.txt') ]
    busco_full_table_specific          = BUSCO_SPECIFIC.out.full_table  // channel: [ val(meta), path('full_table.*.txt') ]
    merquryfk_completeness_stats       = MERQURYFK_MERQURYFK.out.stats // channel: [ val(meta), path('*.completeness.stats') ]
    quast_results                      = QUAST.out.results         // channel: [ val(meta), path("${prefix}") ]
}
