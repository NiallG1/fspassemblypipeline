include { SEQKIT_STATS                         } from '../../../modules/nf-core/seqkit/stats/main'
include { GETSEQKITK                            } from '../../../modules/local/getseqkitk/main'
// include { SEQKIT_STATS as SEQKIT_STATS_MERGED  } from '../../../modules/nf-core/seqkit/stats/main' I can't work on this if the preprocessing subworkflow is not updated to ouput merged reads.
include { KMERGENIE                            } from '../../../modules/nf-core/kmergenie/main'
include { GETKMERGENIEK                        } from '../../../modules/local/getkmergeniek/main'
include { FASTK_FASTK                          } from '../../../modules/nf-core/fastk/fastk/main'
include { SPADES as SPADES_MANUAL              } from '../../../modules/nf-core/spades/main'
include { SPADES as SPADES_KMERGENIE           } from '../../../modules/nf-core/spades/main'
include { MEGAHIT as MEGAHIT_MANUAL            } from '../../../modules/nf-core/megahit/main'
include { MEGAHIT as MEGAHIT_KMERGENIE         } from '../../../modules/nf-core/megahit/main'
include { MINIA as MINIA_MANUAL                } from '../../../modules/nf-core/minia/main'
include { MINIA as MINIA_KMERGENIE             } from '../../../modules/nf-core/minia/main'
include { ABYSS_ABYSSPE as ABYSS_MANUAL } from '../../../modules/nf-core/abyss/abysspe/main'
include { ABYSS_ABYSSPE as ABYSS_KMERGENIE } from '../../../modules/nf-core/abyss/abysspe/main'
include { SPARSEASSEMBLER as SPARSEASSEMBLER_MANUAL } from '../../../modules/local/sparseassembler/main'
include { SPARSEASSEMBLER as SPARSEASSEMBLER_KMERGENIE } from '../../../modules/local/sparseassembler/main'
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
    GETSEQKITK   ( SEQKIT_STATS.out.stats )
//    SEQKIT_STATS_MERGED
    FASTK_FASTK  ( ch_fastp_reads )

    KMERGENIE    ( ch_fastp_reads )

    GETKMERGENIEK ( KMERGENIE.out.html )

    // Channel 1: Manual strategy (uses config values)
    ch_reads_manual_strategy = ch_fastp_reads
        .map { meta, reads ->
            [meta + [kmer_strategy: 'manual'], reads]
        }

    // Channel 2: KmerGenie strategy (adds predicted kmer to default list)
    ch_reads_kmergenie_strategy = ch_fastp_reads
        .map { meta, reads -> [meta.id, meta, reads] }
        .join(GETKMERGENIEK.out.kmer_txt.map { meta, kmer_file ->
            [meta.id, kmer_file.text.trim() as Integer]
        })
        .map { id, meta, reads, kmergenie_kmer ->
            // Build k-mer list: add predicted kmer to defaults if valid
            def default_kmers_spades = [21, 33, 55, 77]  // spades recommended defaults
            def max_kmer_spades = 127  // spades max kmer limit

            def default_kmers_megahit = [21, 29, 39, 59, 79, 99, 119, 141] // megahit recommended defaults
            def max_kmer_megahit = 141 // megahit max kmer limit
            def max_step_megahit = 28  // Maximum allowed gap between consecutive kmers in megahit

            // Validate predicted kmer for spades: must be odd, in range [15, 127], and not already in list
            def is_valid_spades = (kmergenie_kmer >= 15 &&
                            kmergenie_kmer <= max_kmer_spades &&
                            kmergenie_kmer % 2 == 1 &&
                            !(kmergenie_kmer in default_kmers_spades))

            // Validate predicted kmer for megahit: must be odd, in range [15, 141], and not already in list
            def is_valid_megahit = (kmergenie_kmer >= 15 &&
                            kmergenie_kmer <= max_kmer_megahit &&
                            kmergenie_kmer % 2 == 1 &&
                            !(kmergenie_kmer in default_kmers_megahit))

            // Build final k-mer lists for spades
            def kmer_list_spades = is_valid_spades ?
                (default_kmers_spades + [kmergenie_kmer]).sort() :
                default_kmers_spades

            // Build final k-mer lists for megahit
            def kmer_list_megahit = is_valid_megahit ?
                (default_kmers_megahit + [kmergenie_kmer]).sort() :
                default_kmers_megahit

            // For assemblers that only take a single k-mer, use the predicted k-mer if it's valid (between 15 and 127), or fall back to a fixed value (25)
            def single_kmer = (kmergenie_kmer >= 15 && kmergenie_kmer <= 127 && kmergenie_kmer % 2 == 1) ? kmergenie_kmer : 25

            // Enrich metadata with k-mer strategy and lists
            def enriched_meta = meta + [
                kmer_strategy: 'kmergenie',
                predicted_kmer: kmergenie_kmer,
                kmer_list_spades: kmer_list_spades.join(','),
                kmer_list_megahit: kmer_list_megahit.join(','),
                single_kmer: single_kmer
            ]

            log.info """
            Sample: ${id}
            Predicted kmer: ${kmergenie_kmer}
            SPAdes - Valid: ${is_valid_spades}, List: ${kmer_list_spades.join(',')}
            MEGAHIT - Valid: ${is_valid_megahit}, List: ${kmer_list_megahit.join(',')}
            Single k-mer: ${single_kmer}
            """.stripIndent()

            [enriched_meta, reads]
        }

    // Spades needs a tuple with 4 elements as inputs, so we need to map the channel to add empty lists for the other 2 inputs (see PREPROCESSING subworkflow for example)
    // SPADES: [ meta, illumina, pacbio, nanopore ]
    // ch_input_reads_spades = ch_fastp_reads.map { meta, reads -> [ meta, reads, [], [] ] }

    SPADES_MANUAL       ( ch_reads_manual_strategy.map { meta, reads -> [meta, reads, [], []] },
    [],
    []
    )

    SPADES_KMERGENIE       ( ch_reads_kmergenie_strategy.map { meta, reads -> [meta, reads, [], []] },
    [],
    []
    )

    // Megahit needs a tuple with 3 elements as input. I can't use ch_fastp_reads directly because R1 and R2 paths there are in a single list element. So I need to map the channel to split R1 and R2 into separate list elements.
    // MEGAHIT: [ meta, reads1, reads2 ]
    ch_input_reads_megahit_manual = ch_reads_manual_strategy.map { meta, reads -> [ meta, reads[0], reads[1] ] }
    MEGAHIT_MANUAL      ( ch_input_reads_megahit_manual )

    ch_input_reads_megahit_kmergenie = ch_reads_kmergenie_strategy.map { meta, reads -> [ meta, reads[0], reads[1] ] }
    MEGAHIT_KMERGENIE ( ch_input_reads_megahit_kmergenie )

    MINIA_MANUAL        ( ch_reads_manual_strategy )
    MINIA_KMERGENIE     ( ch_reads_kmergenie_strategy )

    ch_abyss_input_manual = ch_reads_manual_strategy.map { meta, reads -> [ meta, reads, [] ] }
    ABYSS_MANUAL ( ch_abyss_input_manual, params.abyss_kmer )

    // create an input channel with just the kmer value for abyss and sparseassembler, as they require it as input (can't be passed from the extra args)
    ch_kmergenie_single_kmer = ch_reads_kmergenie_strategy.map { meta, reads -> meta.single_kmer }

    ch_abyss_input_kmergenie = ch_reads_kmergenie_strategy.map { meta, reads -> [ meta, reads, [] ] }
    ABYSS_KMERGENIE ( ch_abyss_input_kmergenie, ch_kmergenie_single_kmer )

    SPARSEASSEMBLER_MANUAL ( ch_reads_manual_strategy, params.sparseassembler_kmer, params.sparseassembler_genome_size, params.sparseassembler_expected_coverage )
    SPARSEASSEMBLER_KMERGENIE ( ch_reads_kmergenie_strategy, ch_kmergenie_single_kmer, params.sparseassembler_genome_size, params.sparseassembler_expected_coverage )

    // input channel for renaming the assemblies. I need to change the meta.id to include the assembler and avoid conflicts in the output names.
    // def ch_draft_assemblies_input = SPADES_MANUAL.out.scaffolds

    def createAssemblyMeta = { meta, assembly, assembler ->
        def strategy = meta.kmer_strategy
        def new_meta = meta + [
            assembly_id: "${meta.id}_${assembler}_${strategy}",
            assembler: assembler,
            id: meta.id
        ]
        return [new_meta, assembly, "${meta.id}_${strategy}_${assembler}.fa"]
    }

    def ch_draft_assemblies_input = SPADES_MANUAL.out.scaffolds
        .mix(SPADES_KMERGENIE.out.scaffolds)
        .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'spades') }
        .mix( MEGAHIT_MANUAL.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'megahit') } )
        .mix( MEGAHIT_KMERGENIE.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'megahit') } )
        .mix( MINIA_MANUAL.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'minia') } )
        .mix( MINIA_KMERGENIE.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'minia') } )
        .mix( ABYSS_MANUAL.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'abyss') } )
        .mix( ABYSS_KMERGENIE.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'abyss') } )
        .mix( SPARSEASSEMBLER_MANUAL.out.scaffolds
            .concat(SPARSEASSEMBLER_MANUAL.out.contigs)
            .unique { meta, assembly -> meta.id }
            .map { meta, assembly -> createAssemblyMeta(meta, assembly, 'sparseassembler') }
        )
        .mix( SPARSEASSEMBLER_KMERGENIE.out.scaffolds
            .concat(SPARSEASSEMBLER_KMERGENIE.out.contigs)
            .unique { meta, assembly -> meta.id }
            .map { meta, assembly -> createAssemblyMeta(meta, assembly, 'sparseassembler') }
        )

    RENAME_ASSEMBLIES ( ch_draft_assemblies_input )

    BUSCO_BUSCO ( RENAME_ASSEMBLIES.out.renamed_assemblies, params.busco_mode, params.busco_lineage, params.busco_lineages_path ?:[], params.busco_config_file ?:[], params.busco_clean_intermediates )

    // Map samples to their best BUSCO database based on taxonomy
    // BUSCO_SPECIFIC runs BUSCO with a specific lineage for each sample. The lineage is determined based on the sample metadata and the list of available lineages.
    // we attempt to find the most specific available BUSCO database for each sample by walking down the taxonomic hierarchy (family → order → class → phylum).
    // If no matching database is found at any level, we fall back to the lineage specified by params.busco_lineage.

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
    getseqkitk_kmer                    = GETSEQKITK.out.seqkitkmer_txt              // channel: [ val(meta), path('*.txt') ]
    kmergenie_html                     = KMERGENIE.out.html             // channel: [ val(meta), path('*.html') ]
    getkmergeniek_k                    = GETKMERGENIEK.out.kmer_txt              // channel: [ val(meta), path('*.k') ]
    fastk_ktab                         = FASTK_FASTK.out.ktab             // channel: [ val(meta), path('*.ktab') ]
    fastk_hist                         = FASTK_FASTK.out.hist             // channel: [ val(meta), path('*.hist') ]
    spades_scaffolds_manual            = SPADES_MANUAL.out.scaffolds             // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    spades_scaffolds_kmergenie         = SPADES_KMERGENIE.out.scaffolds          // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    megahit_contigs_manual             = MEGAHIT_MANUAL.out.contigs              // channel: [ val(meta), path('*.contigs.fa.gz') ]
    megahit_contigs_kmergenie          = MEGAHIT_KMERGENIE.out.contigs           // channel: [ val(meta), path('*.contigs.fa.gz') ]
    minia_contigs_manual               = MINIA_MANUAL.out.contigs                // channel: [ val(meta), path('*.contigs.fa') ]
    minia_contigs_kmergenie            = MINIA_KMERGENIE.out.contigs                // channel: [ val(meta), path('*.contigs.fa') ]
    abyss_scaffolds_manual             = ABYSS_MANUAL.out.scaffolds      // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    abyss_scaffolds_kmergenie          = ABYSS_KMERGENIE.out.scaffolds   // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    sparseassembler_scaffolds_manual          = SPARSEASSEMBLER_MANUAL.out.scaffolds    // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    sparseassembler_contigs_manual            = SPARSEASSEMBLER_MANUAL.out.contigs      // channel: [ val(meta), path('*.contigs.fa.gz') ]
    sparseassembler_scaffolds_kmergenie         = SPARSEASSEMBLER_KMERGENIE.out.scaffolds // channel: [ val(meta), path('*.scaffolds.fa.gz') ]
    sparseassembler_contigs_kmergenie           = SPARSEASSEMBLER_KMERGENIE.out.contigs   // channel: [ val(meta), path('*.contigs.fa.gz') ]
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
