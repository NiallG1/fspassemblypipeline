include { SEQKIT_STATS                                     } from '../../../modules/nf-core/seqkit/stats/main'
include { GETSEQKITK                                       } from '../../../modules/local/getseqkitk/main'
include { KMERGENIE                                        } from '../../../modules/nf-core/kmergenie/main'
include { GETKMERGENIEK                                    } from '../../../modules/local/getkmergeniek/main'
//include { SPADES as SPADES_MANUAL                          } from '../../../modules/nf-core/spades/main'
//include { SPADES as SPADES_KMERGENIE                       } from '../../../modules/nf-core/spades/main'
//include { SPADES as SPADES_READS_LENGTH                    } from '../../../modules/nf-core/spades/main'
include { MEGAHIT as MEGAHIT_MANUAL                        } from '../../../modules/nf-core/megahit/main'
include { MEGAHIT as MEGAHIT_KMERGENIE                     } from '../../../modules/nf-core/megahit/main'
include { MEGAHIT as MEGAHIT_READS_LENGTH                  } from '../../../modules/nf-core/megahit/main'
include { MINIA as MINIA_MANUAL                            } from '../../../modules/nf-core/minia/main'
include { MINIA as MINIA_KMERGENIE                         } from '../../../modules/nf-core/minia/main'
include { MINIA as MINIA_READS_LENGTH                      } from '../../../modules/nf-core/minia/main'
include { ABYSS_ABYSSPE as ABYSS_MANUAL                    } from '../../../modules/nf-core/abyss/abysspe/main'
include { ABYSS_ABYSSPE as ABYSS_KMERGENIE                 } from '../../../modules/nf-core/abyss/abysspe/main'
include { ABYSS_ABYSSPE as ABYSS_READS_LENGTH              } from '../../../modules/nf-core/abyss/abysspe/main'
include { SPARSEASSEMBLER as SPARSEASSEMBLER_MANUAL        } from '../../../modules/local/sparseassembler/main'
include { SPARSEASSEMBLER as SPARSEASSEMBLER_KMERGENIE     } from '../../../modules/local/sparseassembler/main'
include { SPARSEASSEMBLER as SPARSEASSEMBLER_READS_LENGTH  } from '../../../modules/local/sparseassembler/main'

workflow GENOME_ASSEMBLY_MERGED {

    take:
    ch_fastp_reads_merged // channel: [ val(meta), path(reads) ]
    ch_fastp_reads_unmerged // channel: [ val(meta), path(reads) ] - unmerged reads after merge attempt - needed for spades
    ch_fastp_reads // channel: [ val(meta), path(reads) ] - original paired reads, needed by abyss

    main:

    // Add reads_type to the meta
    ch_merged_reads = ch_fastp_reads_merged.map { meta, reads ->
        def new_meta = meta + [ reads_type: 'merged', single_end: true ]
        tuple(new_meta, reads)
    }

// ==================== K-mer strategies for genome assembly =======================

    def ch_reads_manual_strategy = Channel.empty()
    def ch_reads_kmergenie_strategy = Channel.empty()
    def ch_reads_reads_length_strategy = Channel.empty()

    // Channel 1: Manual strategy (uses config values)
    if (!params.skip_manual_strategy) {
        ch_reads_manual_strategy = ch_merged_reads
            .map { meta, reads ->
                [meta + [kmer_strategy: 'manual'], reads]
            }
    }
    // Channel 2: KmerGenie strategy (adds predicted kmer to default list)
    // KMERGENIE and GETKMERGENIEK only run if skip_kmergenie_strategy is false.
    // The channel consequently is only populated if skip_kmergenie_strategy is false.

    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_kmergenie_html = channel.empty()
    def ch_getkmergeniek_k = channel.empty()

    if (!params.skip_kmergenie_strategy) {
        KMERGENIE(ch_merged_reads)
        GETKMERGENIEK(KMERGENIE.out.html)

        // Capture output for emits (avoids using conditionals in emits section)
        ch_kmergenie_html = KMERGENIE.out.html
        ch_getkmergeniek_k = GETKMERGENIEK.out.kmer_txt

        ch_reads_kmergenie_strategy = ch_merged_reads
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
                [enriched_meta, reads]
        }
    }

    // Channel 3: reads_length strategy (adds kmer calculated from median reads length to default list)
    // SEQKIT_STATS, GETSEQKITK only run if skip_reads_length_strategy is false.
    // The channel consequently is only populated if skip_reads_length_strategy is false.


    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_seqkit_stats = channel.empty()
    def ch_getseqkitk_kmer = channel.empty()

    if (!params.skip_reads_length_strategy) {
        SEQKIT_STATS(ch_merged_reads)
        GETSEQKITK(SEQKIT_STATS.out.stats)

        // Capture output for emits (avoids using conditionals in emits section)
        ch_seqkit_stats = SEQKIT_STATS.out.stats
        ch_getseqkitk_kmer = GETSEQKITK.out.seqkitkmer_txt

        ch_reads_reads_length_strategy = ch_merged_reads
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(GETSEQKITK.out.seqkitkmer_txt.map { meta, kmer_file ->
                [meta.id, kmer_file.text.trim() as Integer]
            })
            .map { id, meta, reads, seqkit_kmer ->
                // Build k-mer list: add predicted kmer to defaults if valid
                def default_kmers_spades = [21, 33, 55, 77]  // spades recommended defaults
                def max_kmer_spades = 127  // spades max kmer limit

                def default_kmers_megahit = [21, 29, 39, 59, 79, 99, 119, 141] // megahit recommended defaults
                def max_kmer_megahit = 141 // megahit max kmer limit

                // Validate predicted kmer for spades: must be odd, in range [15, 127], and not already in list
                def is_valid_spades = (seqkit_kmer >= 15 &&
                                seqkit_kmer <= max_kmer_spades &&
                                seqkit_kmer % 2 == 1 &&
                                !(seqkit_kmer in default_kmers_spades))

                // Validate predicted kmer for megahit: must be odd, in range [15, 141], and not already in list
                def is_valid_megahit = (seqkit_kmer >= 15 &&
                                seqkit_kmer <= max_kmer_megahit &&
                                seqkit_kmer % 2 == 1 &&
                                !(seqkit_kmer in default_kmers_megahit))

                // Build final k-mer lists for spades
                def kmer_list_spades = is_valid_spades ?
                    (default_kmers_spades + [seqkit_kmer]).sort() :
                    default_kmers_spades

                // Build final k-mer lists for megahit
                def kmer_list_megahit = is_valid_megahit ?
                    (default_kmers_megahit + [seqkit_kmer]).sort() :
                    default_kmers_megahit

                // For assemblers that only take a single k-mer, use the predicted k-mer if it's valid (between 15 and 127), or fall back to a fixed value (25)
                def single_kmer = (seqkit_kmer >= 15 && seqkit_kmer <= 127 && seqkit_kmer % 2 == 1) ? seqkit_kmer : 25

                // Enrich metadata with k-mer strategy and lists
                def enriched_meta = meta + [
                    kmer_strategy: 'reads_length',
                    predicted_kmer: seqkit_kmer,
                    kmer_list_spades: kmer_list_spades.join(','),
                    kmer_list_megahit: kmer_list_megahit.join(','),
                    single_kmer: single_kmer
                ]
                [enriched_meta, reads]
            }
    }
// =================== End of k-mer strategies for genome assembly =======================


// =================== Genome assembly with different assemblers and k-mer strategies =======================

    // Add assembler name to meta.id for all assemblies to avoid conflicts in downstream processes that use meta.id for naming outputs (e.g. busco, quast, merquryfk)
    def createAssemblyMeta = { meta, assembly, assembler ->
        def strategy = meta.kmer_strategy
        def reads_type = meta.reads_type
        def new_meta = meta + [
            assembly_id: "${meta.id}_${reads_type}_${assembler}_${strategy}",
            assembler: assembler,
            id: meta.id
        ]
        return [new_meta, assembly, "${meta.id}_${reads_type}_${strategy}_${assembler}.fa"]
    }

    // Create channel with new meta for downstream processes (after assembly). This channel combines all assemblies from different assemblers and strategies, and maps them to the new meta with updated id.
    def ch_draft_assemblies_input = Channel.empty()

    // ======= Spades assemblies - nested conditionals (assembler × strategy) ======
    // Spades is only run if skip_spades is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with Spades assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.
    // Spades needs a tuple with 4 elements as inputs, so we need to map the channel to add empty lists for the other 2 inputs (see PREPROCESSING subworkflow for example)
    // SPADES: [ meta, illumina, pacbio, nanopore ]


    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    // def ch_spades_scaffolds_manual = channel.empty()
    // def ch_spades_scaffolds_kmergenie = channel.empty()
    // def ch_spades_scaffolds_reads_length = channel.empty()

    // if (!params.skip_spades) {

    //     if (!params.skip_manual_strategy) {
    //         SPADES_MANUAL(
    //             ch_reads_manual_strategy.map { meta, reads -> [meta, reads, [], []] },
    //             [],
    //             []
    //         )
    //         // Capture output for emits (avoids using conditionals in emits section)
    //         ch_spades_scaffolds_manual = SPADES_MANUAL.out.scaffolds
    //         // Mix into draft assemblies channel with new meta
    //         ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
    //             SPADES_MANUAL.out.scaffolds
    //                 .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'spades') }
    //         )
    //     }

    //     if (!params.skip_kmergenie_strategy) {
    //         SPADES_KMERGENIE(
    //             ch_reads_kmergenie_strategy.map { meta, reads -> [meta, reads, [], []] },
    //             [],
    //             []
    //         )

    //         // Capture output for emits (avoids using conditionals in emits section)
    //         ch_spades_scaffolds_kmergenie = SPADES_KMERGENIE.out.scaffolds

    //         // Mix into draft assemblies channel with new meta
    //         ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
    //             SPADES_KMERGENIE.out.scaffolds
    //                 .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'spades') }
    //         )
    //     }

    //     if (!params.skip_reads_length_strategy) {
    //         SPADES_READS_LENGTH(
    //             ch_reads_reads_length_strategy.map { meta, reads -> [meta, reads, [], []] },
    //             [],
    //             []
    //         )

    //         // Capture output for emits (avoids using conditionals in emits section)
    //         ch_spades_scaffolds_reads_length = SPADES_READS_LENGTH.out.scaffolds

    //         // Mix into draft assemblies channel with new meta
    //         ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
    //             SPADES_READS_LENGTH.out.scaffolds
    //                 .map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'spades') }
    //         )
    //     }
    // }

    // ======= Megahit assemblies - nested conditionals (assembler × strategy) ======
    // Megahit is only run if skip_megahit is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with Megahit assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.
    // Megahit needs a tuple with 3 elements as input. I can't use ch_paired_reads directly because R1 and R2 paths there are in a single list element. So I need to map the channel to split R1 and R2 into separate list elements.
    // MEGAHIT: [ meta, reads1, reads2 ]


    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_megahit_contigs_manual = channel.empty()
    def ch_megahit_contigs_kmergenie = channel.empty()
    def ch_megahit_contigs_reads_length = channel.empty()

    if (!params.skip_megahit) {

        if (!params.skip_manual_strategy) {
            ch_megahit_input_manual = ch_reads_manual_strategy.map { meta, reads -> [meta, reads, []] }
            MEGAHIT_MANUAL(ch_megahit_input_manual)

            // Capture output for emits (avoids using conditionals in emits section)
            ch_megahit_contigs_manual = MEGAHIT_MANUAL.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MEGAHIT_MANUAL.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'megahit') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            ch_megahit_input_kmergenie = ch_reads_kmergenie_strategy.map { meta, reads -> [ meta, reads, [] ] }
            MEGAHIT_KMERGENIE(ch_megahit_input_kmergenie)

            // Capture output for emits (avoids using conditionals in emits section)
            ch_megahit_contigs_kmergenie = MEGAHIT_KMERGENIE.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MEGAHIT_KMERGENIE.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'megahit') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            ch_megahit_input_reads_length = ch_reads_reads_length_strategy.map { meta, reads -> [ meta, reads, [] ] }
            MEGAHIT_READS_LENGTH(ch_megahit_input_reads_length)

            // Capture output for emits (avoids using conditionals in emits section)
            ch_megahit_contigs_reads_length = MEGAHIT_READS_LENGTH.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MEGAHIT_READS_LENGTH.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'megahit') }
            )
        }
    }

    // ======= MINIA assemblies - nested conditionals (assembler × strategy) ======
    // Minia is only run if skip_minia is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with Minia assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.

    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_minia_contigs_manual = channel.empty()
    def ch_minia_contigs_kmergenie = channel.empty()
    def ch_minia_contigs_reads_length = channel.empty()

    if (!params.skip_minia) {

        if (!params.skip_manual_strategy) {
            MINIA_MANUAL(ch_reads_manual_strategy)

            // Capture output for emits (avoids using conditionals in emits section)
            ch_minia_contigs_manual = MINIA_MANUAL.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MINIA_MANUAL.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'minia') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            MINIA_KMERGENIE(ch_reads_kmergenie_strategy)

            // Capture output for emits (avoids using conditionals in emits section)
            ch_minia_contigs_kmergenie = MINIA_KMERGENIE.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MINIA_KMERGENIE.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'minia') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            MINIA_READS_LENGTH(ch_reads_reads_length_strategy)

            // Capture output for emits (avoids using conditionals in emits section)
            ch_minia_contigs_reads_length = MINIA_READS_LENGTH.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                MINIA_READS_LENGTH.out.contigs.map { meta, contigs -> createAssemblyMeta(meta, contigs, 'minia') }
            )
        }
    }

    // ======= ABYSS assemblies - nested conditionals (assembler × strategy) ======
    // ABYSS is only run if skip_abyss is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with ABYSS assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.

    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_abyss_scaffolds_manual = channel.empty()
    def ch_abyss_scaffolds_kmergenie = channel.empty()
    def ch_abyss_scaffolds_reads_length = channel.empty()

    if (!params.skip_abyss) {

        if (!params.skip_manual_strategy) {
            // Join trimmed PE reads with merged SE reads
            ch_abyss_input_manual = ch_fastp_reads
                .join(ch_fastp_reads_merged, by: 0)
                .map { meta, reads_pe, merged ->
                    def new_meta = meta + [kmer_strategy: 'manual', reads_type: 'merged']
                    tuple(new_meta, reads_pe, merged)
                }
            ABYSS_MANUAL(ch_abyss_input_manual, params.abyss_kmer)

            // Capture output for emits (avoids using conditionals in emits section)
            ch_abyss_scaffolds_manual = ABYSS_MANUAL.out.scaffolds

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                ABYSS_MANUAL.out.scaffolds.map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'abyss') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            // Join trimmed PE + merged SE + k-mer value
            ch_abyss_input_kmergenie = ch_fastp_reads
                .map { meta, reads -> [meta.id, meta, reads] }
                .join(ch_fastp_reads_merged.map { meta, merged -> [meta.id, merged] })
                .join(ch_reads_kmergenie_strategy.map { meta, reads -> [meta.id, meta.single_kmer] })
                .map { id, meta, reads_pe, merged, single_kmer ->
                    def new_meta = meta + [
                        kmer_strategy: 'kmergenie',
                        reads_type: 'merged',
                        single_kmer: single_kmer
                    ]
                    tuple(new_meta, reads_pe, merged, single_kmer)
                }

            ABYSS_KMERGENIE(
                ch_abyss_input_kmergenie.map { meta, reads_pe, merged, kmer -> [meta, reads_pe, merged] },
                ch_abyss_input_kmergenie.map { meta, reads_pe, merged, kmer -> kmer }
            )

            // Capture output for emits (avoids using conditionals in emits section)
            ch_abyss_scaffolds_kmergenie = ABYSS_KMERGENIE.out.scaffolds

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                ABYSS_KMERGENIE.out.scaffolds.map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'abyss') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            // Join trimmed PE + merged SE + k-mer value
            ch_abyss_input_reads_length = ch_fastp_reads
                .map { meta, reads -> [meta.id, meta, reads] }
                .join(ch_fastp_reads_merged.map { meta, merged -> [meta.id, merged] })
                .join(ch_reads_reads_length_strategy.map { meta, reads -> [meta.id, meta.single_kmer] })
                .map { id, meta, reads_pe, merged, single_kmer ->
                    def new_meta = meta + [
                        kmer_strategy: 'reads_length',
                        reads_type: 'merged',
                        single_kmer: single_kmer
                    ]
                    tuple(new_meta, reads_pe, merged, single_kmer)
                }

            ABYSS_READS_LENGTH(
                ch_abyss_input_reads_length.map { meta, reads_pe, merged, kmer -> [meta, reads_pe, merged] },
                ch_abyss_input_reads_length.map { meta, reads_pe, merged, kmer -> kmer }
            )

            // Capture output for emits (avoids using conditionals in emits section)
            ch_abyss_scaffolds_reads_length = ABYSS_READS_LENGTH.out.scaffolds

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                ABYSS_READS_LENGTH.out.scaffolds.map { meta, scaffolds -> createAssemblyMeta(meta, scaffolds, 'abyss') }
            )
        }
    }

    // ======= SPARSEASSEMBLER assemblies - nested conditionals (assembler × strategy) ======
    // SPARSEASSEMBLER is only run if skip_sparseassembler is false. Within that, each strategy is only run if its corresponding skip parameter is false.
    // The channel with SPARSEASSEMBLER assemblies is populated accordingly and mixed into the common ch_draft_assemblies_input channel.
    // SPARSEASSEMBLER is a special case because it can output contigs or scaffolds depending on the parameters used and quality of the reads

    // Initialise channels for outputs as empty to avoid use of conditionals in emit section.
    def ch_sparseassembler_scaffolds_manual = channel.empty()
    def ch_sparseassembler_contigs_manual = channel.empty()
    def ch_sparseassembler_scaffolds_kmergenie = channel.empty()
    def ch_sparseassembler_contigs_kmergenie = channel.empty()
    def ch_sparseassembler_scaffolds_reads_length = channel.empty()
    def ch_sparseassembler_contigs_reads_length = channel.empty()

    if (!params.skip_sparseassembler) {

        if (!params.skip_manual_strategy) {
            SPARSEASSEMBLER_MANUAL(
                ch_reads_manual_strategy,
                params.sparseassembler_kmer,
                params.sparseassembler_genome_size,
                params.sparseassembler_expected_coverage
            )

            // Capture output for emits (avoids using conditionals in emits section)
            ch_sparseassembler_scaffolds_manual = SPARSEASSEMBLER_MANUAL.out.scaffolds
            ch_sparseassembler_contigs_manual = SPARSEASSEMBLER_MANUAL.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPARSEASSEMBLER_MANUAL.out.scaffolds
                    .concat(SPARSEASSEMBLER_MANUAL.out.contigs)
                    .unique { meta, assembly -> meta.id }
                    .map { meta, assembly -> createAssemblyMeta(meta, assembly, 'sparseassembler') }
            )
        }

        if (!params.skip_kmergenie_strategy) {
            // Create k-mer channel for SPARSEASSEMBLER (needs single kmer value)
            ch_sparseassembler_kmergenie_single_kmer = ch_reads_kmergenie_strategy.map { meta, reads -> meta.single_kmer }
            SPARSEASSEMBLER_KMERGENIE(
                ch_reads_kmergenie_strategy,
                ch_sparseassembler_kmergenie_single_kmer,
                params.sparseassembler_genome_size,
                params.sparseassembler_expected_coverage
            )

            // Capture output for emits (avoids using conditionals in emits section)
            ch_sparseassembler_scaffolds_kmergenie = SPARSEASSEMBLER_KMERGENIE.out.scaffolds
            ch_sparseassembler_contigs_kmergenie = SPARSEASSEMBLER_KMERGENIE.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPARSEASSEMBLER_KMERGENIE.out.scaffolds
                    .concat(SPARSEASSEMBLER_KMERGENIE.out.contigs)
                    .unique { meta, assembly -> meta.id }
                    .map { meta, assembly -> createAssemblyMeta(meta, assembly, 'sparseassembler') }
            )
        }

        if (!params.skip_reads_length_strategy) {
            // Create k-mer channel for SPARSEASSEMBLER (needs single kmer value)
            ch_sparseassembler_reads_length_single_kmer = ch_reads_reads_length_strategy.map { meta, reads -> meta.single_kmer }
            SPARSEASSEMBLER_READS_LENGTH(
                ch_reads_reads_length_strategy,
                ch_sparseassembler_reads_length_single_kmer,
                params.sparseassembler_genome_size,
                params.sparseassembler_expected_coverage
            )

            // Capture output for emits (avoids using conditionals in emits section)
            ch_sparseassembler_scaffolds_reads_length = SPARSEASSEMBLER_READS_LENGTH.out.scaffolds
            ch_sparseassembler_contigs_reads_length = SPARSEASSEMBLER_READS_LENGTH.out.contigs

            // Mix into draft assemblies channel with new meta
            ch_draft_assemblies_input = ch_draft_assemblies_input.mix(
                SPARSEASSEMBLER_READS_LENGTH.out.scaffolds
                    .concat(SPARSEASSEMBLER_READS_LENGTH.out.contigs)
                    .unique { meta, assembly -> meta.id }
                    .map { meta, assembly -> createAssemblyMeta(meta, assembly, 'sparseassembler') }
            )
        }
    }

    emit:
    // K-mer strategy outputs
    seqkit_stats                                = ch_seqkit_stats
    getseqkitk_kmer                             = ch_getseqkitk_kmer
    kmergenie_html                              = ch_kmergenie_html
    getkmergeniek_k                             = ch_getkmergeniek_k

    // SPAdes outputs
    // spades_scaffolds_manual                     = ch_spades_scaffolds_manual
    // spades_scaffolds_kmergenie                  = ch_spades_scaffolds_kmergenie
    // spades_scaffolds_reads_length               = ch_spades_scaffolds_reads_length

    // MEGAHIT outputs
    megahit_contigs_manual                      = ch_megahit_contigs_manual
    megahit_contigs_kmergenie                   = ch_megahit_contigs_kmergenie
    megahit_contigs_reads_length                = ch_megahit_contigs_reads_length

    // // Minia outputs
    minia_contigs_manual                        = ch_minia_contigs_manual
    minia_contigs_kmergenie                     = ch_minia_contigs_kmergenie
    minia_contigs_reads_length                  = ch_minia_contigs_reads_length

    // ABySS outputs
    abyss_scaffolds_manual                      = ch_abyss_scaffolds_manual
    abyss_scaffolds_kmergenie                   = ch_abyss_scaffolds_kmergenie
    abyss_scaffolds_reads_length                = ch_abyss_scaffolds_reads_length

    // SparseAssembler outputs
    sparseassembler_scaffolds_manual            = ch_sparseassembler_scaffolds_manual
    sparseassembler_contigs_manual              = ch_sparseassembler_contigs_manual
    sparseassembler_scaffolds_kmergenie         = ch_sparseassembler_scaffolds_kmergenie
    sparseassembler_contigs_kmergenie           = ch_sparseassembler_contigs_kmergenie
    sparseassembler_scaffolds_reads_length      = ch_sparseassembler_scaffolds_reads_length
    sparseassembler_contigs_reads_length        = ch_sparseassembler_contigs_reads_length

    draft_assemblies                            = ch_draft_assemblies_input
}
