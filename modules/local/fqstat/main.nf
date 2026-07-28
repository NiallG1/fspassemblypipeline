process FQSTAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37069a3c615b14eca87a86992d4b8278b399c2ef763da1e139cf5afad886d1af/data' :
        'community.wave.seqera.io/library/perl:5.32.1--ede865d0edc4e459' }"

    input:
    tuple val(meta), path(fastq_files)

    output:
    tuple val(meta), path("${meta.id}*.stats"), emit: stats

    script:
    def fastq_file_list = fastq_files instanceof List ? fastq_files : [fastq_files]
    def fastq_args = fastq_file_list.collect { "\"${it}\"" }.join(' ')
    """
    for fq in ${fastq_args}; do
        zcat "\$fq" | perl ${projectDir}/bin/fq_n50.pl > "${meta.id}_\$(basename "\$fq").stats"
    done
    """

    stub:
    def stat_stubs = [fastq_files].flatten().collect { fq ->
        "touch \"${meta.id}_${fq.getName()}.stats\""
    }.join('\n')
    """
    ${stat_stubs}
    """
}
