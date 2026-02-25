process KMERGENIE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmergenie:1.7051--py311r44hc84137b_9':
        'biocontainers/kmergenie:1.7051--py311r44hc84137b_9' }"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("*_report.html"), emit: html
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $reads1 > ${prefix}_reads.txt
    echo $reads2 >> ${prefix}_reads.txt

    kmergenie \\
        $args \\
        -o ${prefix} \\
        -t $task.cpus \\
        ${prefix}_reads.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmergenie: \$(echo \$(kmergenie --version 2>&1) | sed 's/KmerGenie //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args
    
    touch ${prefix}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmergenie: \$(echo \$(kmergenie --version 2>&1) | sed 's/KmerGenie //')
    END_VERSIONS
    """
}
