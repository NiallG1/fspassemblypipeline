process CONVERTFCSRPT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/coreutils:8.25--2471b967344e8d86"

    input:
    tuple val(meta), path(rpt)

    output:
    tuple val(meta), path("*.tsv"), emit: fcs_report_reformatted

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cut --complement -f 5,12,18,24,30 ${rpt} \\
        | tail -n +2 \\
        | sed '1s/^#//' \\
        > ${prefix}.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(cut --version | head -n1 | sed 's/cut (GNU coreutils) //')
    END_VERSIONS
    """

}