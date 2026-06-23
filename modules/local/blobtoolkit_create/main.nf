process BLOBTOOLKIT_CREATEBLOBDIR {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_gcc_linux-64_gxx_linux-64_sysroot_linux-64_pruned:c50e1fdcdb18252f"

    input:
    tuple val(meta), path(fasta), path(bam), path(yaml), path(busco), path(index), path(taxonomy)

    output:
    tuple val(meta), path("${meta.id}"), emit: blobdir
    tuple val("${task.process}"), val('blobtoolkit'), eval("btk --version | cut -d' ' -f2 | sed 's/v//'"), topic: versions, emit: versions_blobtoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    //if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
    //    exit 1, "BLOBTOOLKIT_BLOBDIR module does not support Conda. Please use Docker / Singularity / Podman instead."
    //}

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_arg = busco ? "--busco ${busco}" : ""

    """
    blobtools create \\
        --fasta ${fasta} \\
        --meta ${yaml} \\
        --threads ${task.cpus} \\
        ${busco_arg} \\
        ${args} \\
        ${prefix}

    blobtools add \\
        --cov ${bam} \\
        --threads ${task.cpus} \\
        ${prefix}

    blobtools add \\
       --text ${taxonomy} \\
       --text-delimiter '\t' \\
       --text-cols 'seq_id=identifiers,taxonomy=taxonomy' \\
       --text-header \\
        --key plot.cat=taxonomy \\
        ${prefix}

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/meta.json
    touch ${prefix}/identifiers.json
    """
}