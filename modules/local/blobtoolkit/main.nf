process BLOBTOOLKIT_CREATEBLOBDIR {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta1), path(busco, stageAs: 'lineage??/*')
    tuple val(meta3), path(yaml)
    path(taxdump, stageAs: 'taxdump/taxdump.json')

    output:
    tuple val(meta), path(prefix), emit: blobdir
    tuple val("${task.process}"), val('blobtoolkit'), eval("btk --version | cut -d' ' -f2 | sed 's/v//'"), topic: versions, emit: versions_blobtoolkit

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "BLOBTOOLKIT_BLOBDIR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def busco_args = (busco instanceof List ? busco : [busco]).collect { file -> "--busco " + file } .join(' ')
    def hits_blastp = blastp ? "--hits ${blastp}" : ""
    """
    blobtools replace \\
        --bedtsvdir windowstats \\
        --meta ${yaml} \\
        --taxdump \$(dirname ${taxdump}) \\
        --taxrule buscogenes \\
        ${busco_args} \\
        ${hits_blastp} \\
        --threads ${task.cpus} \\
        $args \\
        ${prefix}

   
    """
   
    stub:
    def args = task.ext.args ?: ''
    
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo $args
    
    touch ${prefix}.bam
    """
}
