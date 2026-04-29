process BLOBTOOLKIT_CREATEBLOBDIR {
    tag "$meta.id"
    label 'process_medium'

    
    conda "${moduleDir}/environment.yml"
    container "docker.io/genomehubs/blobtoolkit:4.4.6"

    input:
    tuple val(meta), path(fasta), path(bam), path(busco) //add fasta
    //tuple val(meta), path(bam)  // add bam file
    //tuple val(meta), path(busco, stageAs: 'lineage??/*') //add busco, why stageas?
    //tuple val(meta3), path(yaml) // will need to create a yaml somewhere?
    
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
    //def busco_args = (busco instanceof List ? busco : [busco]).collect { file -> "--busco " + file } .join(' ')
    
   
    """
    blobtools create \\
        --fasta ${fasta} \\
        --threads ${task.cpus} \\
        --cov ${bam} \\
        $args \\
        ${prefix}

    """
   
    stub:
    def args = task.ext.args ?: ''
  
    """
    echo $args
    
    touch ${prefix}.bam
    """
}
