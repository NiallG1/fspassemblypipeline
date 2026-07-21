<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-fspassemblypipeline_logo_dark.png">
    <img alt="nf-core/fspassemblypipeline" src="docs/images/nf-core-fspassemblypipeline_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/fspassemblypipeline)
[![GitHub Actions CI Status](https://github.com/nf-core/fspassemblypipeline/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/fspassemblypipeline/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/fspassemblypipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/fspassemblypipeline/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/fspassemblypipeline/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/fspassemblypipeline)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23fspassemblypipeline-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/fspassemblypipeline)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/fspassemblypipeline** is a comprehensive bioinformatics pipeline for fungal genome assembly from Illumina short-read sequencing data. The pipeline ingests raw paired-end reads and performs quality control, read preprocessing (trimming, merging, removing clean reads less than 30bp and deduplication), k-mer profiling (sequencing depth and genome size estimation), de novo genome assembly (with multiple assembler and multiple kmer strategy options), genome assembly quality assessment (BUSCO, QUAST), and contamination detection. It is designed to handle challenging samples such as those with degraded DNA from fungal herbarium specimens, as implemented for the Fungarium Sequencing Project at Royal Botanic Gardens, Kew (https://www.kew.org/science/our-science/projects/sequencing-kews-fungarium).

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->

## Pipeline steps

1. Read QC ([`Falco`](https://github.com/smithlabcode/falco) and [`MultiQC`](http://multiqc.info/))
2. Read preprocessing ([`fastp`](https://github.com/OpenGene/fastp) short reads trimming and merging)
3. K-mer counting ([`FASTK`](https://github.com/thegenemyers/FASTK) and k-mer profiling [`genescopeFK`] (https://github.com/thegenemyers/GENESCOPE.FK))
3. Genome assembly ([`SPAdes`](https://github.com/ablab/spades), [`MEGAHIT`](https://github.com/voutcn/megahit), [`Minia`](https://github.com/GATB/minia), [`ABySS`](https://github.com/bcgsc/abyss), [`SparseAssembler`](https://github.com/yechengxi/SparseAssembler))
4. Assembly assessment ([`BUSCO`](https://busco.ezlab.org/), [`QUAST`](http://quast.sourceforge.net/), [`MerquryFK`](https://github.com/thegenemyers/MerquryFK))
5. Contamination detection and filtering ([`Tiara`](https://github.com/ibe-uw/tiara),[`FCS-GX`](https://github.com/ncbi/fcs-gx), [`BlobTools`](https://github.com/drl/blobtools))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

### Required inputs

#### Pre-dowloaded BUSCO lineages

To run BUSCO using a lineage closely related to each of the samples, we need to input a list of available busco lineages, and to download them beforehand. This can be achieved as follows:

```
cd to/where/you/want/to/store/busco/databases

conda activate busco

busco --list > busco_lineages.txt

gawk '/fungi_odb12/{flag=1; indent=length($0)-length(ltrim($0)); print "fungi_odb12"; next}
     flag && /- [a-z_]*_odb12/ {
         current_indent=length($0)-length(ltrim($0))
         if(current_indent <= indent) flag=0
         else print gensub(/.*- ([a-z_]*_odb12).*/, "\\1", "g")
     }
     function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }' busco_lineages.txt > fungi_busco_lineages.txt
```

In the example above we are extracting the names of all the BUSCO lineages that belong to the fungi kingdom. The target group can be different and its taxonomic level doesn't matter. The user can for example target `eukaryota` or something more specific like `basidiomycota` using the same code. Note that `odb12` extension refers to a specific version of BUSCO lineages, and it can be changed when newer versions will be available.

`fungi_busco_lineages.txt` and the extension must to be provided through the `nextflow.config`:

```
    busco_db_extension         = 'odb12'
    lineages_list_file         = 'path/to/fungi_busco_lineages.txt'
```

Using the list of lineages of interest we can then easily download all of them in one go:

```
for i in $(cat fungi_busco_lineages.txt); do
  echo "downloading $i database"
  busco --download_path . --download $i
done
```

This speeds up the pipeline as it will not have to download busco lineages on the fly, and will avoid connection problems during the run.

We also need to provide the path to where busco lineages are downloaded in `nextflow.config`:

```
    busco_lineages_path        = 'path/to/lineages/parent/directory'
```

Note that BUSCO automatically downloads lineages in a directory called `lineages`. In `nextflow.config` we need to provide the path to the parent directory of `lineages`. This needs to be the full absolute path.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

### Run the pipeline

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/fspassemblypipeline \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Preprocessing subworkflow

The preprocessing subworkflow is implemented in `subworkflows/local/preprocessing/main.nf` and is exposed for local debugging via `preprocessing_only.nf`.
It performs raw read QC, adapter trimming, read merging, QC compilation, and k-mer profiling before downstream assembly and analysis.

### Overview

The preprocessing subworkflow runs the following steps:
- `FALCO` raw read QC
- `fastp` trimming/filtering while keeping complete trimmed R1/R2 output
- `fastp` merge of trimmed reads with merged and unmerged outputs
- `FALCO` QC on trimmed and merged reads
- Falco QC statistics compilation across raw, trimmed, and merged stages
- `FQSTAT` read statistics and summary report generation
- `FASTK` k-mer histogram generation
- `GENESCOPEFK` k-mer profile summarization
- final k-mer summary table generation

### Required inputs

The subworkflow accepts an input samplesheet via `--input`.
It supports raw paired-end reads files that end in .fq.gz.

Optional inputs:
- `--fastp_adapter_fasta path/to/adapters.fa` for custom adapter sequences

### Run the certain subworkflows only

This is a useful entrypoint for running only the preprocessing stage independently or for local debugging.

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/fspassemblypipeline/usage) and the [parameter documentation](https://nf-co.re/fspassemblypipeline/parameters).


## Decontamination subworkflow

The Decontamination subworkflow is implemented in `subworkflows/local/contamination_detection/main.nf`.
It runs FCS-GX and Tiara for taxanomic labelling of contigs and scaffolds, reformats the FCS-GX output for processing, creates a yaml file for each sample. The output is the handed to the blobtools subworkflow.

### Overview

The contamination detection subworkflow runs the following steps:
- `Tiara` assigns domain level taxonomy and organelle/motrochondrial DNA labels to contigs.
- `FCS-GX` Assigns species level taxonomy to contigs.
- `convertrpt` Reformats output of FCS-GX for downstream processing.
- `Comparison` Compares the domain level assignments of Tiara & FCS-GX and uses this to create labels for the final blobplot.
- `Create_yaml` Creates a yaml file from the samplesheet to generate the blobplot

### Required inputs

The subworkflow accepts an input samplesheet via `--input`.
It supports genome assembly fasta files that end in .fa.gz.
requires an installation of the FCS-GX database.

Optional inputs:

## blobtools subworkflow

The blobtools subworkflow is implemented in `subworkflows/local/blobtools/main.nf`.
It takes the taxonomic labels created in the contamination_detection subworkflow and the GC content and coverage information and
plots this on a graph allowing for visualisation of contamination.

### Overview

The blobtools subworkflow runs the following steps:
- `Create_yaml` Creates a yaml file from the samplesheet to generate the blobplot.
- `SAMTOOLS_CSI` Indexes the .bam files to produce .bam.csi files as required by blobtools.
- `blobtoolkit_create` Creates the blobdir from the output of contamination detection and the provided samplesheet.

### Required inputs

needs a .bam coverage file and the output from contamination detection.

Optional inputs:

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/fspassemblypipeline/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/fspassemblypipeline/output).

## Credits

nf-core/fspassemblypipeline was originally written by Lia Obinu, Niall Garvey, Wu Huang, Chris Wyatt, Fernando Duarte Frutos.

We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#fspassemblypipeline` channel](https://nfcore.slack.com/channels/fspassemblypipeline) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/fspassemblypipeline for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
