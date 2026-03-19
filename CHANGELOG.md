# nf-core/fspassemblypipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/fspassemblypipeline, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- 09/02/2026 - Initialised GENOME_ASSEMBLY subworkflow.
- 09/02/2026 - Added SEQKIT_STATS module.
- 10/02/2026 - Added FASTK_FASTK module.
- 10/02/2026 - Added SPADES module.
- 10/02/2026 - Added MEGAHIT module.
- 10/02/2026 - Added MINIA module.
- 12/02/2026 - Added RENAME_ASSEMBLIES local module.
- 12/02/2026 - Added BUSCO_BUSCO module.
- 12/02/2026 - Added MERQURYFK_MERQURYFK module.
- 12/02/2026 - Added QUAST module.
<<<<<<< kmergenie_module
- 25/02/2026 - Added kmergenie local module.
- 25/02/2026 - Added getkmergeniek local module.
- 27/02/2026 - Created nf-core module for kmergenie
- 03/03/2026 - Added kmergenie to the pipeline as nf-core module
=======
- 17/03/2026 - Added sparseassembler as local module
>>>>>>> dev

### `Fixed`

- 12/02/2026 - Renamed assemblies to avoid conflicts in downstream modules.
- 12/02/2026 - New output directory structure.
<<<<<<< kmergenie_module
- 12/02/2026 - User can set extra params from `nextflow.config`.
- 13/03/2026 - fixed kmergenie nf-core module (missing log)
- 17/03/2026 - removed tests and meta.yml from getkmergeniek local module as it's not needed
=======
- 12/02/2026 - User cans set extra params from `nextflow.config`.
>>>>>>> dev

### `Dependencies`

### `Deprecated`
