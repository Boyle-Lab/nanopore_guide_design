# Nanopore Guide Design
## Nextflow implementation for generating, filtering, and prioritizing candidate guide RNA sequences for Nanopore sequencing

## Requirements
- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- conda distribution (Anaconda or Miniconda) [https://www.anaconda.com](https://www.anaconda.com)
- [BEDtools](https://bedtools.readthedocs.io/en/latest/)
- [CHOPCHOP](https://bitbucket.org/valenlab/chopchop/src/master/) (use the README [here](https://github.com/dpear/guidesign) for installation/configuration)
- CRISPRon scripts (need to create repo for these)

## Configuration
Create conda environments

```conda env create -f chopchop.yml && conda env create -f crispron.yml```

Update ```nextflow.config```

```
output = "/path/to/output"
chopchop_path = "/path/to/chopchop.py"
bedtools_path = "/path/to/bedtools"
crispron_path = "/path/to/crispron"
```

Update ```guide_design.nf```

```
process parallel_chopchop
conda '/path/to/conda/env/chopchop'

process parallel_crispron
conda `/path/to/conda/env/crispron'
```

## Usage
To run Nextflow pipeline from the repository directory:

```nextflow run guide_design.nf -output-dir /path/to/output```

To test different paramaters (eg. distance) run:

```nextflow run guide_design.nf -output-dir /path/to/output --distance 2000```
