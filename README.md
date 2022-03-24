# Berrywood
Berrywood is a NextFlow pipeline used to characterize SARS-CoV-2 genomes and provide run-level characterization data.

### Prerequisites

What dependencies you need to run Berrywood


- [NextFlow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)

### Using the pipeline
The pipeline is designed to start from assembled genomes. All fasta files must be in the same directory. Then start the pipeline using:

```
$ nextflow run ~/berrywood/berrywood.nf --assemblies /path/to/assemblies
```

#### To rename the output pdf

```
$ nextflow run ~/berrywood/berrywood.nf --title "Title name in quotes" --assemblies /path/to/assemblies
```

### Output files
The default directory for Berrywood output data is berrywood_results, unless changed by using the ```--outdir``` flag:
```
$ nextflow run ~/berrywood/berrywood.nf --outdir berrywood_results_2 --assemblies /path/to/assemblies
```

![Berrywood output](/assets/berrywood_output.PNG)

Included in the directory:
- Berrywood-report.pdf
    - a summary report of the lingeages, mutations, pass, and failed samples in the dataset
- berrywood_date.csv
    - a csv file of the mutation, lineage, pass, and fail data
- logs/multi_fasta
    - the multi.fasta file produced by the pipeline
- results/berrywood_results/csvs
    - all csv files produced by Berrywood
- results/nextclade_resutls
    - nextclade results
- results/pangolin_results
    - pangolin results
- results/vadr_results
    - annotations
        - sqa file from vadr
    - csvs
        - Berrywood csv file made from VADR output

### Author


* **Rachael St. Jacques** - *Bioinformatics Principal Scientist* - [stjacqrm](https://github.com/stjacqrm)
