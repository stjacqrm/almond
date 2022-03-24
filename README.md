# Almond
Almond is a NextFlow pipeline used to analyze sequencing competency runs.

### Prerequisites

What dependencies you need to run Almond


- [NextFlow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/)

### Using the pipeline
The pipeline is designed to start from assembled genomes. All fasta files must be in the same directory. Then start the pipeline using:

```
$ nextflow run ~/almond/almond.nf --assemblies /path/to/assemblies
```

#### To rename the output pdf

```
$ nextflow run ~/almond/almond.nf --title "Title name in quotes" --assemblies /path/to/assemblies
```

### Output files
The default directory for Almond output data is almond_results, unless changed by using the ```--outdir``` flag:
```
$ nextflow run ~/almond/almond.nf --outdir almond_results_2 --assemblies /path/to/assemblies
```

### Author


* **Rachael St. Jacques** - *Bioinformatics Principal Scientist* - [stjacqrm](https://github.com/stjacqrm)
