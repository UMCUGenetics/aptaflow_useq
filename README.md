# This is a customized version of the aptaflow pipeline

[Aptaflow-pipeline](https://github.com/hovercat/aptaflow)

A pipeline for processing HT-Selex NGS data...

The following modifications have been made in order to be able to run it properly on the inhouse cluster and generated data.

## Lastest versions of packages are installed

```
# Adding of necessary channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels r

# Creating conda environment
conda create --name aptaflow

# Install packages
conda install -c bioconda -c conda-forge nextflow cutadapt fastp multiqc pandas pear
conda install -c r r-base r-dplyr r-tidyr r-ggplot2 r-xlsx r-argparse r-showtext r-openxlsx  bioconductor-biostrings
```

* Add flash (instead of fastp) for merging. Fastp crashes with inhouse data, flash does work and can also be processed by MultiQC (as opposed to Pear)

* Use openxlsx instead of r-xlsx for reporting (which crashes wrt missing java dependency)
   ```
   conda install  -c bioconda -c conda-forge r-openxlsx
   ```

* Add biostrings to environment instead of dynamic install in analyse_selex_composition.r script
   ```
   conda install -c bioconda bioconductor-biostrings
   ```

* Modify process merged_fastq_to_fasta (adjust header lines, was at least required for our inhouse data?)
   ```
   sed -n 'p;n;p;n;n' selex_round.fastq | sed 's/@M/>M/g' > ${round_id}.fasta
   ```
   >
   ```
   sed -n 'p;n;p;n;n' selex_round.fastq | sed 's/@/>/g' > ${round_id}.fasta
   ```
   
* Adjust process trim_selex_primers so that adapter is not forced to be directly at the start (though might be too lenient now?)
   ```
   -g ^${params.primers.p5_f}...${params.primers.p3_f} \
   -G ^${params.primers.p5_r}...${params.primers.p3_r}  \
   ```
   >
   ```
   -g ${params.primers.p5_f}...${params.primers.p3_f} \
   -G ${params.primers.p5_r}...${params.primers.p3_r}  \
   ```
   
## Run pipeline

Fastqs should be unzipped and named as follows:

```
ROUND_SAMPLE_*_R1_001.fastq
ROUND_SAMPLE_*_R2_001.fastq
```

Where ROUND is defined in the config (R01...R10 for example)

```
outdir=/path/to/outdir
/path/to/22.x/nextflow /path/to/aptaflow_uses/aptaflow.nf \
			-c ${outdir}/fle10x.config \
			--merge_mode "flash" \
			--selex_name "<selex_run_name>" \
			--input.raw_reads "/path/to/fastq/*_{R1,R2}_001.fastq" \
			--output.out_dir ${outdir} \
			-w ${outdir}/work/ \
			-resume
```
   
   
# Original readme below...   

More stable/robust pipeliens doing the things below are here:
- Data Preprocessing: https://github.com/hovercat/selex-ngs-prep
- Basic Analysis: https://github.com/hovercat/selex-assess
- K-mer socre based enrichment testing: https://github.com/hovercat/selex-kmer
- Clustering using blast similarity graph: https://github.com/hovercat/selex-blaster (works so so)

# AptaFlow for HT-SELEX NGS Data

This is the bioinformatic pipeline used for the Analysis of *Enterococcus faecalis*, performed in 2020 at TU Wien, Austria.

## Workflow Diagram

Compact workflow: [Workflow Diagram](figures/bioinformatic_wf_compact.jpg)

Detailed workflow: [Workflow Diagram](figures/bioinformatic_wf_detailed.svg)

## Getting Started

To replicate the results for the HT-SELEX against E. faecalis (EF05) the config-file ef05.config in the config folder can be used.
Please note that in the config file the path to the raw NGS files has to be changed first.

### Prerequisites & Installation

The workflow has been developed and run in an [Anaconda](https://www.anaconda.com/understanding-conda-and-pip/) environment.
The conda environment is used for management of dependencies.

First of all, be sure to have at least Python 3.6 installed and install conda:
```
# Installing conda
pip install conda

# Updating conda base environment from default channel
conda update -n base -c defaults conda

# Adding of necessary channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels r
```

Next a conda environment, in which the dependencies will be installed, is created.
```
# Creating conda environment
conda create --name aptaflow
```

Activate the environment and install the dependencies.
```
conda activate aptaflow
conda install -c bioconda -c conda-forge nextflow=19.10.0 cutadapt=2.4 fastp=0.20.0 multiqc=1.9 pandas=1.1.0
conda install -c r r-base=3.6.1 r-dplyr=0.8.0.1 r-tidyr=0.8.3 r-ggplot2=3.1.1 r-xlsx=0.6.1 r-argparse=2.0.1 r-showtext=0.9
```

Packages used:

- nextflow 19.10.0.5170
- cutadapt 2.4
- fastp 0.20.0
- multiqc 1.9
- pandas 1.1.0
- dplyr 0.8.0.1
- r-base 3.6.1
- r-tidyr 0.8.3
- r-ggplot2 3.1.1
- r-xlsx 0.6.1
- r-argparse 2.0.1
- r-showtext 0.9

We had problems using the workflow on a fresh install of Linux (OpenSUSE) due to the xlsx-package using rJava.
For it to properly work the environment variable LD_LIBRARY_PATH has to be set.
```
export LD_LIBRARY_PATH='$JAVA_HOME/lib:$JAVA_HOME/lib/server';
```
Be sure that $JAVA_HOME is pointing at your java install directory.
You can check by echoing the $JAVA_HOME global variable:
```
echo $JAVA_HOME
```

### Running the workflow
Make sure that the scripts in bin/ are executable.
These scripts are called from within the aptaflow.nf script.

If they are not executable, use this command to make them executable.
```
chmod +x bin/*
```

Execute the pipeline like this:
```
nextflow run aptaflow.nf -c configs/YOUR_CONFIG.config
```
Or like this:
```
./aptaflow.nf run -c configs/YOUR_CONFIG.config
```


### Config Files

Make sure to stick to the config file below and change only necessary variables,
such as the round names, data location, random region length and so on.

There are also the config files used in our research for reproducibility.
```
// Parameters for nextflow script aptaflow.nf
params {

    // Provide a selex name for the output directory
    selex_name = "EF05"

    // provide rounds in sorted order
    rounds=["R00", "R02", "R03", "R04", "R05", "R06", "R07", "R08", "R09", "R10", "R11"]

    // random region length
    n=40

    // I/O
    input {
        raw_reads = 'data/*{R1,R2}_001.fastq'
        fwd_suffix="_L001_R1_001.fastq"
        rev_suffix="_L001_R2_001.fastq"
        round_delimiter = "_" // Will be cleaved e.g. R0_XYZ.fastq -> R0
    }
    output {
        out_dir = null
    }

    // System Specs
    specs {
        cpus = 24
    }

    // Primers
    primers {
        p5_f="TAGGGAAGAGAAGGACATATGAT"
        p3_f="TTGACTAGTACATGACCACTTGA"
        p5_r="TCAAGTGGTCATGTACTAGTCAA"
        p3_r="ATCATATGTCCTTCTCTTCCCTA"
    }

    // Trimming
    cutadapt {
        action = "trim"
        max_error = 0.2
        check_primer_contamination = false
        primer_contamination {
            max_error = 0.1
            min_overlap = 12
        }
        N = 40
        primer_length = 23
        max_deviation = 3
    }

    // Filtering and Merging
    fastp {
        filter_min_phred = 30    
    }
    
    // ViennaRNA (Currently not activated)
    vienna_rna {
    	T=21
    	mathews2004="VIENNA_RNA/misc/dna_mathews2004.par"
        for_top=10 // percent [0.00 - 100.00] of every round
    }
}
```

## Affiliated research paper

[Research paper](#) associated with this GitHub Page.

## License

No licensing agreement here yet.


