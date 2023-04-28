# scChICflow

Workflow for processing of [single-cell sortChIC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9925381/) data.

**Author: @vivekbhr**

## DAG (Directed Acyclic Graph) of the Workflow

![DAG](./dag.png)


## Installation and configuration

We assume that the users have python (>=3.8) installed via a conda packlage manager, such as [miniconda](https://docs.conda.io/en/latest/miniconda.html). 


**1. Download this repository**

Move to an appropriate folder and run:

```
git clone https://github.com/vivekbhr/scChICflow.git
```

**2. Go to the cloned directory and set up the tools needed for the workflow**

```
cd scChICflow
conda env create -f env.yaml -n chicflow
```

Note: Setup of this conda environment has been tested with linux, and conda v23. 

### 3. configure the config.yaml

The workflow needs user to specify:
  
  1) path to the (indexed) genome fasta file
  2) path to BWA index of the genome (basename)
  3) path to  cell barcodes (`testdata/chic_384barcodes.txt` file)
  4) path to a bed file with blacklisted regions in the genome (optional)


Copy the `config.yaml` from the folder to your output folder (where you intend to run the pipeline) and replace the information with your relevant information.


## Executing the workflow 

**Test a full run of the workflow with the provided test fasta files.**

This should only take ~10 minutes. Copy the test fastq files provided with the repo (testdata) folder and run the workflow like this:

(provided all files in `config.yaml` are accessible)

```
conda activate chicflow
<scChICflow_folder>/scChICflow -i <testdata_folder> -o <output_folder> -c <your_config.yaml> -j <jobs> -cl
```

here **j** is the number of parallel jobs you want to run, **-cl** means submit to cluster (default is to run locally)

## Expected results



## Notes

### Technical Notes
  - After running the pipeline, **LOG** file are stored in the **<output>/log/** directory and the workflow top-level log is in scChICflow.log file.
  - Currently the -o option is not very flexible and and pipeline works only when it's executed in the output directory.
  - cluster configuration, such as memory and cluster submission command are placed in `cluster_config.yaml`, and can be modified to suite the users internal infrastructure.

## Notes on processing
  - Quality-trimming of the data seems to remove the NLA3 sequences from the 5'-end of R1, which leads to
    rejection of reads during methylation tagging in the NLA-TAPS library. Turn 'trim: False' in the config.yaml to
    turn off trimming for NLA-TAPS library
