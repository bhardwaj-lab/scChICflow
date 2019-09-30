# openTAPS
Workflow for single-cell TAPS analysis


**Author: @vivekbhr**

## DAG (Directed Acyclic Graph) of the Workflow

![DAG](./dag.png)


## How to run

1. Clone the repo

```
git clone https://github.com/vivekbhr/openTAPS.git
```

2. Go to the cloned directory and set up conda env for the workflow

```
cd openTAPS
conda env create -f env.yaml -n taps
```

3. configure the config.yaml

The workflow needs
1) path to the (indexed) genome fasta file
2) path to BWA index (base)
3) path to taps cell barcodes
4) Name of "method", ('chic', 'chic-taps' or 'nla-taps')

Copy the config.yaml from the folder to your output folder (where you intend to run the pipeline) and replace the information with your relevant information.

4. Test workflow with test fasta.

This should only take 3-4 minutes. Copy the test fastq files provided with the repo (testdata) folder and run the workflow like this.

```
conda activate taps
<openTAPS_folder>/openTAPS -i <testdata_folder> -o <output_folder> -c <your_config.yaml> -j <jobs> -cl
```

here **j** is the number of parallel jobs you want to run, **-cl** means submit to cluster (default is to run locally)

After running the pipeline, **LOG** file are stored in the **/log/** directory and the workflow top-level log is in openTAPS.log file.
