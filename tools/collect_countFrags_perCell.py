import pandas as pd
import os
import glob
import argparse

# Load {sample}.per_barcode.tsv files and write to single tsv
def parse_arguments():
    parser = argparse.ArgumentParser(description="Aggregate .per_barcode.tsv files into a single file.")
    parser.add_argument("barcodes", type=str, help="Path to the barcodes file")
    parser.add_argument("input_dir", type=str, help="Directory containing input .per_barcode.tsv files")
    parser.add_argument("output_file", type=str, help="Path to the output aggregated file")
    return parser.parse_args()

def aggregate_files(barcodes, input_dir, output_file):
    all_data = pd.read_csv(barcodes, sep='\t', header=None)
    all_data.columns = ['cell_barcode']
    
    for file in glob.glob(os.path.join(input_dir, "*.per_barcode.tsv")):
        sample_name = os.path.basename(file).replace(".per_barcode.tsv", "").split(".")[-1]

        df = pd.read_csv(file, sep="\t", header=None)
        df.columns = ["cell_barcode", sample_name]
        all_data = all_data.merge(df, on="cell_barcode", how="outer")

    all_data.to_csv(output_file, sep="\t", index=False)

if __name__ == "__main__":

    args = parse_arguments()

    barcodes = args.barcodes
    input_dir = args.input_dir
    output_file = args.output_file

    if not os.path.exists(input_dir):
        raise FileNotFoundError(f"Input directory {input_dir} does not exist.")
    
    aggregate_files(barcodes, input_dir, output_file)
