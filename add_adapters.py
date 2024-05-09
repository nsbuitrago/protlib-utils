import argparse
from Bio import SeqIO
import pandas as pd


def add_adapters(adapter_path: str, lib_path: str, output_path: str):
    """
    Add adapters to library

    :param adapter_path: str
        Path to adapter fasta
    :param lib_path: str
        Path to library csv
    :param output_path: str
    """
    adapters = list(SeqIO.parse(adapter_path, "fasta"))
    adapters = (str(adapters[0].seq).upper(), str(adapters[1].seq).upper())
    breakpoint()
    lib = pd.read_csv(lib_path, index_col="name")
    lib["sequence"] = lib["sequence"].apply(lambda x: adapters[0] + x + adapters[1])

    lib.to_csv(args.output)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("lib", type=str, help="Path to library csv")
    arg_parser.add_argument("adapters", type=str, help="Path to adapter fasta")
    arg_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="output.csv",
        help="Path to output csv",
    )
    args = arg_parser.parse_args()

    add_adapters(args.adapters, args.lib, args.output)
