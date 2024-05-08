import argparse
from Bio import SeqIO
import pandas as pd


def add_adapters(seq: str, adapters: tuple[str, str]) -> str:
    return adapters[0] + seq + adapters[1]


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

    adapters = list(SeqIO.parse(args.adapters, "fasta"))
    adapters = (str(adapters[0].seq), str(adapters[1].seq))
    lib = pd.read_csv(args.lib, index_col="id")
    lib["sequence"] = lib["sequence"].apply(lambda x: add_adapters(x, adapters))

    lib.to_csv(args.output)
