import argparse
import pathlib
import os
import pandas as pd


def pool_libs(lib_dir: str, output: str):
    """
    Pool libraries in a directory into a single csv file

    Parameters

    :param lib_dir: str
        Path to the directory containing csv files to be pooled
    :param output: str
        Path to the output csv file

    :return: None
    """
    pooled_lib = pd.DataFrame()

    for file in os.listdir(lib_dir):
        if pathlib.Path(file).suffix == ".csv":
            lib = pd.read_csv(os.path.join(lib_dir, file), index_col=0)
            pooled_lib = pd.concat([pooled_lib, lib], axis=0)

    pooled_lib.to_csv(output)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "lib_dir", type=str, help="Path to the directory containing the csv files"
    )
    arg_parser.add_argument("output", type=str, help="Path to the output csv file")

    args = arg_parser.parse_args()

    pool_libs(args.lib_dir, args.output)
