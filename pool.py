import argparse
import pathlib
import os
import pandas as pd

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "lib_dir", type=str, help="Path to the directory containing the csv files"
    )
    arg_parser.add_argument("output", type=str, help="Path to the output csv file")

    args = arg_parser.parse_args()

    pooled_lib = pd.DataFrame()

    for file in os.listdir(args.lib_dir):
        if pathlib.Path(file).suffix == ".csv":
            lib = pd.read_csv(os.path.join(args.lib_dir, file), index_col=0)
            pooled_lib = pd.concat([pooled_lib, lib], axis=0)

    pooled_lib.to_csv(args.output)
