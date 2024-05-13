from heatmap import get_aa_pos_prob
import numpy as np
import argparse


def get_lib_entropy(lib_path: str, lib_type: str) -> int:
    """
    Compute the entropy of a given library.

    :param lib_path str: Path to the library csv file
    :param lib_type str: Type of library (DNA or PROTEIN)

    :return int: Entropy of the library
    """

    prob_matrix = get_aa_pos_prob(lib_path, lib_type)
    entropy = -1 * (prob_matrix * np.log2(prob_matrix)).sum().sum()
    return entropy


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Generate heatmap for library")
    arg_parser.add_argument("lib_path", type=str, help="Path to library csv file")
    arg_parser.add_argument(
        "lib_type",
        type=str,
        help="Type of library: DNA or PROTEIN",
        choices=["DNA", "PROTEIN"],
        default="DNA",
    )

    args = arg_parser.parse_args()

    entropy = get_lib_entropy(args.lib_path, args.lib_type)
    print(f"Library Entropy (bits): {entropy}")
