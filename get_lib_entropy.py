from heatmap import get_aa_pos_prob
import pandas as pd
import numpy as np
import argparse
import unittest


def get_lib_entropy(lib: pd.DataFrame, lib_type: str) -> int:
    """
    Compute the entropy of a given library.

    :param lib_path str: Path to the library csv file
    :param lib_type str: Type of library (DNA or PROTEIN)

    :return int: Entropy of the library
    """

    prob_matrix = get_aa_pos_prob(lib, lib_type)
    # ignore 0 to avoid log(0) error and 1 to avoid log(1) = 0
    valid_prob_matrix = prob_matrix[(prob_matrix > 0) & (prob_matrix < 1)]
    entropy = -1 * (valid_prob_matrix * np.log2(valid_prob_matrix)).sum().sum()
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
    lib = pd.read_csv(args.lib_path)
    entropy = get_lib_entropy(lib, args.lib_type)

    print(f"Library Entropy (bits): {entropy}")


class TestLibEntropy(unittest.TestCase):
    def test_get_lib_entropy(self):
        lib_a = pd.DataFrame({"sequence": ["MAGICAL", "MAGIKAL"], "count": [1, 1]})
        lib_b = pd.DataFrame({"sequence": ["MAGICAL", "LIGGAND"], "count": [1, 1]})

        entropy_a = get_lib_entropy(lib_a, "PROTEIN")
        entropy_b = get_lib_entropy(lib_b, "PROTEIN")

        self.assertGreater(entropy_b, entropy_a)
