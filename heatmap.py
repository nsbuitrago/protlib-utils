import argparse
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from Bio.Seq import Seq


def get_aa_pos_prob(lib: pd.DataFrame, lib_type: str = "DNA") -> pd.DataFrame:
    """
    Compute the matrix of amino acid or DNA probability by position for a given
    library.

    :param lib_path str: Path to the library csv file
    :param lib_type str: Type of library (DNA or PROTEIN)
        The returned probability matrix will use either DNA or AA based on the
        lib_type.

    :return pd.DataFrame: AA or DNA probability matrix
    """

    lib_size = len(lib)
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    if lib_type == "DNA":
        lib["sequence"] = lib["sequence"].apply(lambda x: str(Seq(x).translate()))

    lib = lib.groupby("sequence").size().reset_index(name="count")

    # Initialize a dictionary to hold the frequency data
    frequency_matrix = {aa: [0] * len(lib["sequence"].iloc[0]) for aa in amino_acids}

    for index, row in lib.iterrows():
        sequence = row["sequence"]
        read_count = row["count"]
        for position, amino_acid in enumerate(sequence):
            if amino_acid in frequency_matrix:
                frequency_matrix[amino_acid][position] += read_count

    # convert to a pd dataframe for easier plotting
    frequency_matrix = pd.DataFrame(frequency_matrix)

    # Normalize the frequencies by the total read counts to get probabilities
    return frequency_matrix.div(lib_size, axis=0)


def generate_heatmap(lib_path: str, lib_type: str = "DNA"):
    """
    Generate and plot a heatmap of amino acid probabilities by position for
    a given library.

    :param lib_path str:
        Path to the library csv file
    :param lib_type str:
        Type of library: DNA or PROTEIN
    """

    prob_matrix = (lib_path, lib_type)

    # Plotting the heatmap
    plt.figure(figsize=(12, 8))
    sb.heatmap(prob_matrix.transpose(), cmap="Blues", annot=False)
    plt.title("Amino Acid Frequency by Position")
    plt.xlabel("Position")
    plt.ylabel("Amino Acid")
    plt.show()


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
    generate_heatmap(lib, args.lib_type)
