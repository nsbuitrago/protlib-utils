import argparse
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
from Bio.Seq import Seq

if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Generate heatmap for library")

    arg_parser.add_argument("lib_path", type=str, help="Path to library csv file")

    args = arg_parser.parse_args()

    lib = pd.read_csv(args.lib_path)
    lib_size = len(lib)
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    # translate DNA -> protein
    lib["protein"] = lib["sequence"].apply(lambda x: str(Seq(x).translate()))
    # get count of unique sequence
    lib = lib.groupby("protein").size().reset_index(name="count")

    # Initialize a dictionary to hold the frequency data
    frequency_matrix = {aa: [0] * len(lib["protein"].iloc[0]) for aa in amino_acids}

    # Calculate the frequency of each amino acid at each position
    for index, row in lib.iterrows():
        sequence = row["protein"]
        read_count = row["count"]
        for position, amino_acid in enumerate(sequence):
            if amino_acid in frequency_matrix:
                frequency_matrix[amino_acid][position] += read_count

    # Convert the frequency matrix to a DataFrame for easier plotting
    frequency_lib = pd.DataFrame(frequency_matrix)

    # Normalize the frequencies by the total read counts to get probabilities
    frequency_lib = frequency_lib.div(lib_size, axis=0)

    # Plotting the heatmap
    plt.figure(figsize=(12, 8))
    sb.heatmap(frequency_lib.transpose(), cmap="Blues", annot=False)
    plt.title("Amino Acid Frequency by Position")
    plt.xlabel("Position")
    plt.ylabel("Amino Acid")
    plt.show()
