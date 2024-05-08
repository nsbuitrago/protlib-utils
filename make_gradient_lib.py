import argparse
import os
import pandas as pd
import random

from Bio import SeqIO


def mutate(seq: str, mut_freq: int) -> str:
    """
    Mutate a sequence by making M random mutations

    :param seq (str): Sequence to mutate
    :param mut_freq (int): Number of mutations to make

    :return (str): Mutated sequence
    """

    codons = [seq[i : i + 3] for i in range(0, len(seq), 3)]
    max_codon_idx = len(codons) - 1
    stop_codons = ["TAA", "TAG", "TGA"]

    for i in range(mut_freq):
        codon_idx = random.randint(0, max_codon_idx)
        codon = codons[codon_idx]
        mut_codon = "TAA"
        nt_idx = random.randint(0, 2)
        nt = codon[nt_idx]

        while mut_codon in stop_codons:
            mutation = random.choice("ATGC".replace(nt, ""))
            mut_codon = codon[:nt_idx] + mutation + codon[nt_idx + 1 :]

        codons[codon_idx] = mut_codon

    return "".join(codons)


def make_lib(seq: str, lib_size: int, mut_freq: int) -> set:
    """
    Make a library from random sampling of L sequences with M
    mutations from the parent sequence

    :param seq (str): Parent sequence
    :param lib_size (int): Library size (L)
    :param mut_freq (int): Mutation frequency (M)

    :return (set): Library of sequences
    """
    lib = set()

    for i in range(lib_size):
        mut_seq = mutate(seq, mut_freq)
        while mut_seq in lib:
            mut_seq = mutate(seq, mut_freq)

        lib.add(mut_seq)

    return lib


def dump_args(args: argparse.Namespace):
    """
    Write (dump) job arguments to a YAML file

    :param args (argparse.Namespace): Job arguments

    :return: None
    """

    output_path = os.path.join(args.output, "job_args.yaml")
    with open(output_path, "w+") as f:
        f.write(
            f"FASTA: {args.fasta}\n"
            f"Mutation frequencies: {args.mut_freqs}\n"
            f"Library size: {args.lib_size}\n"
            f"Enable pooling: {args.pool}\n"
            f"Output path: {args.output}"
        )


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("fasta", type=str, help="Path to fasta file")
    arg_parser.add_argument(
        "-m",
        "--mut_freqs",
        nargs="+",
        type=int,
        help="Mutation frequencies",
        required=True,
    )
    arg_parser.add_argument(
        "-l", "--lib_size", type=int, help="Library size", required=True
    )
    arg_parser.add_argument(
        "-p", "--pool", action="store_true", help="Pool libraries into a single file"
    )
    arg_parser.add_argument(
        "-o", "--output", type=str, help="Output directory", default="output"
    )
    params = arg_parser.parse_args()
    breakpoint()

    parents = SeqIO.parse(params.fasta, "fasta")

    if not os.path.isdir(params.output):
        os.makedirs(params.output)

    pooled_lib = pd.DataFrame()

    for parent in parents:
        for freq in params.mut_freqs:
            target_seq = str(parent.seq).upper()
            lib = make_lib(target_seq, params.lib_size, freq)

            ids = [f"{parent.id}_s{freq}.{i}" for i in range(params.lib_size)]
            lib_df = pd.DataFrame(lib, index=ids, columns=["sequence"])
            lib_df.index.name = "id"

            if params.pool:
                pooled_lib = pd.concat([pooled_lib, lib_df])
            else:
                csv_path = os.path.join(params.output, f"{parent.id}_m{freq}_lib.csv")
                lib_df.to_csv(csv_path)
                dump_args(params)

    if pooled_lib.shape[0] > 0:
        csv_path = os.path.join(params.output, "pooled_lib.csv")
        pooled_lib.to_csv(csv_path)
        dump_args(params)
