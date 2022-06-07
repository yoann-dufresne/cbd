
import argparse
from random import sample


def restricted_float(x):
    """ Verify if an argument is a float and is in between 0. and 1. included """
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("Not a floating-point literal")

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("Not in range [0.0, 1.0]")
    return x


def parse_args():
    parser = argparse.ArgumentParser(description='This script generate a tsv file containing fake kmer counts.')
    parser.add_argument('--kmer-size', '-k', type=int, default=31, help='Size of the kmers.')
    parser.add_argument('--filling-ratio', '-r', type=restricted_float, help='The ratio of kmer present in the output.')
    parser.add_argument('--kmer-number', '-n', type=int, default=10, help='The number of kmers in the output. This argument is not used if the ratio is set.')
    
    args = parser.parse_args()

    if args.filling_ratio is not None:
        args.kmer_number = int(pow(4, args.kmer_size) * args.filling_ratio)

    return args


def integer_to_kmer(i, k):
    str_repr = ""
    for _ in range(k):
        str_repr += "ACGT"[i%4]
        i = i // 4

    return str_repr


def generate_kmers(k, num_kmers):
    for i in sample(range(pow(4, k)), num_kmers):
        chaine = integer_to_kmer(i, k)
        print(chaine)


def main():
    args = parse_args()
    generate_kmers(args.kmer_size, args.kmer_number)


if __name__ == "__main__":
    main()
