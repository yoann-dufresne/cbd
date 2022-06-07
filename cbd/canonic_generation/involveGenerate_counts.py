
import argparse
from random import sample

#Remake of yoann's generate_counts.py : now able to manage canonical nucleotides in ACGT and ACTG format

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
    parser.add_argument('--format-encode', '-f', type=str, default="ACGT", help='Encoding format of the sequence, need ACGT or ACTG')   #format choice : ACGT or ACTG
    
    args = parser.parse_args()

    if args.filling_ratio is not None:
        args.kmer_number = int(pow(4, args.kmer_size) * args.filling_ratio)

    return args


def integer_to_kmer(i, k, format_encode):
    str_repr = ""
    if format_encode == "ACGT":
        for _ in range(k):
            str_repr += "ACGT"[i%4] #creation of k-mers, ACGT format
            i = i // 4
    else:
        for _ in range(k):
            str_repr += "ACTG"[i%4] #creation of k-mers, ACTG format
            i = i // 4
        str_repr = "".join(reversed(str_repr))
        #print("str_repr : {}".format(str_repr))
    return str_repr


def generate_kmers(k, num_kmers, format_encode):
    for i in sample(range(pow(4, k)), num_kmers):   #generate num_kmers randomlu between 0 and 4**k excluded
        #randomer = randint(0, pow(4, k)-1)
       # print("randomer : {}".format(randomer))
        cano = getCanonical(i, k, format_encode) 
        print("{}\t1".format(integer_to_kmer(cano, k, format_encode)))

def reverseACGT(mer, kmerSize): #same ACGT reverse as in the library, python version
    res = ~mer

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16)
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32)

    return (res >> (2 * (32 - kmerSize)))

def reverseACTG(x, sizeKmer):   #same ACTG reverse as in the library, python version
    res = x

    res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
    res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
    res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
    res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2*( 32 - sizeKmer))) ;

def getCanonical(i, k, format_encode):  #return the canonical version of the k-mers, return itself if it is canonical already
    if format_encode == "ACGT":
        rev = reverseACGT(i, k)
        if rev < i :
            return rev
        else:
            return i
    if format_encode == "ACTG":
        rev = reverseACTG(i, k)
        if rev < i :
            return rev
        else:
            return i

def main():
    args = parse_args()
    generate_kmers(args.kmer_size, args.kmer_number, args.format_encode)

    
if __name__ == "__main__":
    main()
