from Bio import SeqIO
import random
import sys
def main():
    print(test(500,1000))

def test(size,number):
    list=[]
    for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
        for i in range(number):
            x=random.randrange(len(seq_record))
            list.append(str(seq_record.seq[x:x+size]))
    return list



if __name__ == "__main__":
    main()

