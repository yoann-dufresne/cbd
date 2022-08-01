from Bio import SeqIO
import random
import sys
def main():
    printlist(test(500,1000))

def test(size,number):
    list=[]
    for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
        for i in range(number):
            x=random.randrange(len(seq_record))
            list.append(str(seq_record.seq[x:x+size]))
    return list
def printlist(list):
    for a in list:
        print(a+"\n")
    return 0


if __name__ == "__main__":
    main()

