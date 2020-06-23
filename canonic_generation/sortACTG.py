import sys
#sort element of ACTG format
def sortForACTG(path):
    file = open(path, "r")
    list = []
    size = 0
    for line in file:   #reading line by line
        for word in line.split():
            if word != "1":
                size = len(word)
                list.append(encode(word))   #collect each k-mers in compressed version (uint64_t)
    list.sort() #We can sort numbers no matter what the format is
    for i in list:
        print("{}\t1".format(decode(i, size)))  #print of uncompressed version (str)
    file.close()

def encode(word):   #Code to compress a str k-mer (from str to uint64_t)
    hash = 0
    for letter in word:
        hash <<= 2
        charval = 0
        c = letter
        if c == "C":
            charval = 1
        elif c == "T":
            charval = 2
        elif c == "G":
            charval = 3
        hash += charval
 
    return hash

def decode(seq, size):  #COde to uncompress a uint64_t k-mer (from uint64_t to str)
    res = ""
    lastIndex = size-1
    for i in range(size):
        if seq & 0x3 == 0:
            res += 'A'
        elif seq & 0x3 == 1:
            res += 'C'
        elif seq & 0x3 == 2:
            res += 'T'
        elif seq & 0x3 == 3:
            res += 'G'
        seq >>= 2
    res = "".join(reversed(res))
    return res

def main():
    sortForACTG(sys.argv[1])

    
if __name__ == "__main__":
    main()

    
