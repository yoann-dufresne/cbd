import sys

# Python script to sort little sequence (10 millions elements) into the A<C<T<G order
# It is a naive sort
# You don't need it to run the library

def sortForACTG(path):
    file = open(path, "r")
    list = []
    size = 0
    for line in file:
        for word in line.split():
            if word != "1":
                size = len(word)
                list.append(encode(word))
    list.sort()
    for i in list:
        print("{}\t1".format(decode(i, size)))
    file.close()

def encode(word):
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

def decode(seq, size):
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

    
