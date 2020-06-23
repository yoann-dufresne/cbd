[![Build Status](https://travis-ci.com/yoann-dufresne/ConwayBromageLib.svg?branch=master)](https://travis-ci.com/yoann-dufresne/ConwayBromageLib)

# Conway and Bromage succinct data structures for assembling large genomnes

## Current content

### encode function
**encode** allow us to translate a nucleotide sequence into an int thanks to a given code :
A = 0
C = 1
G = 2
T = 3
We use an hash to create and return unique int. This will help us to complete the sd_vector
The function returns an int that matches with the location of a one in the future sd_vector

### decode function
Returns the string representation of a value representing a k-mer.<br>
<br>
Example of code :
```
int value = 0;
int size = 4; //the size of the kmer
string kmer1 = decode(value, size); // kmer = "AAAA"
string kmer2 = decode(100, 4);      // kmer = "CGCA"
```

### fromFileToSdVector function
The goal of this function is to get sort nucleotide sequences which are store in a tubular file to transcribe them
into a sd_vector (from the sdsl library).
We call **encode** in **fromFileToSdVector** and we create the sd_vector from the location of ones information
At the end, we return an sd_vector which contain 0 and 1. Ones correspond to the location of elements which are in the file (after encode to understand where the location is). <br>
<br>
Example of code : <br>
```
sd_vector<> sdv = fromFileToSdVector("./sorted_kmers.txt"); // the size of sdv is 4^k
```

### successors function
Returns the successors of a k-mer. A k-mer has at maximum 8 successors and minimum 0 successors. <br>
<br>
Example of code : <br>
```
//Let's say that the file in entry (in ACGT encoding) contains the following canonical, and sorted, p-mers : AAG and TAA
sd_vector<> sdv = fromFileToSdVectorChooser("./sorted_kmers.txt", "ACGT");
int Kmer = 0; // representation of AA
vector<uint64_t> successorsOfAA = successors(Kmer, sdv, true); //true = ACGT encoding
//contains {2, 12} which correspond to {AG, TA}
```

### previous function
Returns the predecessors of a k-mer. A k-mer has at maximum 4 predecessors and minimum 0 predecessors.<br> 
<br>
Example of code : <br>
```
//We supposed that the text file is filled with 4-mers.<br>
sd_vector<> sdv = fromFileToSdVector("./sorted_kmers.txt");
vector<string> predecessorsOfCATC = previous("CATC", sdv);
//previous of CATC : { ACAT CCAT GCAT TCAT }
vector<string> predecessorsOfTTTT = previous("TTTT", sdv);
//previous of TTTT : { ATTT CTTT GTTT TTTT }
```

