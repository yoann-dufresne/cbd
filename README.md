[![Build Status](https://travis-ci.com/yoann-dufresne/ConwayBromageLib.svg?branch=master)](https://travis-ci.com/yoann-dufresne/ConwayBromageLib)

# Conway and Bromage succinct data structures for assembling large genomnes

## Object KmerManipulator
Abstract class to manage ACGT or ACTG encoding for k-mers<br>
Mother class of KmerManipulatorACGT and KmerManipulatorACTG

### KmerManipulatorACGT::encode and KmerManipulatorACTG::encode
**encode** allow us to translate a nucleotide sequence into an uint64_t thanks to ACGT or ACTG encoding<br>
We use an hash to create and return unique uint64_t.<br>
The function returns an uint64_t which corresponds to the given string sequence<br>
Example of code :<br>
```
KmerManipulatorACGT encoder(4);	//4-mers in ACGT format
encoder.encode("GGTA");		//172
```

### KmerManipulatorACGT::decode and KmerManipulatorACTG::decode
Returns the string representation of a value representing a k-mer.<br>
Example of code :
```
KmerManipulatorACTG decoder(4);	//4-mers in ACTG format
decoder.decode(248);	//GGTA
```

### KmerManipulatorACGT::getCanonical and KmerManipulatorACTG::getCanonical
Calculate the canonical version of a compressed k-mer (a uint64_t).<br>
The canonical version depends on encoding format.<br>
Example of code :
```
KmerManipulatorACGT canonicaler(4);	//4-mers in ACGT format
canonicaler.getCanonical(172);		//172
KmerManipulatorACTG canonicaler(4);	//4-mers in ACTG format
canonicaler.getCanonical(172);		//144
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

