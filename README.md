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
### ConwayBromage::ConwayBromage(istream& kmerFlux, KmerManipulator* km)
Take an istream which represents the k-mer flux and the KmerManipulator which contains all the informations about encoding format.
This constructor will build an object containing all the k-mers in the flux with very low memory consumption.
Example of use: 
```
int kmerSize = 31;
ifstream f("./k-mers.txt", ios::in);
KmerManipulatorACGT km(kmerSize);
ConwayBromage cb(f, &km);        
f.close();
```
#### Requirements
- The k-mers in the istream must me **canonical** and have a **size X <= 32**.
- The isstream must have at each line (if it's a file for example) **the k-mer and its counts**.
Example of a classic txt-file with k-mers of size 4: 
AAAG	1
AAGA	1
AAGC	2
AAGC	1
...

### bool ConwayBromage::isPresent(uint64_t Kmer)
Takes a (k-1)-mer and returns true if it is present among the k-mers.
Example of use:
```
bool kmer_exists = cb.isPresent(230);
```

### uint8_t ConwayBromage::successors(uint64_t Kmer)
Returns the successors of a (k-1)-mer. A (k-1)-mer has at minimum 0 successors and maximum 8 successors.
Informations : 
- The method doesn't check if the Kmer exists (so you have to do it on your own).
- The result can sometimes have duplicate.
Example of use:
```
uint8_t successors = cb.successors(78);
```
#### Explanation on how it works
Let's say the function takes as a parameter the integer which represents the (k-1)-mer GTT and we assume we are in ACGT encoding.
First, the method will generate the following k-mers : GTTA, GTTC, GTTG, GTTT, AGTT, CGTT, GGTT, TGTT.
Then, it will check if the canonical version of these k-mers are present.
Let's suppose that only the canonical version of these k-mers are present : GTTA, GTTT, CGTT, TGTT.
GTTA and GTTT are, what we call **next** k-mers so we will store the information in the 4 left bits of the result (uint8_t).
CGTT and TGTT are **previous** k-mers so, this time, the information will be stored in the 4 right bits of the result.
Thus, the function will return 1001 0101 which corresponds to 149 in base 10.
It means that, in this example, the successors of the (k-1)-mers GTT are TTA, TTT, CGT and TGT.
