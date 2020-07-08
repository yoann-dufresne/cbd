[![Build Status](https://travis-ci.com/yoann-dufresne/ConwayBromageLib.svg?branch=master)](https://travis-ci.com/yoann-dufresne/ConwayBromageLib)

# Conway and Bromage succinct data structures for assembling large genomnes

## Object KmerManipulator
Abstract class to manage ACGT or ACTG encoding for k-mers<br>
Mother class of KmerManipulatorACGT and KmerManipulatorACTG

### KmerManipulatorACGT::encode and KmerManipulatorACTG::encode
**encode** allow us to translate a nucleotide sequence into an uint64_t thanks to ACGT or ACTG encoding<br>
We use an hash to create and return unique uint64_t.<br>
The function returns an uint64_t which corresponds to the given string sequence<br>
Example of use :<br>
```
KmerManipulatorACGT encoder(4);	//4-mers in ACGT format
encoder.encode("GGTA");		//172
```

### KmerManipulatorACGT::decode and KmerManipulatorACTG::decode
Returns the string representation of a value representing a k-mer.<br>
Example of use :
```
KmerManipulatorACTG decoder(4);	//4-mers in ACTG format
decoder.decode(248);	//GGTA
```

### KmerManipulatorACGT::getCanonical and KmerManipulatorACTG::getCanonical
Calculate the canonical version of a compressed k-mer (a uint64_t).<br>
The canonical version depends on encoding format.<br>
Example of use :
```
KmerManipulatorACGT canonicaler(4);	//4-mers in ACGT format
canonicaler.getCanonical(172);		//172
KmerManipulatorACTG canonicaler(4);	//4-mers in ACTG format
canonicaler.getCanonical(172);		//144
```
### ConwayBromage::ConwayBromage(istream& kmerFlux, KmerManipulator* km)
Take an istream which represents the k-mer flux and the KmerManipulator which contains all the informations about encoding format.<br>
This constructor will build an object containing all the k-mers in the flux with very low memory consumption.<br>
Example of use: <br>
```
int kmerSize = 31;
ifstream f("./k-mers.txt", ios::in);
KmerManipulatorACGT km(kmerSize);
ConwayBromage cb(f, &km);        
f.close();
```
#### Requirements
- The k-mers in the istream must me **canonical** and have a **size <= 32**.<br>
- The istream must have at each line (if it's a file for example) **the k-mers and its counts**.<br>
<br>
Example of a classic txt-file with k-mers of size 4: <br>
AAAG	1 <br>
AAGA	3 <br>
AAGC	2 <br>
AAGC	1 <br>
...

### bool ConwayBromage::isPresent(uint64_t Kmer)
Takes a (k-1)-mer and returns true if it's present among the stored k-mers.<br>
Example of use:<br>
```
bool kmer_exists = cb.isPresent(230);

KmerManipulatorACGT km(3);
uint64_t intGTT = km.decode("GTT");
bool GTT_exists = cb.isPresent(intGTT);
```

### uint8_t ConwayBromage::successors(uint64_t Kmer)
Returns the successors of a (k-1)-mer. A (k-1)-mer has at minimum 0 successors and maximum 8 successors.<br>
Informations : <br>
- The method doesn't check if the Kmer exists (so you have to do it on your own).<br>
- The result can sometimes have duplicate.<br>
<br>

#### Explanation on how it works
Let's say the function takes as a parameter the integer which represents the (k-1)-mer **GTT** and we assume we are in ACGT encoding.<br>
First, the method will generate the following k-mers : GTTA, GTTC, GTTG, GTTT, AGTT, CGTT, GGTT, TGTT.<br>
Then, it will check if the canonical version of these k-mers are present.<br>
Let's suppose that only the canonical version of these k-mers are present : GTTA, GTTT, CGTT, TGTT.<br>
GTTA and GTTT are, what we call **next** k-mers so we will store the information in the 4 left bits of the result (uint8_t).<br>
CGTT and TGTT are **previous** k-mers so, this time, the information will be stored in the 4 right bits of the result.<br><br>
Thus, the function will return 1001 0101 which corresponds to 149 in base 10.<br>
It means that, in this example, **the successors of the (k-1)-mer GTT are TTA, TTT, CGT and TGT**.<br>

Example of use: <br>
```
uint8_t successors = cb.successors(78);

KmerManipulatorACGT km(3);
uint64_t intGTT = km.decode("GTT");
uint8_t successorsOfGTT = cb.successors(intGTT);
```
