# ConwayBromageLib : an implementation of Conway and Bromage succinct data structure for assembling large genomes
[![Build Status](https://travis-ci.com/yoann-dufresne/ConwayBromageLib.svg?branch=master)](https://travis-ci.com/yoann-dufresne/ConwayBromageLib)
## What is ConwayBromageLib ?
The ConwayBromageLib (CBL) is a C++11 library that implements the succinct data structure described in the article of Thomas C. Conway and Andrew J. Bromage.<br>
CBL is especially developped for assembling larges genomes. It stores k-mers and permits query operations. It represents an edge-centric De Bruijn Graph in a memory-efficient way by also supporting time-efficient operations.

Article : https://www.researchgate.net/publication/49765043_Succinct_Data_Structures_for_Assembling_Large_Genomes
## Method
ConwayBromageLib is based on the library SDSL (https://github.com/simongog/sdsl-lite) which provides succinct data structures. More precisely, CBL use the class '''sd_vector''', implementing a sparse bit vector (vector made of zeros and ones), to store the k-mers. <br>
<br>
CBL takes a list of canonical k-mers as input. The bit vector use to store these k-mers is of size 4^k. These latter are represented by an index/position in the bit vector. If an element of a specific position is set to one, then it means that the k-mer representing this position is present at least once in the sequence. Otherwise (set to 0) it is absent.<br>
<br>
A succinct data structure is a data structure  objects in memory-space close to the information-theoretic lower bounds and perform efficient query operations.

## What can I do with CBL ?
CBL stores genomes sequences. One it is done, we can apply some query operations on the bitvector : <br>
The query operation **contains** to know if a given k-mer (an element of size k, made of nucleotides) is present in the bitvector (the sequence) or not.<br>
The query operation **successors** to find out which successors of a given k-mer are present.<br> A successor of a k-mer x is any k-mer of the form a+x[1:k-1] or x[2:k]+a with a nucleotide.
## Requirements
CBL requires :<br>
``g++`` version 8.2.0 or higher.<br>
``cmake`` version 3.12 or higher (to build unit tests).
Mac OS and Linux (Ubuntu) are supported.
## Installation
To download the library, please use the following command : 
```
git clone --recurse-submodules https://github.com/yoann-dufresne/ConwayBromageLib.git
```
If you have already used ``git clone`` without the ``--recurse-submodules`` flag, please use in the project :
```
git submodule init
git submodule update
```
## How to use CBL ?
The current ``CMakeLists.txt`` **is made for unit tests**.<br>
To use CBL, please create a new C++ file.<br>
Example of code, main.cpp :<br>
```
#include "ConwayBromageLib.h"

using namespace std;

int main(){
    ifstream f("my_demo_test.txt"); //file which contains genome sequence
    KmerManipulatorACGT k(4);       // size of p-mers ((k+1)-mers)
    ConwayBromage cb(f, &k);        //build of the bitvector with 4-mers

    KmerManipulatorACGT k3(3);      //size of k-mers
    uint64_t my_mer = k3.encode("ACG"); //encode "ACG" into an unique number

    cb.contains(my_mer);    //return 1 if my_mer is present, otherwise 0
    cb.successors(my_mer);  //return a uint8_t that represent the 8 potential successors :
                            //1 = is present in the sequence, 0 = is absent.
    return 0;
}
```
To compile modify the ``CMakeLists.txt``, for example with main.cpp above : 
```
cmake_minimum_required(VERSION 3.12)    #The minimal version of cmake we need
                                        #Here we need 3.12 or higher version
project(example)               #The name of the projet, it will be the executable name
                               #Can be renamed, renamed it in all the CM akeLists.txt
add_compile_options(-Wall -Wextra)      #basic compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -mavx2 -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -O3 -DNDEBUG") #Compilation flag
#Variables which represents compilation elements
#One for .cpp elements
set(SRCS
        main.cpp  #Our main.cpp file example, can be modify according to your work
        ConwayBromageLib.cpp)
#One for .h
set(HEADERS
        ConwayBromageLib.h)
link_directories(example ./sdsl-lite/lib ./sdsl-lite/external/libdivsufsort/lib)    #Where we put .cpp from sdsl library
add_executable(example ${SRCS} ${HEADERS}) #Create executable example
add_subdirectory(external)      #Add subdirectories like sdsl and lest, call the CMakeLists.txt from external directory
target_include_directories(example PUBLIC ./external/lest/include)   #Where we put .h from lest libraries : It needs to be after add_executable !!!
target_link_libraries(example sdsl divsufsort divsufsort64)    #Specific flags from sdsl library

```
## Bug reporting
If there are bugs or mistakes in CBL please send us an issue on [CBL issue reports](https://github.com/yoann-dufresne/ConwayBromageLib/issues)
## External ressources
CBL uses others libraries to work.
The [SDSL](https://github.com/simongog/sdsl-lite) that is a library implementing succinct data structures. We use it to implement bitvector (sd_vector in SDSL syntax).
## KmerManipulator class
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

### bool ConwayBromage::contains(uint64_t Kmer)
Takes a (k-1)-mer and returns true if it's present among the stored k-mers.<br>
Example of use:<br>
```
bool kmer_exists = cb.contains(230);

KmerManipulatorACGT km(3);
uint64_t intGTT = km.decode("GTT");
bool GTT_exists = cb.contains(intGTT);
```

### uint8_t ConwayBromage::successors(uint64_t Kmer)
Returns the successors of a (k-1)-mer. A (k-1)-mer has at minimum 0 successors and maximum 8 successors.<br>
Informations : <br>
- The method doesn't check if the Kmer exists (so you have to do it on your own).<br>
- The result can sometimes have duplicates.<br>
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
