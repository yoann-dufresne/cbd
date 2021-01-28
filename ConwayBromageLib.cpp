#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include "ConwayBromageLib.h"
#include <immintrin.h>  //need -mavx2

using namespace std;
using namespace sdsl;

//abstract class KmerManipulator
KmerManipulator::KmerManipulator(uint64_t size): m_size(size) {}
KmerManipulator::~KmerManipulator() noexcept {}
//class KmerManipulatorACTG
/**
 * The constructor.
 * @param size - k-mer size. Must be inferior or equal to 32.
 */
KmerManipulatorACTG::KmerManipulatorACTG(uint64_t size): KmerManipulator(size), m_format("ACTG"){}
KmerManipulatorACTG::~KmerManipulatorACTG() noexcept {}

/** 
 * Encodes a k-mer in ACTG format.
 * @param word - the string we want to encode into a uint64_t.
 * @return a uint64_t which represent the encoding version of the word.
 */
uint64_t KmerManipulatorACTG::encode(const string &word) {
    //Value in ASCII table
    //A = 65 : 0100 0001 <-> 0 (encoded value)
    //C = 67 : 0100 0011 <-> 1  
    //T = 84 : 0101 0100 <-> 2 
    //G = 71 : 0100 0111 <-> 3 
    static const uint64_t caracterToValue[8] = {0, 0, 0, 1, 2, 0, 0, 3}; //declared and initialized only during the first call to this method
    //some 0 values in the array are never used because (word[i] & 0x7) will always produce one of the followind index : 1, 3, 4, or 7
    uint64_t res = 0;
    for(uint i = 0 ; i < m_size ; i++){   //We go through the sequence to encode each nucleotides
        res = (res << 2) + caracterToValue[word[i] & 0x7];  
    }
    return res;
}

/**
 * Returns a string representing the k-mer.
 * @param kmer - a value.
 * @return a string representing the decoding version of the uint64_t.
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * KmerManipulatorACTG decoder(4); //4-mers in ACTG format
 * decoder.decode(248);            //GGTA
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
string KmerManipulatorACTG::decode(uint64_t kmer) {
    static const char valueToCaracter[4] = {'A','C','T','G'}; //declared and initialized only during the first call to this method
    string res(m_size, ' ');
    uint64_t lastIndex = res.size()-1;
    for(int i = 0; i < m_size; i++){
        res[lastIndex-i] = valueToCaracter[kmer & 0x3];
        kmer >>= 2;
    }
    return res;
}

/**
 * Returns the canonical version of a k-mer in ACTG format.
 * @param kmer - a uint64_t which represents the compressed kmer we want to study.
 * @return an uint64_t which is the kmer itself if it is already canonical, or its canonical version.
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * KmerManipulatorACGT canonicaler(4); //4-mers in ACGT format
 * canonicaler.getCanonical(172);      //172
 * KmerManipulatorACTG canonicaler(4); //4-mers in ACTG format
 * canonicaler.getCanonical(172);      //144
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
uint64_t KmerManipulatorACTG::getCanonical(const uint64_t kmer) {
    uint64_t reverseCompl = reverseComplement(kmer);
    return((kmer < reverseCompl)?kmer:reverseCompl);
}

/** 
 * Returns the reverse complement.
 * @param kmer - an uin64_t which represent the compressed version of a k-mer.
 * @return a uin64_t which represent the compressed version of the reverse complement of the given k-mer.
 */
uint64_t KmerManipulatorACTG::reverseComplement(const u_int64_t kmer) {
    u_int64_t res = kmer;
    res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
    res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
    res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
    res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;
    return (res >> (2*( 32 - m_size))) ;
}

/**
 * Returns the reverse complement of a nucleotide.
 * @param nucleotide - an int between 0 and 3 (included)
 * @return the complement (either 0, 1, 2, or 3)
 */
uint8_t KmerManipulatorACTG::reverseComplementOfNucleotide(const uint8_t nucleotide){
    if(nucleotide == 0) //A
        return 2; //T
    if(nucleotide == 1) //C
        return 3; //G
    if(nucleotide == 2) //T
        return 0; //A
    // nucleotide = 3 <-> G
    return 1; //C
}

/**
 * Returns the corresponding caracter to the nucleotide's value.
 * @param nucleotide - an encode nucleotide
 * @return a char that correspond to the nucleotide value (A, C, T or G)
 */
char KmerManipulatorACTG::decodeNucleotide(const uint8_t nucleotide){
    if(nucleotide == 0) 
        return 'A';
    if(nucleotide == 1)
        return 'C';
    if(nucleotide == 2)
        return 'T';
    //nucleotide = 3
    return 'G';
}

/**
 * Returns the size attribute.
 * @return size k of k-mers that it is representing
 */
int KmerManipulatorACTG::getSize(){
    return m_size;
}

//class KmerManipulatorACGT
/**
 * The constructor.
 * @param size - k-mer size. Must be inferior or equal to 32.
 */
KmerManipulatorACGT::KmerManipulatorACGT(uint64_t size): KmerManipulator(size), m_format("ACGT") {}
KmerManipulatorACGT::~KmerManipulatorACGT() noexcept {}

/** 
 * Encodes a k-mer in ACGT format
 * @param word - the string we want to encode into a uint64_t
 * @return a uint64_t which represent the encoding version of the word
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * KmerManipulatorACGT encoder(4);  //4-mers in ACGT format
 * encoder.encode("GGTA");          //172
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
uint64_t KmerManipulatorACGT::encode(const string &word) {
    //value in ASCII table
    //A = 65 : 0100 0001 <-> 0 (encoded value)
    //C = 67 : 0100 0011 <-> 1   
    //G = 71 : 0100 0111 <-> 2 
    //T = 84 : 0101 0100 <-> 3
    static const int caracterToValue[8] = {0, 0, 0, 1, 3, 0, 0, 2}; //declared and initialized only during the first call to this method
    //some 0 values in the array are never used because (word[i] & 0x7) will always produce one of the followind index : 1, 3, 4, or 7
    uint64_t res = 0;
    for(uint i = 0 ; i < m_size ; i++){   //We go through the sequence to encode each nucleotides
        res = (res << 2) + caracterToValue[word[i] & 0x7];
    }
    return res;   
}

/**
 * Returns a string representing the k-mer.
 * @param kmer - a value.
 * @return a string representing the decoding version of the uint64_t.
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * KmerManipulatorACTG decoder(4); //4-mers in ACTG format
 * decoder.decode(248);            //GGTA
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
string KmerManipulatorACGT::decode(uint64_t kmer) {
    static const char valueToCaracter[4] = {'A','C','G','T'}; //declared and initialized only during the first call to this method
    string res(m_size, ' ');
    uint64_t lastIndex = res.size()-1;
    for(int i = 0; i < m_size; i++){
        res[lastIndex-i] = valueToCaracter[kmer & 0x3];
        kmer >>= 2;
    }
    return res;
}

/**
 * Returns the canonical version of a k-mer in ACGT format.
 * @param kmer - an uint64_t which represents the compressed kmer we want to study.
 * @return a uint64_t which is the compressed kmer itself if it is already canonical, or its canonical version.
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * KmerManipulatorACGT canonicaler(4); //4-mers in ACGT format
 * canonicaler.getCanonical(172);      //172
 * KmerManipulatorACTG canonicaler(4); //4-mers in ACTG format
 * canonicaler.getCanonical(172);      //144
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
uint64_t KmerManipulatorACGT::getCanonical(const uint64_t kmer) {
    uint64_t reverseCompl = reverseComplement(kmer);
    return((kmer < reverseCompl)?kmer:reverseCompl);
}

/** 
 * Returns the reverse complement.
 * @param kmer - a uin64_t which represent the compressed version of a k-mer.
 * @return a uin64_t which represent the compressed version of the reverse complement of the given k-mer.
 */
uint64_t KmerManipulatorACGT::reverseComplement(const uint64_t kmer) {
    uint64_t res = ~kmer;
    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    return (res >> (2 * (32 - m_size)));
}

/**
 * Returns the reverse complement of a nucleotide.
 * @param nucleotide - an int between 0 and 3 (included)
 * @return the complement (either 0, 1, 2, or 3)
 */
uint8_t KmerManipulatorACGT::reverseComplementOfNucleotide(const uint8_t nucleotide){
    if(nucleotide == 0) //A
        return 3; //T
    if(nucleotide == 1) //C
        return 2; //G
    if(nucleotide == 2) //G
        return 1; //C
    //nucleotide = 3 <-> T
    return 0; //A
}

/**
 * Returns the corresponding caracter to the nucleotide's value.
 * @param nucleotide - an encode nucleotide
 * @return a char that correspond to the nucleotide value (A, C, T or G)
 */
char KmerManipulatorACGT::decodeNucleotide(const uint8_t nucleotide){
    if(nucleotide == 0)
        return 'A';
    if(nucleotide == 1)
        return 'C';
    if(nucleotide == 2)
        return 'G';
    // nucleotide = 3
    return 'T';
}

/**
 * Returns the size k of the k-mers that the object is representing.
 * @return k-mers length
 */
int KmerManipulatorACGT::getSize(){
    return m_size;
}

/**
 * First ConwayBromage constructor : Stores memory-efficiently k-mers coming from an istream in a specific format.  
 * @param kmerFlux - An istream of k-mers. For example a file represented by an ifstream.
 * @param km - Object which contains all the information about encoding and decoding in a specific format.
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * int kmerSize = 31;
 * ifstream f("./k-mers.txt", ios::in);
 * KmerManipulatorACGT km(kmerSize);
 * ConwayBromage cb(f, &km);        
 * f.close();
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * @warning The k-mers in the istream must me canonical and have a size k <= 32.
 * @warning The istream must have at each line (if it's a file for example) a unique k-mers.
 * @warning The k-mers in the istream must be sorted either in lexicographical order (A<C<G<T) or in A<C<T<G. It depends on the encoding format.
 */
ConwayBromage::ConwayBromage(istream& kmerFlux, KmerManipulator* km){
    m_kmerManipulator = km;

    string line("");
    uint64_t numberOfKmer = 0;
    while(getline(kmerFlux, line)) //Counts the number of k-mer in the file
        numberOfKmer++;
 
    kmerFlux.clear();
    kmerFlux.seekg(0, ios::beg);    //Return to the beginning of the file
    uint64_t one = 1;
    uint64_t sdvSize = one << ((2*m_kmerManipulator->getSize())-1); //Creation of the total length to create the sd_vector_builder
    sd_vector_builder builder(sdvSize, numberOfKmer);
    uint64_t previousKmer(0);
    while(getline(kmerFlux, line)){
        uint64_t k = m_kmerManipulator->encode(line);
        //first chack : canonical elements expected
        if(k >  m_kmerManipulator->reverseComplement(k)){
            cout << "The file is not completely canonical" << endl;
            exit(1);
        }
        //second check : ascending order expected
        if(k < previousKmer){
            cout << "The file is not sort in the ascending order" << endl;
            exit(1);
        }
        builder.set(k);   
        previousKmer = k;
    }
    
    m_sequence = builder;
    m_limit = (m_sequence.size() >> 1);
    
    //initialization of the cache
    int nbOfBitsToShift = 2 * (m_kmerManipulator->getSize()-1);
    m_numberOfBitsToShift = nbOfBitsToShift;

    m_RC[0] = km->reverseComplementOfNucleotide(0);
    m_RC[1] = km->reverseComplementOfNucleotide(1);
    m_RC[2] = km->reverseComplementOfNucleotide(2);
    m_RC[3] = km->reverseComplementOfNucleotide(3);

    m_RC_shifted[0] = m_RC[0] << nbOfBitsToShift;
    m_RC_shifted[1] = m_RC[1] << nbOfBitsToShift;
    m_RC_shifted[2] = m_RC[2] << nbOfBitsToShift;
    m_RC_shifted[3] = m_RC[3] << nbOfBitsToShift;
    
    m_nucleotides_shifted[0] = 0 << nbOfBitsToShift;
    m_nucleotides_shifted[1] = 1 << nbOfBitsToShift;
    m_nucleotides_shifted[2] = 2 << nbOfBitsToShift;
    m_nucleotides_shifted[3] = 3 << nbOfBitsToShift;

    for(int i = 0; i < 4; i++){
        //for next and previous pmers
        char nucleotideOfI = m_kmerManipulator->decodeNucleotide(i);
        if(nucleotideOfI == 'A'){ 
            m_correspondingBitValueForNextPmers[i] = 128;
            m_correspondingBitValueForPrevPmers[i] = 8;  
        } else if(nucleotideOfI == 'C') {     
            m_correspondingBitValueForNextPmers[i] = 64; 
            m_correspondingBitValueForPrevPmers[i] = 4; 
        } else if(nucleotideOfI == 'G') {
            m_correspondingBitValueForNextPmers[i] = 32; 
            m_correspondingBitValueForPrevPmers[i] = 2; 
        } else {//'T'
            m_correspondingBitValueForNextPmers[i] = 16; 
            m_correspondingBitValueForPrevPmers[i] = 1; 
        }
    }
}

/**
 * Second ConwayBromage constructor : Build the ConwayBromage object based on the sd_vector in parameter and the KmerManipulator.
 * @param sdv - An sd_vector which represents the sequence.
 * @param km - A KmerManipulator.
 */
ConwayBromage::ConwayBromage(sdsl::sd_vector<> const& sdv, KmerManipulator* km){
    m_sequence = sdv; //copy of the sd_vector
    m_kmerManipulator = km;
    m_limit = (m_sequence.size() >> 2) - 1;
    
    //initialization of the cache
    int nbOfBitsToShift = 2 * (m_kmerManipulator->getSize()-1);
    m_numberOfBitsToShift = nbOfBitsToShift;

    m_RC[0] = km->reverseComplementOfNucleotide(0);
    m_RC[1] = km->reverseComplementOfNucleotide(1);
    m_RC[2] = km->reverseComplementOfNucleotide(2);
    m_RC[3] = km->reverseComplementOfNucleotide(3);

    m_RC_shifted[0] = m_RC[0] << nbOfBitsToShift;
    m_RC_shifted[1] = m_RC[1] << nbOfBitsToShift;
    m_RC_shifted[2] = m_RC[2] << nbOfBitsToShift;
    m_RC_shifted[3] = m_RC[3] << nbOfBitsToShift;
    
    m_nucleotides_shifted[0] = 0 << nbOfBitsToShift;
    m_nucleotides_shifted[1] = 1 << nbOfBitsToShift;
    m_nucleotides_shifted[2] = 2 << nbOfBitsToShift;
    m_nucleotides_shifted[3] = 3 << nbOfBitsToShift;

    for(int i = 0; i < 4; i++){
        //for next and previous pmers
        char nucleotideOfI = m_kmerManipulator->decodeNucleotide(i);
        if(nucleotideOfI == 'A'){ 
            m_correspondingBitValueForNextPmers[i] = 128;
            m_correspondingBitValueForPrevPmers[i] = 8;  
        } else if(nucleotideOfI == 'C') {     
            m_correspondingBitValueForNextPmers[i] = 64; 
            m_correspondingBitValueForPrevPmers[i] = 4; 
        } else if(nucleotideOfI == 'G') {
            m_correspondingBitValueForNextPmers[i] = 32; 
            m_correspondingBitValueForPrevPmers[i] = 2; 
        } else {//'T'
            m_correspondingBitValueForNextPmers[i] = 16; 
            m_correspondingBitValueForPrevPmers[i] = 1; 
        }
    }
}

/**
 * Check if the given (k-1)-mer is present. The (k-1)-mer can either be canonical or not.
 * @param Kmer - An uint64_t representing the (k-1)-mer.
 * @return True if the (k-1)-mer is present. False otherwise. 
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * bool kmer_exists = cb.contains(230);
 *
 * KmerManipulatorACGT km(3);
 * uint64_t intGTT = km.encode("GTT");
 * bool GTT_exists = cb.contains(intGTT);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
bool ConwayBromage::contains(uint64_t Kmer) const{
    if(Kmer > m_limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << m_limit << endl;
        return false; 
    }  
    
    uint64_t PmerPrev, RC_PmerNext, RC_PmerPrev;
    uint64_t PmerNext = Kmer << 2; //<-> XA where X is the Kmer
    uint64_t RC_Kmer = m_kmerManipulator->reverseComplement(Kmer) >> 2;
    uint64_t RC_Kmer_ShiftedOf2Bits = RC_Kmer << 2;

    for(int i = 0; i < 4; i++){ 
        //next
        RC_PmerNext = RC_Kmer + m_RC_shifted[i]; //YX where X is the reverse complement of the kmer and Y is the reverse complement of the last nucleotide of the pmer
        if(m_sequence[(PmerNext < RC_PmerNext)?PmerNext:RC_PmerNext]) return true;
        //previous
        PmerPrev = Kmer + m_nucleotides_shifted[i]; //equals to AX then CX then GX then TX where X is the Kmer
        RC_PmerPrev = RC_Kmer_ShiftedOf2Bits + m_RC[i];
        if(m_sequence[(PmerPrev < RC_PmerPrev)?PmerPrev:RC_PmerPrev]) return true;
        
        PmerNext++; //equals to XA then XC then XG then XT
    }
    return false;
}

/**
 * Returns the successors of a (k-1)-mer. The (k-1)-mer can either be canonical or not. A (k-1)-mer has at minimum 0 successors and maximum 8 successors.
 * ### How it works ? 
 * Let's say the function takes as a parameter the integer which represents the (k-1)-mer GTT and we assume we are in ACGT encoding.
 * First, the method will generate the following k-mers : GTTA, GTTC, GTTG, GTTT, AGTT, CGTT, GGTT, TGTT.
 * Then, it will check if the canonical version of these k-mers are present.
 * Let's suppose that only the canonical version of these k-mers are present : GTTA, GTTT, CGTT, TGTT.
 * GTTA and GTTT are, what we call next k-mers so we will store the information in the 4 left bits of the result (uint8_t).
 * CGTT and TGTT are previous k-mers so, this time, the information will be stored in the 4 right bits of the result.
 * Thus, the function will return 1001 0101 which corresponds to 149 in base 10.
 * It means that, in this example, the successors of the (k-1)-mer GTT are TTA, TTT, CGT and TGT.
 *
 * @param Kmer : a (k-1)-mer.
 * @return a uint8_t which carry information about the presence/absence of the 8 potential successors of the (k-1)-mer.
 * Read of the return : 
 * read the uint8_t under bit form : the first four elements are "next" k-mers, the last four are "previous"
 * Always the same reading direction : A then C then G then T regardless of the encoding
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * uint8_t successors = cb.successors(78);
 * KmerManipulatorACGT km(3);
 * uint64_t intGTT = km.encode("GTT");
 * uint8_t successorsOfGTT = cb.successors(intGTT);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * @warning The method doesn't check the presence of the (k-1)-mer (so you have to do it on your own with the method contains for example).
 */
uint8_t ConwayBromage::successors(uint64_t Kmer) const{
    if(Kmer > m_limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be equal or inferior to 4^(P-1) i.e " << m_limit << endl;
        return 0; //empty
    }
    uint8_t res = 0;
    //we build the eight possible successors    
    uint64_t PmerPrev, RC_PmerNext, RC_PmerPrev;
    uint64_t PmerNext = Kmer << 2; //<-> XA where X is the Kmer
    uint64_t RC_Kmer = m_kmerManipulator->reverseComplement(Kmer) >> 2;
    uint64_t RC_Kmer_ShiftedOf2Bits = RC_Kmer << 2;

    for(int i = 0; i < 4; i++){ 
        RC_PmerNext = RC_Kmer + m_RC_shifted[i]; //YX where X is the reverse complement of the kmer and Y is the reverse complement of the last nucleotide of the pmer
        if(m_sequence[(PmerNext < RC_PmerNext)?PmerNext:RC_PmerNext]){ //if the pmer is present 
            res ^= m_correspondingBitValueForNextPmers[i];
        }
        PmerPrev = Kmer + m_nucleotides_shifted[i]; //equals to AX then CX then GX then TX where X is the Kmer
        RC_PmerPrev = RC_Kmer_ShiftedOf2Bits + m_RC[i];
        if(m_sequence[(PmerPrev < RC_PmerPrev)?PmerPrev:RC_PmerPrev]){
            res ^= m_correspondingBitValueForPrevPmers[i];
        }
        PmerNext++; //equals to XA then XC then XG then XT
    }
    
    return res;
}

/**
 * Returns the size k of the k-mer that are stored.
 * @return an int
 */
int ConwayBromage::getKmerSize(){
    return m_kmerManipulator->getSize();
}

/**
 * Returns the compressed sequence.
 * @return an sd_vector
 */
sdsl::sd_vector<> ConwayBromage::getSequence(){
    return m_sequence;
}

/**
 * Returns the kmer manipulator.
 * @return a KmerManipulator object.
 */
KmerManipulator* ConwayBromage::getKmerManipulator(){
    return m_kmerManipulator;
}


