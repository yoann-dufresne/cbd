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
//For KmerManipulator and children comments, please tcheck the ConwayBromageLib.h
KmerManipulator::KmerManipulator(uint64_t size): m_size(size) {}
KmerManipulator::~KmerManipulator() noexcept {}
//class KmerManipulatorACTG
KmerManipulatorACTG::KmerManipulatorACTG(uint64_t size): KmerManipulator(size), m_format("ACTG"){}
KmerManipulatorACTG::~KmerManipulatorACTG() noexcept {}

/** 
 * Encodes a k-mer in ACTG format
 * @param word - the string we want to encode into a uint64_t
 * @return a uint64_t which represent the encoding version of the word
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
 * Returns a string which is the word which corresponds to the uint64_t value
 * @param kmer - the value.
 * @return a string representing the decoding version of the uint64_t.
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
 * Returns the canonical version of a  k-mer, depends on encoding : see daughter classes
 * @param kmer - a uint64_t which represents the compressed kmer we want to study
 * @return a uint64_t which is the kmer itself if it is already canonical, or its canonical version
 */
uint64_t KmerManipulatorACTG::getCanonical(const uint64_t kmer) {
    uint64_t reverseCompl = reverseComplement(kmer);
    return((kmer < reverseCompl)?kmer:reverseCompl);
}

//AVX version of getCanonical
__m256i KmerManipulatorACTG::getCanonicalAVX(const __m256i kmer) { 
    __m256i reverseCompl = reverseComplementAVX(kmer);  //creation of a __m256i which contains reverse complement of the __m256i k-mers
    __m256i eq = _mm256_cmpeq_epi64(kmer,reverseCompl); //Case where forward and reverse are the same
    __m256i cmp1 = _mm256_or_si256(_mm256_cmpgt_epi64(kmer, reverseCompl), eq);  //kmer < reverseCompl => 0, else -1
    __m256i cmp2 = _mm256_cmpgt_epi64(reverseCompl, kmer);  //reverseCompl < kmer => 0, else -1
    __m256i and1 = _mm256_and_si256(cmp2, kmer);    //delete elements wich are higher than their reverse
    __m256i and2 = _mm256_and_si256(cmp1, reverseCompl);    //delete reverses which are higher than their forward
    __m256i ult = _mm256_or_si256(and1, and2);  //final construction
    return ult;
}

/** 
 * Returns the reverse complement
 * @param kmer - a uin64_t which represent the compressed version of a k-mer
 * @return a uin64_t which represent the compressed version of the reverse complement of the given k-mer
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

//AVX version of reverseComplement
__m256i KmerManipulatorACTG::reverseComplementAVX(const __m256i kmer){
    //We initialize variables like __mm256i vectors
    __m256i AAA = _mm256_set1_epi64x(0xAAAAAAAAAAAAAAAA);
    __m256i val1 = _mm256_set1_epi64x(0x3333333333333333);
    __m256i val2 = _mm256_set1_epi64x(0x0F0F0F0F0F0F0F0F);
    __m256i val3 = _mm256_set1_epi64x(0x00FF00FF00FF00FF);
    __m256i val4 = _mm256_set1_epi64x(0x0000FFFF0000FFFF);
    __m256i val5 = _mm256_set1_epi64x(0x00000000FFFFFFFF);
    __m256i res = kmer;  // <=> res = kmer
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 2), val1), _mm256_slli_epi64(_mm256_and_si256(res, val1), 2))); // <=> res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 4), val2), _mm256_slli_epi64(_mm256_and_si256(res, val2), 4)));
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 8), val3), _mm256_slli_epi64(_mm256_and_si256(res, val3), 8)));
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 16), val4), _mm256_slli_epi64(_mm256_and_si256(res, val4), 16)));
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 32), val5), _mm256_slli_epi64(_mm256_and_si256(res, val5), 32)));
    res = _mm256_xor_si256(res, AAA);
    res = _mm256_srli_epi64(res, (2 * (32 - m_size)));  //(res >> (2 * (32 - m_size))), we can return it directly
    return res;
}

/**
 * Returns the reverse complement of a nucleotide.
 * @param nucleotide - an encode nucleotide
 * @return the reverse complement of the param nucleotide
 * We use this method for cache creation
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
 * @return an int that is the size of the sequence
 */
int KmerManipulatorACTG::getSize(){
    return m_size;
}

//class KmerManipulatorACGT
KmerManipulatorACGT::KmerManipulatorACGT(uint64_t size): KmerManipulator(size), m_format("ACGT") {}
KmerManipulatorACGT::~KmerManipulatorACGT() noexcept {}

/** 
 * Encodes a k-mer in ACGT format
 * @param word - the string we want to encode into a uint64_t
 * @return a uint64_t which represent the encoding version of the word
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
 * Returns a string which is the word which corresponds to the uint64_t value (in ACGT format)
 * @param kmer - the value.
 * @return a string representing the decoding version of the uint64_t.
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
 * Returns the canonical version of a k-mer in ACGT format
 * @param kmer - a uint64_t which represents the compressed kmer we want to study
 * @return a uint64_t which is the compressed kmer itself if it is already canonical, or its canonical version
 */
uint64_t KmerManipulatorACGT::getCanonical(const uint64_t kmer) {
    uint64_t reverseCompl = reverseComplement(kmer);
    return((kmer < reverseCompl)?kmer:reverseCompl);
}

// AVX version of getCanonical
__m256i KmerManipulatorACGT::getCanonicalAVX(const __m256i kmer) {
    __m256i reverseCompl = reverseComplementAVX(kmer);  //creation of a __m256i which contains reverse complement of the __m256i k-mers
    __m256i eq = _mm256_cmpeq_epi64(kmer,reverseCompl); //Case where reverse and forward are the same
    __m256i cmp1 = _mm256_or_si256(_mm256_cmpgt_epi64(kmer, reverseCompl), eq);  //kmer < reverseCompl => 0, else -1
    __m256i cmp2 = _mm256_cmpgt_epi64(reverseCompl, kmer);  //reverseCompl < kmer => 0, else -1
    __m256i and1 = _mm256_and_si256(cmp2, kmer);    //delete elements wich are higher than their reverse
    __m256i and2 = _mm256_and_si256(cmp1, reverseCompl);    //delete reverses which are higher than their forward
    __m256i ult = _mm256_or_si256(and1, and2);  //final construction
    return ult;
}

/** 
 * Returns the reverse complement
 * @param kmer - a uin64_t which represent the compressed version of a k-mer
 * @return a uin64_t which represent the compressed version of the reverse complement of the given k-mer
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

//AVX version of reverseComplement
__m256i KmerManipulatorACGT::reverseComplementAVX(const __m256i kmer) {
    //We initialize variables like __mm256i vectors
    __m256i allAtOne = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
    __m256i val1 = _mm256_set1_epi64x(0x3333333333333333);
    __m256i val2 = _mm256_set1_epi64x(0x0F0F0F0F0F0F0F0F);
    __m256i val3 = _mm256_set1_epi64x(0x00FF00FF00FF00FF);
    __m256i val4 = _mm256_set1_epi64x(0x0000FFFF0000FFFF);
    __m256i val5 = _mm256_set1_epi64x(0x00000000FFFFFFFF);
    __m256i res = _mm256_andnot_si256(kmer, allAtOne);  // <=> res = ~kmer
    //translation in AVX instruction of the original reverseComplement : CARE ABOUT UNSIGNED INT !!!
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 2), val1), _mm256_slli_epi64(_mm256_and_si256(res, val1), 2))); // <=> res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 4), val2), _mm256_slli_epi64(_mm256_and_si256(res, val2), 4)));
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 8), val3), _mm256_slli_epi64(_mm256_and_si256(res, val3), 8)));
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 16), val4), _mm256_slli_epi64(_mm256_and_si256(res, val4), 16)));
    res = (_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi64(res, 32), val5), _mm256_slli_epi64(_mm256_and_si256(res, val5), 32)));
    res = _mm256_srli_epi64(res, (2 * (32 - m_size)));  //(res >> (2 * (32 - m_size))), we can return it directly
    return res;
}

/**
 * Returns the reverse complement of a nucleotide.
 * @param nucleotide - an encode nucleotide
 * @return the reverse complement of the param nucleotide
 * We use this method for cache creation
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
 * Returns the size attribute.
 * @return an int that is the size of the sequence
 */
int KmerManipulatorACGT::getSize(){
    return m_size;
}

/**
 * First ConwayBromage constructor : Build the ConwayBromage object based on the sd_vector in parameter and the KmerManipulator.
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
 * Second ConwayBromage constructor : Transform sequences which are contain in a file in a sd_vector
 * @param kmerFlux - An istream of kmer. For example a file represented by an ifstream.
 * @param km - Object which contains all the information about encoding and decoding.
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
    uint64_t sdvSize = one << (2*m_kmerManipulator->getSize()); //Creation of the total length to create the sd_vector_builder

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
 * Return the value stored at the block of index i.
 * @param i - the index.
 * @return true if the block of index i is set to one.
 */
bool ConwayBromage::operator[](uint64_t i) const{
    return m_sequence[i];
}
//AVX version of the [] operator
bool ConwayBromage::operator[](__m256i i) const {
    if(m_sequence[i[0]]) return true;
    if(m_sequence[i[1]]) return true;
    if(m_sequence[i[2]]) return true;
    if(m_sequence[i[3]]) return true;
    return false;
}

/**
 * Return the size of the compressed sequence. It's always 4^K with K the size of the K-mers 
 * @return the size
 */
uint64_t ConwayBromage::size() const{
    return m_sequence.size();
}

/**
 * Check if the given (k-1)-mer is present. The (k-1)-mer can either be canonical or not.
 * @param Kmer - An uint64_t representing the (k-1)-mer.
 * @return true if the (k-1)-mer is present. False otherwise. 
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

//AVX version of isPresent
bool ConwayBromage::isPresentAVX(uint64_t Kmer) const {
    int KmerSize = m_kmerManipulator->getSize()-1;
    uint64_t next = Kmer << 2;
    int numberOfBitsToShift = KmerSize << 1;
    /*
     * Creation of 2 __m256i vectors of 4 elements (can't create an higher one)
     * _mm256_setr_epi64x is faster than _mm256_set_epi64x
     */
    if(operator[](m_kmerManipulator->getCanonicalAVX(_mm256_setr_epi64x(next, next+1, next+2, next+3)))) return true;
    if(operator[](m_kmerManipulator->getCanonicalAVX(_mm256_setr_epi64x(Kmer, (Kmer + (1 << numberOfBitsToShift)),
                                                                        (Kmer + (2 << numberOfBitsToShift)), (Kmer + (3 << numberOfBitsToShift)))))) return true;
    return false;
}

/**
 * Returns the successors of a (k-1)-mer. The (k-1)-mer can either be canonical or not.
 * The function doesn't check if the (k-1)-mer is present.
 * @param Kmer : a (k-1)-mer.
 * @return a uint8_t which carry information about the presence/absence of the 8 potential successors of the (k-1)-mer.
 * Read of the return : 
 * read the uint8_t under bit form : the first four elements are "next" k-mers, the last four are "previous"
 * Always the same reading direction : A then C then G then T regardless of the encoding
 * Check the read.me "Explanation on how it works" in https://github.com/yoann-dufresne/ConwayBromageLib for an example
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
 * Return the size of the sequence.
 * @return an int
 */
int ConwayBromage::getKmerSize(){
    return m_kmerManipulator->getSize();
}

/**
 * Return the compressed sequence.
 * @return an sd_vector
 */
sdsl::sd_vector<> ConwayBromage::getSequence(){
    return m_sequence;
}

/**
 * Return the kmer manipulator.
 * @return a KmerManipulator object
 */
KmerManipulator* ConwayBromage::getKmerManipulator(){
    return m_kmerManipulator;
}


