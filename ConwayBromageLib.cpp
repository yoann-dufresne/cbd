//
// Created by Alexandra and Murat on 05/06/2020.
//

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

/* Successor counter to verify veracity of successors function
 * It is a non optimized function (slow), we use it for tests exclusively
 * @param compressedKmer - The k-mer for which we want to know successors
 * @param currentCompressedSeq - The given which is supposed to contain successors of compressedKer
 * @param format - a boolean to konw if we use ACGT or ACTG encoding
 * @return a vector which contain successors of compressedKMer
 */
/*vector<uint64_t> successorCounter(uint64_t compressedKMer, sd_vector<>currentCompressedSeq, bool format){     //Need upload to POO version
    uint64_t sizeOfSeq = log(currentCompressedSeq.size()) / log(ALPHABET);
    uint64_t sizeOfKMer = sizeOfSeq - 1;
    uint64_t list[8];
    string istr;
    vector<uint64_t> reader;
    if(format){
        istr = decode(compressedKMer, sizeOfKMer);
        list[0] = encode(istr + "A", sizeOfSeq);
        list[1] = encode(istr + "C", sizeOfSeq);
        list[2] = encode(istr + "G", sizeOfSeq);
        list[3] = encode(istr + "T", sizeOfSeq);
        list[4] = encode("A" + istr, sizeOfSeq);
        list[5] = encode("C" + istr, sizeOfSeq);
        list[6] = encode("G" + istr, sizeOfSeq);
        list[7] = encode("T" + istr, sizeOfSeq);
    }else{
        istr = decodeEcoli(compressedKMer, sizeOfKMer);
        list[0] = encodeEcoli(istr + "A", sizeOfSeq);
        list[1] = encodeEcoli(istr + "C", sizeOfSeq);
        list[2] = encodeEcoli(istr + "G", sizeOfSeq);
        list[3] = encodeEcoli(istr + "T", sizeOfSeq);
        list[4] = encodeEcoli("A" + istr, sizeOfSeq);
        list[5] = encodeEcoli("C" + istr, sizeOfSeq);
        list[6] = encodeEcoli("G" + istr, sizeOfSeq);
        list[7] = encodeEcoli("T" + istr, sizeOfSeq);
    }
    for(int i = 0 ; i < 8 ; i++){
        uint64_t lecteur;
        if(!isCanonical(list[i], sizeOfSeq, format)){
            //cout << "not cano : " << decodeEcoli(list[i], 3) << endl;
            uint64_t canoVers = getCanonical(list[i], sizeOfSeq, format);
            //cout << "cano version : " << canoVers << endl;
            if(currentCompressedSeq[canoVers]){
                //cout << "PASSED" << endl;
                if(format){
                    if(i > 3){
                        lecteur = encode(decode(list[i], sizeOfSeq).erase(sizeOfSeq-1, 1), sizeOfKMer);
                    }else{
                        lecteur = encode(decode(list[i], sizeOfSeq).erase(0, 1), sizeOfKMer);
                    }
                }else{
                    if(i > 3){
                        lecteur = encodeEcoli(decodeEcoli(list[i], sizeOfSeq).erase(sizeOfSeq-1, 1), sizeOfKMer);
                    }else{
                        lecteur = encodeEcoli(decodeEcoli(list[i], sizeOfSeq).erase(0, 1), sizeOfKMer);
                    }
                }
                //cout << "lecteur : " << lecteur << endl;
                if(find(reader.begin(), reader.end(), lecteur) == reader.end()){
                   reader.push_back(lecteur);
               }
            }
        }else{
            if(currentCompressedSeq[list[i]]){
                if(format){
                    if(i > 3){
                        lecteur = encode(decode(list[i], sizeOfSeq).erase(sizeOfSeq-1, 1), sizeOfKMer);
                    }else{
                        lecteur = encode(decode(list[i], sizeOfSeq).erase(0, 1), sizeOfKMer);
                    }
                }else{
                    if(i > 3){
                        lecteur = encodeEcoli(decodeEcoli(list[i], sizeOfSeq).erase(sizeOfSeq-1, 1), sizeOfKMer);
                    }else{
                        lecteur = encodeEcoli(decodeEcoli(list[i], sizeOfSeq).erase(0, 1), sizeOfKMer);
                    }
                }
                //cout << "lecteur : "<< lecteur << endl;
                if(find(reader.begin(), reader.end(), lecteur) == reader.end()){
                   reader.push_back(lecteur);
               }
            }
        }
    }
    return reader;
}*/

/* Variation of getCanonical : just tell if it is canonical or not
 * Just use for tests
 * @param kmer - the k-mer we want to test
 * @param kmerSize - the size if this k-mer
 * @param encodingIsACGT - boolean to know if encoding is ACGT or ACTG
 * @return true if it is canonical, else false
 */
/*bool isCanonical (uint64_t kmer, uint64_t kmerSize, bool encodingIsACGT){     //Need upload to POO version
    if(encodingIsACGT){ //ACGT encoding
        uint64_t reverseComplement = reverseComplementLexico(kmer, kmerSize);
        return (kmer <= reverseComplement)?true:false;
    }
    //ACTG encoding
    uint64_t reverseComplement = reverseComplementGATBLibEcoli(kmer, kmerSize);
    return (kmer <= reverseComplement)?true:false;
}*/

/* Translate an 8 bits number into a succession of uint64_t which represent compressed successors of the compressed K-mer
 * @param successors - a int which is a number between 0 and 255 and represent successors   -> maybe take uint8_t ?
 * @param compressedKMer - A compress version of the k-mer for which we have calculate successors
 * @param size - The size of the compressed k-mer
 * @param format - Give the encoding format : true = ACGT, false = ACTG
 * @return a vector of uint64_t which represent compressed successors of the compressed kmer
 */
vector<uint64_t> successorTranslator(int successors, uint64_t compressedKMer, uint64_t size, bool format){    //need upload to POO version
    vector<uint64_t>compressedSucc; //for the return
    if(successors > 256){
        cout << "bad 8 bits numbers" << endl;
    }else{
        int binaryList[8]{0};   // successors is on 8 bits, each cases will contain a bit of successors
        uint64_t binaryForm(1); //element for the int to binary transformation
        //cout << successors << endl;
        int succVal(successors); //copy of successors
        while(succVal > 0){ //from int to binary transformation begining
            for(int i = 0 ; i <= 8 ; i++){  //8 bits -> from 0 to 255 <=> 0 to (2^8)-1
                //binaryForm is like : 2^0, 2^1 till 2^8
                if(binaryForm > succVal){   // successors is not a 2^n form, with n in[0, 8]
                    //if 2^n > succVal then succVal > 2^(n-1) (case equal below)
                    succVal -= (binaryForm >> 1);   //The new succVal   binary >> 1 <=> 2^(n-1)
                                                    //example : 16 >> 1 = 8 -> 2^4 >> 1 = 2^3
                    binaryList[8-i] = 1;    //The case is set at 1
                    break;  //incrementation of i and begin again until succVal = 0
                }else if(binaryForm == succVal){    //case equal
                    succVal -= binaryForm;      //if it is equal, succVal will be zero
                    binaryList[8-i-1] = 1;      //Start from 0
                    break;  //succVal is 0, we will quit the for loop
                }
                binaryForm = binaryForm << 1;   // <=> 2^n << 1 = 2^(n+1) then i incrementation
            }
            binaryForm = 1; //we have quit the for loop, if succVal != O we need to restart it until it is set at 0
        }
        /*cout << "binary form :" << endl;
        for(int i = 0 ; i < 8 ; i++){
            cout << binaryList[i] << " ";
        }
        cout << endl;*/
        uint64_t limit = pow(ALPHABET, size);   //The limit we can't be taller than. Example : 3-mer limit is 64
        uint64_t numberOfBitToShift((size-1) << 1); //The number of bits we want to shift if there is previous among successors
        /* Explanation :
         * if we have AACA, previous successors will have a XAAC form and X can be A, C; T or G.
         * Each letters cost 2 bits. So X <- +2bits -- A <- +2bits -- A <- +2bits -- C (C is at zero position, don't need to count)
         * So, I need to shift 6 bits to go to X and modify it. Here AACA is a 4-mers, so I need (4-1) * 2 = 6 bits.
         * A faster writing is : (4-1) << 1 = 3 << 1 <=> 3*2 = 6
         */
        /*cout << "number of bit to shift : " << numberOfBitToShift << endl;
        cout << "limit : " << limit << endl;*/
        for(int i = 0 ; i < 8 ; i++){   //one case per bit
            if(binaryList[i] == 1){ //The successor is present
                if(i < 4){  //First fourth (from 0 to 3) are next
                    compressedSucc.push_back(((compressedKMer << 2)%limit) + ((!format)?i:((i%2 == 0)?i+1:i-1)));
                    /* Same explanation but we want ACAX form
                    * AACA << 2 <=> AACAA and we quit the limit of a 4-mers.
                    * In this case, we have to take care of the limit (256 for 4-mers) and have to stay in it
                    * That is why we use the rest of the divide by limit
                    * Finally we obtain ACAA which is the first form
                    * We don't need to shit to the left because we just have to add the letter at the end (place 0)
                    */
                }else{ // from 4 to 7 are previous
                    int nexI = i-4;
                    compressedSucc.push_back((compressedKMer >> 2) + (((!format)?nexI:((nexI%2 == 0)?nexI+1:nexI-1)) << numberOfBitToShift));
                    /* Explanation : 2 steps, AACA example
                     * I know it is previous, I need a XAAC form and each letters cost 2 bits
                     * So AACA >> 2 <=> AAAC, we have the first form, let's check the others 3 with addition
                     * Depends on format = ACGT : A=0, C=1, G=2, T=3 ; ACTG : A=1, C=0, T=3, G=2 (according to decoders)
                     * When we have the correct number, we shift it at the begin (where the
                     * X is) thanks to numberOfBitToShift (see it above)
                     */
                }
            }
        }
        /*cout << "final :"  << compressedSucc << endl;
        for(int i = 0 ; i < compressedSucc.size() ; i++){
            cout << compressedSucc[i] << " : " << decodeEcoli(compressedSucc[i], size) << endl;
        }*/
    }
    return compressedSucc;
}

//POO for KmerManipulator
//abstract class KmerManipulator
KmerManipulator::KmerManipulator(uint64_t size): m_size(size) {}
KmerManipulator::~KmerManipulator() noexcept {}
//class KmerManipulatorACTG
KmerManipulatorACTG::KmerManipulatorACTG(uint64_t size): KmerManipulator(size), m_format("ACTG"){}
KmerManipulatorACTG::~KmerManipulatorACTG() noexcept {}
uint64_t KmerManipulatorACTG::encode(const string word) {
    uint64_t hash = 0;
    char c;
    for(uint i = 0 ; i < m_size ; i++){   //We go through the sequence to encode each nucleotides
        hash <<= 2; // We shift 2 bits to the right
        c = word[i];    //Take each nucleotides to encode them
        uint64_t charval = 0; // 'A' = 0
        if(c == 'C'){
            charval = 1;    // 'C' = 1
        }else if(c == 'T'){
            charval = 2;    // 'T' = 2
        }else if(c == 'G'){
            charval = 3;    // 'G' = 3
        }
        hash += charval;    //creation of the hash for the given sequence
    }
    return hash;    //return the final hash of the sequence
}
string KmerManipulatorACTG::decode(uint64_t kmer) {
    string res(m_size, ' ');
    uint64_t lastIndex = res.size()-1;
    for(int i(0); i < m_size; i++){
        switch(kmer & 0x3){ //compares the decimal value of the first two bits
            case 0: res[lastIndex-i] = 'A'; break;
            case 1: res[lastIndex-i] = 'C'; break;
            case 2: res[lastIndex-i] = 'T'; break;
            case 3: res[lastIndex-i] = 'G'; break;
        }
        kmer >>= 2;
    }
    return res;
}
uint64_t KmerManipulatorACTG::getCanonical(const uint64_t kmer) {
    uint64_t reverseCompl = reverseComplement(kmer);
    return((kmer < reverseCompl)?kmer:reverseCompl);
}

__m256i KmerManipulatorACTG::getCanonicalAVX(const __m256i kmer) {  //WORK IN PROGRESS
    return _mm256_set1_epi64x(0);
}

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

__m256i KmerManipulatorACTG::reverseComplementAVX(const __m256i kmer){  //WORK IN PROGRESS
    return _mm256_set1_epi64x(0);
}
//class KmerManipulatorACGT
KmerManipulatorACGT::KmerManipulatorACGT(uint64_t size): KmerManipulator(size), m_format("ACGT") {}
KmerManipulatorACGT::~KmerManipulatorACGT() noexcept {}
uint64_t KmerManipulatorACGT::encode(const std::string word) {
    uint64_t hash = 0;
    char c;
    for(uint i = 0 ; i < m_size ; i++){   //We go through the sequence to encode each nucleotides
        hash <<= 2; // We shift 2 bits to the right
        c = word[i];    //Take each nucleotides to encode them
        uint64_t charval = 0; // 'A' = 0
        if(c == 'C'){
            charval = 1;    // 'C' = 1
        }else if(c == 'G'){
            charval = 2;    // 'G' = 2
        }else if(c == 'T'){
            charval = 3;    // 'T' = 3
        }
        hash += charval;    //creation of the hash for the given sequence
    }
    return hash;    //return the final hash of the sequence
}
std::string KmerManipulatorACGT::decode(uint64_t kmer) {
    string res(m_size, ' ');
    uint64_t lastIndex = res.size()-1;
    for(int i(0); i < m_size; i++){
        switch(kmer & 0x3){ //compares the decimal value of the first two bits
            case 0: res[lastIndex-i] = 'A'; break;
            case 1: res[lastIndex-i] = 'C'; break;
            case 2: res[lastIndex-i] = 'G'; break;
            case 3: res[lastIndex-i] = 'T'; break;
        }
        kmer >>= 2;
    }
    return res;
}
uint64_t KmerManipulatorACGT::getCanonical(const uint64_t kmer) {
    uint64_t reverseCompl = reverseComplement(kmer);
    return((kmer < reverseCompl)?kmer:reverseCompl);
}

__m256i KmerManipulatorACGT::getCanonicalAVX(const __m256i kmer) {
    __m256i tab = kmer;
    __m256i reverseCompl = reverseComplementAVX(kmer);  //call reverseComplement AVX version, 4 by 4 calls
    if(reverseCompl[0] < kmer[0]) tab = _mm256_insert_epi64(tab, reverseCompl[0], 0);   //if the reverse is the canonic one, we replace the ancient one by it
    if(reverseCompl[1] < kmer[1]) tab = _mm256_insert_epi64(tab, reverseCompl[1], 1);
    if(reverseCompl[2] < kmer[2]) tab = _mm256_insert_epi64(tab, reverseCompl[2], 2);
    if(reverseCompl[3] < kmer[3]) tab = _mm256_insert_epi64(tab, reverseCompl[3], 3);
    return tab;
}

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
 * AVX version of reverseComplement
 */
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
 * Construct the ConwayBromage object based on the sd_vector in parameter and the KmerManipulator.
 * @param sdv - An sd_vector which represents the sequence.
 * @param km - A KmerManipulator.
 */
ConwayBromage::ConwayBromage(sdsl::sd_vector<> const& sdv, KmerManipulator* km){
    m_kmerSize = (int)(log(sdv.size())/log(4));
    m_sequence = sdv; //copy of the sd_vector
    m_kmerManipulator = km;
}

/**
 * Transform sequences which are contain in a file in a sd_vector
 * @param kmerFlux - An istream of kmer. For example a file represented by an ifstream.
 * @param km - Object which contains all the information about encoding and decoding.
 */
ConwayBromage::ConwayBromage(istream& kmerFlux, KmerManipulator* km){
    m_kmerManipulator = km;
    string word("");
    string line("");
    kmerFlux >> word;   //Take the first word to analyze size of one k-mer
    kmerFlux.seekg(0, ios::beg);    //Return to the beginning of the file
    int KmerSize(word.size()); //Size of k_mer, it is the 'k'
    int numberOfKmer = 0;
    while(getline(kmerFlux, line)) //Counts the number of k-mer in the file
        numberOfKmer++;
 
    kmerFlux.clear();
    kmerFlux.seekg(0, ios::beg);    //Return to the beginning of the file
    //cout << "K-MER SIZE      : " << KmerSize << endl;
    //cout << "NUMBER OF K-MER : " << numberOfKmer << endl;
    uint64_t sdvSize(pow(4,KmerSize));//1 << (2*KmerSize));//pow(4,KmerSize);   //Creation of the total length to create the sd_vector_builder
    //cout << "SD VECTOR SIZE  : " << sdvSize << endl;
    
    sd_vector_builder builder(sdvSize, numberOfKmer);
    while(kmerFlux >> word){
        uint64_t kmer = m_kmerManipulator->encode(word);
        if(m_kmerManipulator->getCanonical(kmer) != kmer){
            cout << "The file is not completely canonical" << endl;
            exit(1); //EXIT_FAILURE
        }
        builder.set(kmer);
        kmerFlux >> word;
    }
    m_sequence = builder;
    m_kmerSize = KmerSize;
}
/**
 * Return the value stored at the block of index i.
 * @param i - the index.
 * @return true if the block of index i is set to one.
 */
bool ConwayBromage::operator[](uint64_t i) const{
    return m_sequence[i];
}
/**
 * AVX version of the [] operator
 * does not return elements, just said if one or more elements of the sequence are set to 1 or not
 */
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
 * Check if the given Kmer is present in the sequence.
 * @param Kmer - An integer representing the k-mer.
 * @return true if the k-mer is present.
 */
bool ConwayBromage::isPresent(uint64_t Kmer) const{
    int KmerSize = m_kmerSize-1;
    uint64_t limit = m_sequence.size() >> 2;
    if(Kmer >= limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << limit << endl;
        return false; 
    }  
    
    uint64_t next = Kmer << 2; //next contains XA where X is the Kmer
    if(m_sequence[m_kmerManipulator->getCanonical(next)])   return true; //getCanonical(XA) where X is the Kmer
    if(m_sequence[m_kmerManipulator->getCanonical(next+1)]) return true; //getCanonical(XC)
    if(m_sequence[m_kmerManipulator->getCanonical(next+2)]) return true; //getCanonical(XG)
    if(m_sequence[m_kmerManipulator->getCanonical(next+3)]) return true; //getCanonical(XT)
    
    int numberOfBitsToShift = KmerSize << 1;  
    if(m_sequence[m_kmerManipulator->getCanonical(Kmer)])   return true; //getCanonical(AX) where X is the Kmer
    if(m_sequence[m_kmerManipulator->getCanonical(Kmer + (1 << numberOfBitsToShift))]) return true; //getCanonical(CX)
    if(m_sequence[m_kmerManipulator->getCanonical(Kmer + (2 << numberOfBitsToShift))]) return true; //getCanonical(GX)
    if(m_sequence[m_kmerManipulator->getCanonical(Kmer + (3 << numberOfBitsToShift))]) return true; //getCanonical(TX)
    
    return false;
}
/**
 * AVX version of isPresent, for perf tests
 */
bool ConwayBromage::isPresentAVX(uint64_t Kmer) const {
    int KmerSize = m_kmerSize-1;
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
 * Returns the successors of a canonical kmer.
 * @param Kmer : the Kmer that we want to find its successors. Can either be canonical or not.
 * @return a uint8_t representing the (at most 8) successors of the kmer.
 */
uint8_t ConwayBromage::successors(uint64_t Kmer) const{
    uint8_t res = 0;
    int KmerSize = m_kmerSize-1;
    uint64_t limit = m_sequence.size() >> 2;
    if(Kmer >= limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << limit << endl;
        return res; //empty
    }
    //we build the eight possible successors    
    //4 next pmers
    uint64_t Pmer = Kmer << 2; //<-> XA where X is the Kmer
    for(int i = 0; i < 4; i++){
        if(m_sequence[m_kmerManipulator->getCanonical(Pmer)]){ //if the pmer is present 
            //if compressedSeq.size() = 256 (P=4) then compressedSeq.size() >> 2 = 64 ==> AACA%64 ==> ACA (i_next_kmer)
            string PmerStr = m_kmerManipulator->decode(Pmer);
            char lastLetter = PmerStr[PmerStr.size()-1];
            //we set the concerned bit to 1
            if(lastLetter == 'A') 
                res = res ^ 128; 
            else if(lastLetter == 'C')      
                res = res ^ 64;
            else if(lastLetter == 'G') 
                res = res ^ 32;
            else if(lastLetter == 'T')
                res = res ^ 16;
        }
        Pmer++; //equals to XA then XC then XG then XT
    }
    //4 predecessors
    int numberOfBitsToShift = KmerSize << 1;
    for(int i = 0; i < 4; i++){
        Pmer = Kmer + (i << numberOfBitsToShift); //equals to AX then CX then GX then TX where X is the Kmer
        if(m_sequence[m_kmerManipulator->getCanonical(Pmer)]){
            string PmerStr = m_kmerManipulator->decode(Pmer);
            char lastLetter = PmerStr[0];
            //we set the concerned bit to 1
            if(lastLetter == 'A')
                res = res ^ 8;
            else if(lastLetter == 'C')      
                res = res ^ 4;
            else if(lastLetter == 'G')
                res = res ^ 2;
            else if(lastLetter == 'T')
                res = res ^ 1;
        }
    }
    return res;
}


/**
 * Returns the number of 1-bit to the left of position index
 * @param index the position
 * @param v an sd_vector
 * @return a count
 */
uint64_t ConwayBromage::rank1bit(uint64_t index){
    sd_vector<>::rank_1_type sdb_rank(&m_sequence);
    return sdb_rank(index);
}

/**
 * Returns the number of 0-bit to the left of position index
 * @param index the position 
 * @param v an sd_vector
 * @return a count
 */
uint64_t ConwayBromage::rank0bit(uint64_t index){
    sd_vector<>::rank_0_type sdb_rank(&m_sequence);
    return sdb_rank(index);
}

/**
 * Returns the position of the index-th 1-bit.
 * @param index the position 
 * @param v an sd_vector
 * @return an index
 */
uint64_t ConwayBromage::select1bit(uint64_t index){
    sd_vector<>::select_1_type sdb_sel(&m_sequence);
    return sdb_sel(index);
}

/**
 * Returns the position of the index-th 0-bit.
 * @param index the position 
 * @param v an sd_vector
 * @return an index
 */
uint64_t ConwayBromage::select0bit(uint64_t index){
    sd_vector<>::select_0_type sdb_sel(&m_sequence);
    return sdb_sel(index);
}

/**
 * Return the size of the sequence.
 * @return an int
 */
int ConwayBromage::getKmerSize(){
    return m_kmerSize;
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
 */
KmerManipulator* ConwayBromage::getKmerManipulator(){
    return m_kmerManipulator;
}

