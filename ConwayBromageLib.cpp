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
   uint64_t charval(0);
   for(uint i = 0 ; i < m_size ; i++){   //We go through the sequence to encode each nucleotides
       hash <<= 2; // We shift 2 bits to the right
       switch (word[i]) {
           case 'C':
                charval = 1;
                break;
           case 'T':
                charval = 2;
                break;
           case 'G':
                charval = 3;
                break;
           default:
                charval = 0;
       }
       hash += charval;
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
    __m256i reverseCompl = reverseComplementAVX(kmer);  //creation of a __m256i which contains reverse complement of the __m256i k-mers
    __m256i eq = _mm256_cmpeq_epi64(kmer,reverseCompl); //Case where forward and reverse are the same
    __m256i cmp1 = _mm256_or_si256(_mm256_cmpgt_epi64(kmer, reverseCompl), eq);  //kmer < reverseCompl => 0, else -1
    __m256i cmp2 = _mm256_cmpgt_epi64(reverseCompl, kmer);  //reverseCompl < kmer => 0, else -1
    __m256i and1 = _mm256_and_si256(cmp2, kmer);    //delete elements wich are higher than their reverse
    __m256i and2 = _mm256_and_si256(cmp1, reverseCompl);    //delete reverses which are higher than their forward
    __m256i ult = _mm256_or_si256(and1, and2);  //final construction
    return ult;
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

__m256i KmerManipulatorACTG::reverseComplementAVX(const __m256i kmer){
    //We initialize variables like __mm256i vectors
    __m256i AAA = _mm256_set1_epi64x(0xAAAAAAAAAAAAAAAA);
    __m256i val1 = _mm256_set1_epi64x(0x3333333333333333);
    __m256i val2 = _mm256_set1_epi64x(0x0F0F0F0F0F0F0F0F);
    __m256i val3 = _mm256_set1_epi64x(0x00FF00FF00FF00FF);
    __m256i val4 = _mm256_set1_epi64x(0x0000FFFF0000FFFF);
    __m256i val5 = _mm256_set1_epi64x(0x00000000FFFFFFFF);
    __m256i res = kmer;  // <=> res = kmer
    //translation in AVX instruction of the original reverseComplement : CARE ABOUT UNSIGNED INT !!!
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
 */
uint8_t KmerManipulatorACTG::reverseComplementOfNucleotide(const uint8_t nucleotide){
    if(nucleotide == 0) //A
        return 2; //T
    if(nucleotide == 1) //C
        return 3; //G
    if(nucleotide == 2) //T
        return 0; //A
    if(nucleotide == 3) //G
        return 1; //C
}

/**
 * Returns the corresponding caracter to the nucleotide's value.
 */
char KmerManipulatorACTG::decodeNucleotide(const uint8_t nucleotide){
    if(nucleotide == 0)
        return 'A';
    if(nucleotide == 1)
        return 'C';
    if(nucleotide == 2)
        return 'T';
    if(nucleotide == 3)
        return 'G';
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
    __m256i reverseCompl = reverseComplementAVX(kmer);  //creation of a __m256i which contains reverse complement of the __m256i k-mers
    __m256i eq = _mm256_cmpeq_epi64(kmer,reverseCompl); //Case where reverse and forward are the same
    __m256i cmp1 = _mm256_or_si256(_mm256_cmpgt_epi64(kmer, reverseCompl), eq);  //kmer < reverseCompl => 0, else -1
    __m256i cmp2 = _mm256_cmpgt_epi64(reverseCompl, kmer);  //reverseCompl < kmer => 0, else -1
    __m256i and1 = _mm256_and_si256(cmp2, kmer);    //delete elements wich are higher than their reverse
    __m256i and2 = _mm256_and_si256(cmp1, reverseCompl);    //delete reverses which are higher than their forward
    __m256i ult = _mm256_or_si256(and1, and2);  //final construction
    return ult;
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
 * Returns the reverse complement of a nucleotide.
 */
uint8_t KmerManipulatorACGT::reverseComplementOfNucleotide(const uint8_t nucleotide){
    if(nucleotide == 0) //A
        return 3; //T
    if(nucleotide == 1) //C
        return 2; //G
    if(nucleotide == 2) //G
        return 1; //C
    if(nucleotide == 3) //T
        return 0; //A
}

/**
 * Returns the corresponding caracter to the nucleotide's value.
 */
char KmerManipulatorACGT::decodeNucleotide(const uint8_t nucleotide){
    if(nucleotide == 0)
        return 'A';
    if(nucleotide == 1)
        return 'C';
    if(nucleotide == 2)
        return 'G';
    if(nucleotide == 3)
        return 'T';
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
    m_limit = (m_sequence.size() >> 2) - 1;
    
    //initialization of the cache
    int nbOfBitsToShift = 2 * (m_kmerSize-1);
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
 * Transform sequences which are contain in a file in a sd_vector
 * @param kmerFlux - An istream of kmer. For example a file represented by an ifstream.
 * @param km - Object which contains all the information about encoding and decoding.
 */
ConwayBromage::ConwayBromage(istream& kmerFlux, KmerManipulator* km){
    m_kmerManipulator = km;
    m_kmerSize = m_kmerManipulator->getSize(); //Size of k_mer, it is the 'k'

    string line("");
    int numberOfKmer = 0;
    while(getline(kmerFlux, line)) //Counts the number of k-mer in the file
        numberOfKmer++;
 
    kmerFlux.clear();
    kmerFlux.seekg(0, ios::beg);    //Return to the beginning of the file
    uint64_t one = 1;
    uint64_t sdvSize = one << (2*m_kmerSize); //Creation of the total length to create the sd_vector_builder
    //cout << "K-MER SIZE      : " << KmerSize << endl;
    //cout << "NUMBER OF K-MER : " << numberOfKmer << endl;
    //cout << "SD VECTOR SIZE  : " << sdvSize;
    
    sd_vector_builder builder(sdvSize, numberOfKmer);
    while(getline(kmerFlux, line)){
        uint64_t k = m_kmerManipulator->encode(line);
        if(k >  m_kmerManipulator->reverseComplement(k)){
            cout << "The file is not completely canonical" << endl;
            exit(1);
        }
        builder.set(k);   
    }
    
    m_sequence = builder;
    m_kmerSize = KmerSize;
    m_limit = (m_sequence.size() >> 2) - 1;
    
    //initialization of the cache
    int nbOfBitsToShift = 2 * (m_kmerSize-1);
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
    //Try for an AVX version of the new isPresent
    /*if(Kmer > m_limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << m_limit << endl;
        return false;
    }
    uint64_t RC_Kmer = m_kmerManipulator->reverseComplement(Kmer) >> 2;
    uint64_t RC_Kmer_ShiftedOf2Bits = RC_Kmer << 2;
    uint64_t next = Kmer << 2;

    __m256i PmerNext = _mm256_setr_epi64x(next, next+1, next+2, next+3);    //__m256i for all the forward nexts
    __m256i RC_PmerNext = _mm256_setr_epi64x(RC_Kmer + m_RC_shifted[0], RC_Kmer + m_RC_shifted[1],  //__m256i for all the reverseComplement nexts
                                            RC_Kmer + m_RC_shifted[2], RC_Kmer + m_RC_shifted[3]);
    __m256i cmp1 = _mm256_cmpgt_epi64(PmerNext, RC_PmerNext);   //Apply the same thing than in getCanonicalAVX
    __m256i cmp2 = _mm256_cmpgt_epi64(RC_PmerNext, PmerNext);
    __m256i and1 = _mm256_and_si256(cmp2, PmerNext);
    __m256i and2 = _mm256_and_si256(cmp1, RC_PmerNext);
    __m256i ultN = _mm256_or_si256(and1, and2);
    if(operator[](ultN)) return true;   //verify existence
    __m256i PmerPrev = _mm256_setr_epi64x(Kmer + m_nucleotides_shifted[0], Kmer + m_nucleotides_shifted[1], //__m256i for all the forward previous
                                         Kmer + m_nucleotides_shifted[2], Kmer + m_nucleotides_shifted[3]);
    __m256i RC_PmerPrev = _mm256_setr_epi64x(RC_Kmer_ShiftedOf2Bits + m_RC[0], RC_Kmer_ShiftedOf2Bits + m_RC[1],    //__m256i for all the reverseComplement previous
                                            RC_Kmer_ShiftedOf2Bits + m_RC[2], RC_Kmer_ShiftedOf2Bits + m_RC[3]);
    cmp1 = _mm256_cmpgt_epi64(PmerPrev, RC_PmerPrev);
    cmp2 = _mm256_cmpgt_epi64(RC_PmerPrev, PmerPrev);
    and1 = _mm256_and_si256(cmp2, PmerPrev);
    and2 = _mm256_and_si256(cmp1, RC_PmerPrev);
    __m256i ultP = _mm256_or_si256(and1, and2);
    if(operator[](ultP)) return true;
    return false;*/
}

/**
 * Returns the successors of a canonical kmer.
 * @param Kmer : the Kmer that we want to find its successors. Can either be canonical or not.
 * @return a uint8_t representing the (at most 8) successors of the kmer.
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


