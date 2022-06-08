#include <Kmanip.h>

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
uint64_t KmerManipulatorACTG::encode(const std::string &word) {
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
std::string KmerManipulatorACTG::decode(uint64_t kmer) {
    static const char valueToCaracter[4] = {'A','C','T','G'}; //declared and initialized only during the first call to this method
    std::string res(m_size, ' ');
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
uint64_t KmerManipulatorACGT::encode(const std::string &word) {
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
std::string KmerManipulatorACGT::decode(uint64_t kmer) {
    static const char valueToCaracter[4] = {'A','C','G','T'}; //declared and initialized only during the first call to this method
    std::string res(m_size, ' ');
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
