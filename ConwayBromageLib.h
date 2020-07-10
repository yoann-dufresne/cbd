#include <sdsl/sd_vector.hpp>
#include <immintrin.h>  //for AVX/AVX2 use
#ifndef CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
#define CONWAYBROMAGELIB_CONWAYBROMAGELIB_H

const int ALPHABET(4);  // Size of the alphabet : E = {A, C, G, T}
const std::string NUCLEOTIDES [4] = {"A", "C", "G", "T"};

std::vector<uint64_t> successorTranslator(int successors, uint64_t compressedKMer, uint64_t size, bool format);
//Functions for tests only
bool isCanonical (uint64_t kmer, uint64_t kmerSize, bool encodingIsACGT);
std::vector<uint64_t> successorCounter(uint64_t compressedKMer, sdsl::sd_vector<>currentCompressedSeq, bool format);

//POO for KmerManipulator
class KmerManipulator{  //abstract class
public:
    KmerManipulator(uint64_t size); //Constructor
    //abstract function
    /* function to encode k-mer sequences, depends on encoding : see daughter classes
    *@param word - the string we want to encode into a uint64_t
    *@return a uint64_t which represent the encoding version of the word
    */
    virtual uint64_t encode(const std::string &word) = 0;
    /*
    * Returns a string which is the word which corresponds to the uint64_t value, depends on encoding : see daughter classes
    * @param kmer - the value.
    * @return a string representing the decoding version of the uint64_t.
    */
    virtual std::string decode(uint64_t kmer) = 0;
    /*
    *calculate the canonical version of a compressed k-mer, depends on encoding : see daughter classes
    *@param kmer - a uint64_t which represents the compressed kmer we want to study
    *@return a uint64_t which is the compressed kmer itself if it is already canonical, or its canonical version
    */
    virtual uint64_t getCanonical(const uint64_t kmer) = 0;
    /**
     * AVX version of getCanonical
     */
    virtual __m256i getCanonicalAVX(const __m256i kmer) = 0;
    /* Calculate the reverse complement, depends on encoding : see daughter classes
    * @param kmer - a uin64_t which represent the compressed version of a k-mer
    * @return a uin64_t which represent the compressed version of the reverse complement of the given k-mer
    */
    virtual uint64_t reverseComplement(const  uint64_t kmer) = 0;
    /**
     * AVX version of reverseComplement
     */
    virtual __m256i reverseComplementAVX(const __m256i kmer) = 0;
    
    /**
    * Returns the reverse complement of a nucleotide.
    */
    virtual uint8_t reverseComplementOfNucleotide(const uint8_t nucleotide) = 0;
    
    /**
    * Returns the corresponding caracter to the nucleotide's value.
    */
    virtual char decodeNucleotide(const uint8_t nucleotide) = 0;
    
    /**
    * Returns the size attribute.
    */
    virtual int getSize() = 0;
    
    virtual ~KmerManipulator();     //Desctructor
protected:
    uint64_t m_size;
};
//KmerManipulator for ACGT encoding
class KmerManipulatorACGT : public KmerManipulator{
public:
    KmerManipulatorACGT(uint64_t size);
    uint64_t encode(const std::string &word);
    std::string decode(uint64_t kmer);
    uint64_t getCanonical(const uint64_t kmer);
    __m256i getCanonicalAVX(const __m256i kmer);
    uint64_t reverseComplement(const uint64_t kmer);
    __m256i reverseComplementAVX(const __m256i kmer);
    uint8_t reverseComplementOfNucleotide(const uint8_t nucleotide);
    char decodeNucleotide(const uint8_t nucleotide);
    ~KmerManipulatorACGT();
    int getSize();
private:
    std::string m_format;
};
//KmerManipulator for ACTG encoding
class KmerManipulatorACTG : public KmerManipulator{
public:
    KmerManipulatorACTG(uint64_t size);
    uint64_t encode(const std::string &word);
    std::string decode(uint64_t kmer);
    uint64_t getCanonical(const uint64_t kmer);
    __m256i getCanonicalAVX(const __m256i kmer);
    uint64_t reverseComplement(const u_int64_t kmer);
    __m256i reverseComplementAVX(const __m256i kmer);
    uint8_t reverseComplementOfNucleotide(const uint8_t nucleotide);
    char decodeNucleotide(const uint8_t nucleotide);
    ~KmerManipulatorACTG();
    int getSize();
private:
    std::string m_format;
};

class ConwayBromage{
private:
    int m_kmerSize;                     //PmerSize actually
    sdsl::sd_vector<> m_sequence;       //the compressed k-mer sequence
    KmerManipulator* m_kmerManipulator; //stores information about encode/decode
    uint64_t m_limit;
    //cache used in isPresent and successors in order to go faster
    int m_correspondingBitValueForNextPmers[4];
    int m_correspondingBitValueForPrevPmers[4];
    int m_numberOfBitsToShift;
    uint64_t m_RC[4];                      //{reverse complement of nucleotide (0), ... , reverse complement of of nucleotide(3)}
    uint64_t m_RC_shifted[4];              //{reverse complement of of nucleotide(0) << (2*(m_kmerSize-1)), ... , reverse complement of (3) << (2*(m_kmerSize-1))}
    uint64_t m_nucleotides_shifted[4];     //{0 << m_numberOfBitsToShift, ... , 3 << m_numberOfBitsToShift}
    
public:
    //constructors
    ConwayBromage(std::istream& kmerFlux, KmerManipulator* km);
    ConwayBromage(sdsl::sd_vector<> const& sdv, KmerManipulator* km);
    //principal functions
    bool isPresent (uint64_t Kmer) const;
    //AVX prototype of isPresent
    bool isPresentAVX(uint64_t Kmer) const;
    uint8_t successors(uint64_t Kmer) const;
    //operators on the sequence
    bool operator[](uint64_t i) const;
    //AVX version of the operator
    bool operator[](__m256i i) const;
    uint64_t size() const;
    //getters
    int getKmerSize();
    sdsl::sd_vector<> getSequence();
    KmerManipulator* getKmerManipulator();
    /*
     * Functions for isPresent performance tests
     * Use them as friend of ConwayBromage to allow them to use attribute easily
     */
    friend sdsl::int_vector<> ratioForIsPresent(int ratioIn, int nbOfOnes, ConwayBromage cb);
    friend void metricForIsPresent();
};

#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
