#include <sdsl/sd_vector.hpp>
#include <immintrin.h>  //for AVX/AVX2 use
#ifndef CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
#define CONWAYBROMAGELIB_CONWAYBROMAGELIB_H

std::vector<uint64_t> successorTranslator(int successors, uint64_t compressedKMer, uint64_t size, bool format);

class KmerManipulator{  //abstract class
protected:
    uint64_t m_size;    //size of the k-mers
public:
    KmerManipulator(uint64_t size);
    virtual ~KmerManipulator();  
    //operations on k-mers
    virtual uint64_t encode(const std::string &word) = 0;
    virtual std::string decode(uint64_t kmer) = 0;
    virtual uint64_t getCanonical(const uint64_t kmer) = 0;
    virtual uint64_t reverseComplement(const  uint64_t kmer) = 0;
    virtual uint8_t reverseComplementOfNucleotide(const uint8_t nucleotide) = 0;
    virtual char decodeNucleotide(const uint8_t nucleotide) = 0;
    virtual int getSize() = 0;  
};

//KmerManipulator for ACGT encoding : permits to apply operations on k-mers in ACGT format
class KmerManipulatorACGT : public KmerManipulator{
private:
    std::string m_format;   //always equal to "ACGT"
public:
    KmerManipulatorACGT(uint64_t size);
    ~KmerManipulatorACGT();
    //operations on k-mers
    uint64_t encode(const std::string &word);
    std::string decode(uint64_t kmer);
    uint64_t getCanonical(const uint64_t kmer);
    uint64_t reverseComplement(const uint64_t kmer);
    uint8_t reverseComplementOfNucleotide(const uint8_t nucleotide);
    char decodeNucleotide(const uint8_t nucleotide);
    int getSize();
};

//KmerManipulator for ACTG encoding : permits to apply operations on k-mers in ACTG format
class KmerManipulatorACTG : public KmerManipulator{
private:
    std::string m_format; // always equal to "ACTG"
public:
    KmerManipulatorACTG(uint64_t size);
    ~KmerManipulatorACTG();
    //operations on k-mers
    uint64_t encode(const std::string &word);
    std::string decode(uint64_t kmer);
    uint64_t getCanonical(const uint64_t kmer);
    uint64_t reverseComplement(const u_int64_t kmer);
    uint8_t reverseComplementOfNucleotide(const uint8_t nucleotide);
    char decodeNucleotide(const uint8_t nucleotide);
    int getSize();
};

//ConwayBromage : permits to store k-mers and apply two operations on them : 'contains' and 'successors'
class ConwayBromage{
private:
    sdsl::sd_vector<> m_sequence;          //the sparse bit vector which stores the k-mers
    KmerManipulator* m_kmerManipulator;    //carry information about the k-mers
    uint64_t m_limit;
    
    //cache used in contains and successors functions in order to go faster
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
    bool contains (uint64_t Kmer) const;
    uint8_t successors(uint64_t Kmer) const;
    uint64_t size() const;
    //getters
    int getKmerSize();
    sdsl::sd_vector<> getSequence();
    KmerManipulator* getKmerManipulator();
    /*
     * Functions for isPresent performance tests
     * Use them as friend of ConwayBromage to allow them to use attribute easily
     * Not needed for the global functioning
     */
    friend sdsl::int_vector<> ratioForIsPresent(int ratioIn, int nbOfOnes, ConwayBromage cb);
    friend void metricForIsPresent();
};

#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
