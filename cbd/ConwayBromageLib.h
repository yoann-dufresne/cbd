#include <sdsl/sd_vector.hpp>
#include <immintrin.h>  //for AVX/AVX2 use
#include "Kmanip.h"
#include "bm64.h"
#include "bmserial.h"
#include "bmundef.h" /* clear the pre-proc defines from BM */
#include <filesystem>
#ifndef CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
#define CONWAYBROMAGELIB_CONWAYBROMAGELIB_H


class ConwayBromage{  //abstract class
protected:
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
    ConwayBromage(KmerManipulator* km);
    //principal functions
    virtual bool contains (uint64_t Kmer) const=0;
    virtual uint8_t successors(uint64_t Kmer) const=0;

    //getters
    virtual int getKmerSize();
    virtual KmerManipulator* getKmerManipulator();
};

/// ConwayBromage : permits to store k-mers and apply two operations on them : 'contains' and 'successors'
class ConwayBromageSD : public ConwayBromage{
private:
    sdsl::sd_vector<> m_sequence;          //the sparse bit vector which stores the k-mers        
public:
    //constructors
    ConwayBromageSD(std::istream& kmerFlux, KmerManipulator* km);
    ConwayBromageSD(sdsl::sd_vector<> const& sdv, KmerManipulator* km);
    ConwayBromageSD(KmerManipulator* km);
    //principal functions
    bool contains (uint64_t Kmer) const;
    uint8_t successors(uint64_t Kmer) const;

    //getters
    sdsl::sd_vector<> getSequence();
    int serialize(std::string filename);
    static ConwayBromageSD deserialize(std::string filename,KmerManipulator* km);

    /*
     * Functions for isPresent performance tests
     * Use them as friend of ConwayBromage to allow them to use attribute easily
     * Not needed for the global functioning
     */
    friend sdsl::int_vector<> ratioForIsPresent(int ratioIn, int nbOfOnes, ConwayBromageSD cb);
    friend void metricForIsPresent();
};
#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
