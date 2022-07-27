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
class Intermediate{
    private :
        uint64_t mask=0b1111111111111111000000000000000000000000000000000000000000000000; // mask to separe the 16 first bits from the 48 other
        uint64_t mask2=0b0000000000000000111111111111111111111111111111111111111111111111;
    public :
    unsigned int nbv;
    std::vector<bm::bvector<>> vbv;
    Intermediate();
    Intermediate(int nb);
    Intermediate(std::vector<bm::bvector<>>& tmp2,int a);
    void set(uint64_t id);
    int present(uint64_t id)const;
    void optimize();
};

/// ConwayBromage : permits to store k-mers and apply two operations on them : 'contains' and 'successors'
class ConwayBromageBM : public ConwayBromage{
private:
    Intermediate m_sequence;          //the succint bit vector which stores the k-mers
    static int serializeaux(std::ostream& output, bm::bvector<>& sequence);
    static bm::bvector<> deserializeaux(std::istream& bitVector);
public:
    //constructors
    ConwayBromageBM(std::istream& kmerFlux, KmerManipulator* km);
    ConwayBromageBM(Intermediate&  sdv, KmerManipulator* km);
    //principal functions
    bool contains (uint64_t Kmer) const;
    uint8_t successors(uint64_t Kmer) const;
    //getters
    Intermediate getSequence();
    //serializer
    int serialize(std::string dirpath);
        static int read_bvector(std::ifstream& bv_file, bm::bvector<>& bv, bm::serializer<bm::bvector<> >::buffer& sbuf);

    static ConwayBromageBM deserialize(std::string dirpath,KmerManipulator* km);
    /*
     * Functions for isPresent performance tests
     * Use them as friend of ConwayBromage to allow them to use attribute easily
     * Not needed for the global functioning
     */

    friend bm::bvector<> ratioForIsPresent(int ratioIn, int nbOfOnes, ConwayBromageBM cb);
    friend void metricForIsPresent();

};

#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
