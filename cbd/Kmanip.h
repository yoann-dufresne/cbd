#ifndef KMANIP_H
#define KMANIP_H
#include <cstdlib>
#include <string>
using namespace std;
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

///KmerManipulator for ACGT encoding : permits to apply operations on k-mers in ACGT format
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

///KmerManipulator for ACTG encoding : permits to apply operations on k-mers in ACTG format
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
#endif //KMANIP_H