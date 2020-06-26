#include <sdsl/sd_vector.hpp>
#ifndef CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
#define CONWAYBROMAGELIB_CONWAYBROMAGELIB_H

const int ALPHABET(4);  // Size of the alphabet : E = {A, C, G, T}
const std::string NUCLEOTIDES [4] = {"A", "C", "G", "T"};

//rank and select
uint64_t rank1bit(uint64_t index, const sdsl::sd_vector<> &v);
uint64_t rank0bit(uint64_t index, const sdsl::sd_vector<> &v);
uint64_t select1bit(uint64_t index, const sdsl::sd_vector<> &v);
uint64_t select0bit(uint64_t index, const sdsl::sd_vector<> &v);
//Encode/Decode lexicographical version
uint64_t encode(std::string word, uint64_t size);
std::string decode(uint64_t seq, uint64_t size);
//Encode/decode fastest ASCII version
std::string decodeEcoli(uint64_t seq, uint64_t size);
uint64_t encodeEcoli(std::string word, uint64_t size);
//fromFileToSdVector variations
sdsl::sd_vector<>fromFileToSdVector(std::string path); //can now choose the encoding format, canonical k-mer
sdsl::sd_vector<>fromFileToSdVectorChooser(std::string path, std::string format);
//isThisKMerHere
bool isThisKMerHere(uint64_t Kmer, sdsl::sd_vector<> const& compressedSeq, bool encodingIsAGCT);
//Successors of a canonical Kmer
uint64_t getCanonical (uint64_t kmer, uint64_t kmerSize, bool encodingIsACGT);
std::vector<uint64_t> successors(uint64_t Kmer, sdsl::sd_vector<> const& compressedSeq, bool encodingIsACGT);
std::vector<uint64_t> successorTranslator(int successors, uint64_t compressedKMer, uint64_t size, bool format);
//reverseComplement variations
uint64_t reverseComplementLexico (const uint64_t mer, uint64_t kmerSize); //faster for lexicographical order
u_int64_t reverseComplementGATBLibEcoli (const u_int64_t x, uint64_t sizeKmer); //the fastest : fastest ASCII version

//Functions for tests only
bool isCanonical (uint64_t kmer, uint64_t kmerSize, bool encodingIsACGT);
std::vector<uint64_t> successorCounter(uint64_t compressedKMer, sdsl::sd_vector<>currentCompressedSeq, bool format);

//For perf tests
uint64_t encodeEcoliSwitchVers(std::string word, uint64_t size);
std::string decodeEcoliIfVers(uint64_t seq, uint64_t size);

//POO for KmerManipulator
class KmerManipulator{  //abstract class
public:
    KmerManipulator(uint64_t size); //Constructor
    //abstract function
    /* function to encode k-mer sequences, depends on encoding : see daughter classes
    *@param word - the string we want to encode into a uint64_t
    *@return a uint64_t which represent the encoding version of the word
    */
    virtual uint64_t encode(const std::string word) = 0;
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
    /* Calculate the reverse complement, depends on encoding : see daughter classes
    * @param kmer - a uin64_t which represent the compressed version of a k-mer
    * @return a uin64_t which represent the compressed version of the reverse complement of the given k-mer
    */
    virtual uint64_t reverseComplement(const  uint64_t kmer) = 0;
    virtual ~KmerManipulator();     //Desctructor
protected:
    uint64_t m_size;
};
//KmerManipulator for ACGT encoding
class KmerManipulatorACGT : public KmerManipulator{
public:
    KmerManipulatorACGT(uint64_t size);
    uint64_t encode(const std::string word);
    std::string decode(uint64_t kmer);
    uint64_t getCanonical(const uint64_t kmer);
    uint64_t reverseComplement(const uint64_t kmer);
    ~KmerManipulatorACGT();
private:
    std::string m_format;
};
//KmerManipulator for ACTG encoding
class KmerManipulatorACTG : public KmerManipulator{
public:
    KmerManipulatorACTG(uint64_t size);
    uint64_t encode(const std::string word);
    std::string decode(uint64_t kmer);
    uint64_t getCanonical(const uint64_t kmer);
    uint64_t reverseComplement(const u_int64_t kmer);
    ~KmerManipulatorACTG();
private:
    std::string m_format;
};

class ConwayBromage{
private:
    int m_kmerSize;                     //PmerSize actually
    sdsl::sd_vector<> m_sequence;       //the compressed k-mer sequence
    KmerManipulator* m_kmerManipulator; //stores information about encode/decode
    
public:
    //constructors
    ConwayBromage(std::istream& kmerFlux, KmerManipulator* km);
    ConwayBromage(sdsl::sd_vector<> const& sdv, KmerManipulator* km);
    //principal functions
    bool isPresent (uint64_t Kmer) const;
    uint8_t successors(uint64_t Kmer) const;
    //operators on the sequence
    bool operator[](uint64_t i) const;
    uint64_t size() const;
    //rank and select operators
    uint64_t rank1bit  (uint64_t index);
    uint64_t rank0bit  (uint64_t index);
    uint64_t select1bit(uint64_t index);
    uint64_t select0bit(uint64_t index);
};

#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
