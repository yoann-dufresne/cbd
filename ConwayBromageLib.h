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
bool isThisKMerHere(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq, bool format);
//Successors of a canonical Kmer
uint64_t getCanonical (uint64_t kmer, uint64_t kmerSize, bool encodingIsACGT);
std::vector<uint64_t> successors(uint64_t Kmer, sdsl::sd_vector<> const& compressedSeq, bool encodingIsACGT);
//reverseComplement variations
uint64_t reverseComplementLexico (const uint64_t mer, uint64_t kmerSize); //faster for lexicographical order
u_int64_t reverseComplementGATBLibEcoli (const u_int64_t x, uint64_t sizeKmer); //the fastest : fastest ASCII version



#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
