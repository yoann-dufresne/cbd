#include <sdsl/sd_vector.hpp>
#ifndef CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
#define CONWAYBROMAGELIB_CONWAYBROMAGELIB_H

const int ALPHABET(4);  // Size of the alphabet : E = {A, C, G, T}
const std::string NUCLEOTIDES [4] = {"A", "C", "G", "T"};

uint64_t rank1bit(uint64_t index, const sdsl::sd_vector<> &v);
uint64_t rank0bit(uint64_t index, const sdsl::sd_vector<> &v);
uint64_t select1bit(uint64_t index, const sdsl::sd_vector<> &v);
uint64_t select0bit(uint64_t index, const sdsl::sd_vector<> &v);

uint64_t encode(std::string word, uint64_t size);
std::string decode(uint64_t seq, uint64_t size);
sdsl::sd_vector<>fromFileToSdVector(std::string path);
//Prototype of others functions
bool isThisKMerHere(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);
std::vector<uint64_t> next(uint64_t nonCompressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);
std::vector<uint64_t> previous(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);  

uint64_t reverseComplement(std::string seq, uint64_t len);
std::vector<sdsl::sd_vector<>>fromFileToSdVectorWithReverseEcoli(std::string path);
#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
