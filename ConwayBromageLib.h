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
sdsl::sd_vector<>fromFileToSdVector(std::string path);  //original, lexicographical order
sdsl::sd_vector<>fromFileToSdVectorWithReverse(std::string path); //original with reverse
sdsl::sd_vector<>fromFileToSdVectorEcoli(std::string path); //fastest ASCII version
sdsl::sd_vector<>fromFileToSdVectorWithReverseEcoli(std::string path);  //fastest ASCII version with reverse merging

//isThisKMerHere
bool isThisKMerHere(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);

//next and previous, insensitive to the encoding version ?
std::vector<uint64_t> next(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq); 
std::vector<uint64_t> previous(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);

//reverseComplement variations
uint64_t reverseComplement(std::string seq, uint64_t len);  //original
uint64_t reverseComplementLexico (const uint64_t mer, uint64_t kmerSize); //faster for lexicographical order
u_int64_t reverseComplementGATBLibEcoli (const u_int64_t x, uint64_t sizeKmer); //the fastest : fastest ASCII version

//merge
sdsl::sd_vector<> merge (const sdsl::sd_vector<> &a, const sdsl::sd_vector<> &b, int nb_of_1_in_A, int nb_of_1_in_B);


#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
