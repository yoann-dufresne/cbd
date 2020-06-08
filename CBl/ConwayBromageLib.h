//
// Created by Alexandra and Murat on 05/06/2020.
//
#include <sdsl/sd_vector.hpp>

#ifndef CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
#define CONWAYBROMAGELIB_CONWAYBROMAGELIB_H

const int ALPHABET(4);  // Size of the alphabet : E = {A, C, G, T}
const std::string NUCLEOTIDES [4] = {"A", "C", "G", "T"};

uint64_t encode(std::string word, uint64_t size);
std::string decode(uint64_t seq, uint64_t size);    //Proto //Murat
sdsl::sd_vector<>fromFileToSdVector(std::string path);
//Prototype of others functions
bool isThisKMerHere(std::string nonCompressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);
std::vector<std::string> nexts(std::string nonCompressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);   //Murat
std::vector<std::string> previous(std::string nonCompressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);

#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
