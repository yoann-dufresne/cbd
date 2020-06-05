//
// Created by Alexandra and Murat on 05/06/2020.
//
#include <sdsl/sd_vector.hpp>
#ifndef CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
#define CONWAYBROMAGELIB_CONWAYBROMAGELIB_H

const int ALPHABET(4);  // Size of the alphabet : E = {A, C, G, T}

uint64_t encode(std::string word, uint64_t size);
sdsl::sd_vector<>fromFileToSdVector(std::string path);

#endif //CONWAYBROMAGELIB_CONWAYBROMAGELIB_H
