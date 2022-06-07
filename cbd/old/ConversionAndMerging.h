//
// Created by alexa on 18/06/2020.
//

#ifndef CONWAYBROMAGELIB_CONVERSIONANDMERGING_H
#define CONWAYBROMAGELIB_CONVERSIONANDMERGING_H
//merge
sdsl::sd_vector<> merge (const sdsl::sd_vector<> &a, const sdsl::sd_vector<> &b, int nb_of_1_in_A, int nb_of_1_in_B);

//Base64<->Base10
std::string convertToBase64(uint64_t valueInBase10);
uint64_t convertFromBase64ToBase10(std::string valueInBase64);

//reverseComplement variations
uint64_t reverseComplement(std::string seq, uint64_t len);  //original
#endif //CONWAYBROMAGELIB_CONVERSIONANDMERGING_H
