//
// Created by Alexandra on 18/06/2020.
//

#ifndef CONWAYBROMAGELIB_FROMFILETOSDVECTORVARIATIONS_H
//fromFileToSdVector variations
sdsl::sd_vector<>fromFileToSdVectorWithReverse(std::string path); //original with reverse complement
sdsl::sd_vector<>fromFileToSdVectorEcoli(std::string path); //fastest ASCII version, use encoding ACTG
sdsl::sd_vector<>fromFileToSdVectorWithReverseEcoli(std::string path);  //fastest ASCII version with reverse complement merging
sdsl::sd_vector<> fromFileToSdVector_TXTversion(std::string path);  //version with txt file
#define CONWAYBROMAGELIB_FROMFILETOSDVECTORVARIATIONS_H

#endif //CONWAYBROMAGELIB_FROMFILETOSDVECTORVARIATIONS_H
