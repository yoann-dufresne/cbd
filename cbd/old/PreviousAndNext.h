//
// Created by alexa on 18/06/2020.
//

#ifndef CONWAYBROMAGELIB_PREVIOUSANDNEXT_H
#define CONWAYBROMAGELIB_PREVIOUSANDNEXT_H

//next and previous
std::vector<uint64_t> next(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);
std::vector<uint64_t> previous(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq);


#endif //CONWAYBROMAGELIB_PREVIOUSANDNEXT_H
