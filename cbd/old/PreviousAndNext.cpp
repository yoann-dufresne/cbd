//
// Created by alexa on 18/06/2020.
//

/**
 * Returns the successors of a kmer.
 * @param nonCompressedKMer.
 * @param currentCompressedSeq.
 * @return a vector<uint64_t> representing the (at most 4) successors of the kmer.
 */
std::vector<uint64_t> next(uint64_t nonCompressedKMer, sdsl::sd_vector<> const& currentCompressedSeq){
    vector<uint64_t> res;
    uint64_t Pmer_size_in_sdvector = (int)(log(currentCompressedSeq.size())/log(4));
    //we must have nonCompressedKmer < 4^(P-1)
    uint64_t limit = pow(4,Pmer_size_in_sdvector-1);
    if(nonCompressedKMer >= limit) {
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << limit << endl;
        return res;
    }

    uint64_t i_pmer1 = nonCompressedKMer << 2; //<==> AAC << 2 ==> AACA
    uint64_t i_kmer1 = i_pmer1%(currentCompressedSeq.size() >> 2);
    //if currentCompressedSeq.size() = 256 (P=4) then currentCompressedSeq.size() >> 2 = 64
    //<==> AACA%64 ==> ACA (i_kmer1)
    if(currentCompressedSeq[i_pmer1])   res.push_back(i_kmer1);
    if(currentCompressedSeq[i_pmer1+1]) res.push_back(i_kmer1+1);
    if(currentCompressedSeq[i_pmer1+2]) res.push_back(i_kmer1+2);
    if(currentCompressedSeq[i_pmer1+3]) res.push_back(i_kmer1+3);
    return res;
}

/* Look for previous version of a given (k-1)-mer in a generated sequence
 * @param compressedKMer - a uint64_t which represent a compressed (k-1)-mer. We can find it in the generated sequence
 * @param currentCompressedSeq - the generated sequence which have to contains compressedKMer and its potantial previous (k-1)-mers
 * @return a vector of uint64_t which contains all the previous (k-1)-mers which are present in the generated sequence, compressed form
 */
vector<uint64_t> previous(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq){
    vector<uint64_t>prev;
    if(compressedKMer < currentCompressedSeq.size()/4 && compressedKMer >= 0){
        int currentKMerLen = log(currentCompressedSeq.size()) / log(ALPHABET);
        string sub = decode(compressedKMer, currentKMerLen-1);
        uint64_t potentialPrevious[4];
        potentialPrevious[0] = compressedKMer % currentCompressedSeq.size();
        for(int i = 1 ; i < 4 ; i++){
            potentialPrevious[i] = potentialPrevious[i-1] + (currentCompressedSeq.size() / 4);
        }
        for(int i = 0 ; i < 4 ; i++){
            if(currentCompressedSeq[potentialPrevious[i]]){
                uint64_t subPrev = encode(decode(potentialPrevious[i], currentKMerLen).erase(currentKMerLen-1, 1), currentKMerLen-1);
                prev.push_back(subPrev);
            }
        }
    }
    return prev;
}