//
// Created by Alexandra and Murat on 05/06/2020.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include "ConwayBromageLib.h"

using namespace std;
using namespace sdsl;

/**
 * Returns the number of 1-bit to the left of position index
 * @param index the position
 * @param v an sd_vector
 * @return a count
 */
uint64_t rank1bit(uint64_t index, const sd_vector<> &v){
    sd_vector<>::rank_1_type sdb_rank(&v);
    return sdb_rank(index);
}

/**
 * Returns the number of 0-bit to the left of position index
 * @param index the position 
 * @param v an sd_vector
 * @return a count
 */
uint64_t rank0bit(uint64_t index, const sd_vector<> &v){
    sd_vector<>::rank_0_type sdb_rank(&v);
    return sdb_rank(index);
}

/**
 * Returns the position of the index-th 1-bit.
 * @param index the position 
 * @param v an sd_vector
 * @return an index
 */
uint64_t select1bit(uint64_t index, const sd_vector<> &v){
    sd_vector<>::select_1_type sdb_sel(&v);
    return sdb_sel(index);
}

/**
 * Returns the position of the index-th 0-bit.
 * @param index the position 
 * @param v an sd_vector
 * @return an index
 */
uint64_t select0bit(uint64_t index, const sd_vector<> &v){
    sd_vector<>::select_0_type sdb_sel(&v);
    return sdb_sel(index);
}

/**
 * The ecoli_count.txt file. Sort the k-mers with this order : A > C > T > G.
 * @param word the k-mer
 * @param size the size of the k-mer
 * @return an integer representing the k-mer
 */
uint64_t encodeEcoli(string word, uint64_t size){
    uint64_t hash = 0;
    char c;
    for(uint i = 0 ; i < size ; i++){   //We go through the sequence to encode each nucleotides
        hash <<= 2; // We shift 2 bits to the right
        c = word[i];    //Take each nucleotides to encode them
        uint64_t charval = 0; // 'A' = 0
        if(c == 'C'){
            charval = 1;    // 'C' = 1
        }else if(c == 'T'){
            charval = 2;    // 'T' = 2
        }else if(c == 'G'){
            charval = 3;    // 'G' = 3
        }
        hash += charval;    //creation of the hash for the given sequence
    }
    return hash;    //return the final hash of the sequence
}


string decodeEcoli(uint64_t seq, uint64_t size){
    string res(size, ' ');
    uint64_t lastIndex = res.size()-1;
    for(int i(0); i < size; i++){
        switch(seq & 0x3){ //compares the decimal value of the first two bits
            case 0: res[lastIndex-i] = 'A'; break;
            case 1: res[lastIndex-i] = 'C'; break;
            case 2: res[lastIndex-i] = 'T'; break;
            case 3: res[lastIndex-i] = 'G'; break;
        }
        seq >>= 2;
    }
    return res;
}

/* function to encode k-mer sequences
 * Need a string (the sequence) and a uint64_t which is the size of the sequence
 * return a uint64_t which is the encoding version of the sequence
 * We follow the following encoding :
 * 'A' = 0 ; 'C' = 1 ; 'G' = 2 : 'T' = 3
 */
uint64_t encode(string word, uint64_t size){
    uint64_t hash = 0;
    char c;
    for(uint i = 0 ; i < size ; i++){   //We go through the sequence to encode each nucleotides
        hash <<= 2; // We shift 2 bits to the right
        c = word[i];    //Take each nucleotides to encode them
        uint64_t charval = 0; // 'A' = 0
        if(c == 'C'){
            charval = 1;    // 'C' = 1
        }else if(c == 'G'){
            charval = 2;    // 'G' = 2
        }else if(c == 'T'){
            charval = 3;    // 'T' = 3
        }
        hash += charval;    //creation of the hash for the given sequence
    }
    return hash;    //return the final hash of the sequence
}
/**
 * Returns a k-mer which corresponds to the value of seq.
 * @param seq - the value.
 * @param size - the size of the k-mer.
 * @return a string representing the value seq.
 */
string decode(uint64_t seq, uint64_t size){
    string res(size, ' ');
    uint64_t lastIndex = res.size()-1;
    for(int i(0); i < size; i++){
        switch(seq & 0x3){ //compares the decimal value of the first two bits
            case 0: res[lastIndex-i] = 'A'; break;
            case 1: res[lastIndex-i] = 'C'; break;
            case 2: res[lastIndex-i] = 'G'; break;
            case 3: res[lastIndex-i] = 'T'; break;
        }
        seq >>= 2;
    }
    return res;
}

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

/* Verifiy the size of 2 k-mer
 * @param KMerLen - The size of the current K-mer (the one we study)
 * @param currentCompressedSeqLen - The size of the generated sequence
 * @retrurn true if the size of the current K-mer and K-mers of the sequence is equal, else false
 */
bool isTheSameSize(int KMerLen, int currentCompressedSeqLen){
    int currentKMerLen = log(currentCompressedSeqLen) / log(ALPHABET);  //Original size of k-mers of the compressed sequence
    if(KMerLen != currentKMerLen){
        cout << "Invalidated : Your sequence is a " << KMerLen << "-mers, we need a " << currentKMerLen << "-mers" << endl;
        return false;
    }else{
        return true;
    }
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

/* Transform sequences which are contain in a file in a sd_vector
 * Need a string which is the path to the file
 * Return a sd_vector which contains encoding version of sequences of the file
 */
sd_vector<>fromFileToSdVectorEcoli(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(file){   // File is now open
        string word;
        string line("");
        file >> word;   //Take the first word to analyze size of one k-mer
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        int myWordLen(word.size()); //Size of k_mer, it is the 'k'
        int myOneLen(0);
        while(getline(file, line)){ //Counts the number of ones in the file
            myOneLen++;
        }
        file.clear();
        file.seekg(0, ios::beg);    //Return to he beginning of the file
        cout << "length of ones : " << myOneLen << endl;
        cout << "length of a seq : " << myWordLen << endl;
        uint64_t myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create the sd_vector_builder
        cout << "Total length : " << myTotalLen << endl;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones
        while(file >> word){
            //cout << "The seq is : " << word << endl;
            constructSparse.set(encodeEcoli(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
            file >> word;
        }
        sd_vector<>finalSparse(constructSparse);    //Construction of the final sd_vector
        file.close();
        return finalSparse;
    }else{
        cout << "Error while opening" << endl;
    }
    return bit_vector{0};
}

/* Transform sequences which are contain in a file in a sd_vector
 * @param path - The path to the file which contains info
 * @param format - The format we want for the transformation : ACGT or ACTG
 * @return a sd_vector which contains elements of the file transformed according to the format
 */
sd_vector<>fromFileToSdVector(string path, string format){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(file && (format == "ACGT" || format == "ACTG")){   // File is now open
        string word;
        string line("");
        file >> word;   //Take the first word to analyze size of one k-mer
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        int myWordLen(word.size()); //Size of k_mer, it is the 'k'
        int myOneLen(0);
        while(getline(file, line)){ //Counts the number of ones in the file
            myOneLen++;
        }
        file.clear();
        file.seekg(0, ios::beg);    //Return to he beginning of the file
        cout << "length of ones : " << myOneLen << endl;
        cout << "length of a seq : " << myWordLen << endl;
        uint64_t myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create the sd_vector_builder
        cout << "Total length : " << myTotalLen << endl;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones
        if(format == "ACGT"){
            cout << "encoding in ACGT format... " << endl;
            while(file >> word){
                constructSparse.set(encode(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
                file >> word;
            }
        }else{
            cout << "encoding in ACTG format... " << endl;
            while(file >> word){
                constructSparse.set(encodeEcoli(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
                file >> word;
            }
        }
        sd_vector<>finalSparse(constructSparse);    //Construction of the final sd_vector
        file.close();
        return finalSparse;
    }else{
        cout << "Error while opening or bad format : need ACGT or ACTG" << endl;
    }
    return bit_vector{0};
}

/* Verify if a given (k-1)-mer is present in the generated sequence
 * If the given (k-1)-mer is bigger than the generated seuquence size /4, research is impossible
 * @param compressedKMer - a string which represents the (k-1)-mer we study
 * @param currentCompressedSeq - the generated sequence where we want to verify if compressedKMer is in
 * @return true if the given (k-1)-mer is present, else false
 */
bool isThisKMerHere(uint64_t compressedKMer, sd_vector<> const& currentCompressedSeq){
    if(compressedKMer < currentCompressedSeq.size() / 4 && compressedKMer >= 0) {
        int currentKMerLen = log(currentCompressedSeq.size()) / log(ALPHABET);
        string sub = decode(compressedKMer, currentKMerLen - 1);
        for (int i = 0; i < currentCompressedSeq.size(); i++) {
            if (currentCompressedSeq[i]) {
                if (decode(i, currentKMerLen).find(sub) != string::npos) {
                    return true;
                }
            }
        }
    }
    return false;
}

/* Calculate the reverse complement compressed version of a given k-mer
 * @param seq - the k-mer we want to know the reverse complement
 * @param len - the length of the total sequence of the given k-mer
 * @return a uint64_t which is the compressed version of the reverse complement
 */
uint64_t reverseComplement(string seq, uint64_t len){
    int sizeOfSeq = log(len) / log(ALPHABET);
    reverse(begin(seq), end(seq));  //Reverse of string
    uint64_t complem = encode(seq, sizeOfSeq);  //encoding the reverse
    uint64_t final;
    if(complem > len/2){    //Position for complement
        uint64_t position(len-complem);
        final = position-1;
    }else{
        uint64_t position(complem);
        final = len-position-1;
    }
    return final;
}

/* Calculate the reverse complement
 * Come from the wev page : https://www.biostars.org/p/113640/
 * Version for the lexicographical order : A = 00 ; C = 01 ; G = 10 ; T = 11
 * !!!!! Not the fastest version according to the webpage !!!!!
 * @param mer - a uin64_t which represent the compressed version of a k-mer
 * @param kmerSize - a uint64_t which represent the size of the given k-mer
 * @return a uin64_t which represent the compressed version of the reverse complement of the given k-mer (lexicographicalorder)
 */
uint64_t reverseComplementLexico (const uint64_t mer, uint64_t kmerSize)   //same with A = 00 ; C = 01 ; G = 10 ; T = 11
{
    uint64_t res = ~mer;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

    return (res >> (2 * (32 - kmerSize)));
}


/* Calculate the reverse complement
 * Come from the wev page : https://www.biostars.org/p/113640/
 * Version for the fastest ASCII order : A = 00 ; C = 01 ; T = 10 ; G = 11
 * @param x - a uin64_t which represent the compressed version of a k-mer
 * @param sizeKmer - a uint64_t which represent the size of the given k-mer
 * @return a uin64_t which represent the compressed version of the reverse complement of the given k-mer (lexicographicalorder)
 */
u_int64_t reverseComplementGATBLibEcoli (const u_int64_t x, uint64_t sizeKmer)      //GATB library edrezen  case A = 00 ; C = 01 ; G = 11 ; T = 10
{
    u_int64_t res = x;

    res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
    res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
    res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
    res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2*( 32 - sizeKmer))) ;
}

/*
 * Return the merging version of 2 sd_vector
 * @param a an sd_vector of the same size as b
 * @param b an sd_vector
 * @param nb_of_1_in_A the number of 1-bit in the sd_vector a
 * @param nb_of_1_in_B the number of 1-bit in the sd_vector b
 * @return a sd_vector which is the merging version of the 2 sd_vector params
 */
sd_vector<> merge (const sd_vector<> &a, const sd_vector<> &b, int nb_of_1_in_A, int nb_of_1_in_B){
    int len = a.size();
    int nb_of_1 = nb_of_1_in_A + nb_of_1_in_B;
    int nb_of_1_in_merged = 0;
    sd_vector_builder s(len, nb_of_1);
    //filling builder
    for(int i = 0; i < len; i++) {
        if (a[i] == 1 && b[i] == 1)
            s.set(i);
        if (a[i] == 1 || b[i] == 1) {
            s.set(i);
            nb_of_1_in_merged++;
        }
    }
    //creation of the sd_vector
    sd_vector<> merger(s);
    return merger;
}


/* Transform sequences which are contain in a file in a sd_vector
 * Merge with the reverse complement
 * @param path - a string which is the path to the file which contains the generated sequence
 * @return a sd_vector which contains an encoding merging version of the sequence of the file and of the reverse complement of the sequence
 */
sd_vector<>fromFileToSdVectorWithReverse(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(file){   // File is now open
        string word;
        string line("");
        file >> word;   //Take the first word to analyze size of one k-mer
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        uint64_t myWordLen(word.size()); //Size of k_mer, it is the 'k'
        uint64_t myOneLen(0);
        while(getline(file, line)){ //Counts the number of ones in the file
            myOneLen++;
        }
        file.clear();
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        cout << "length of ones : " << myOneLen << endl;
        cout << "length of a seq : " << myWordLen << endl;
        uint64_t myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create sd_vector_builders
        int_vector<> reverser(myOneLen, 0); //to simplify the reverse complement build
        cout << "Total length : " << myTotalLen << endl;
        int i = 0;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones, right version
        sd_vector_builder constructReverse(myTotalLen, myOneLen);   //same total size and ones size, reverse version
        while(file >> word){
            if(word != "1"){
                constructSparse.set(encode(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
                reverser[i] = reverseComplementLexico(encode(word, myWordLen), myWordLen);  //each reverses of each words, compressed version
                i++;
            }
        }
        sort(reverser.begin(), reverser.end()); //sorting of the reverse otherwise it will be impossible to use sd_vector_builder
        for(int i = 0 ; i < myOneLen ; i++){
            constructReverse.set(reverser[i]);  //filled for the reverse version
        }
        vector<sd_vector<>> seq;
        sd_vector<>finalSparseRight(constructSparse);    //Construction of the final sd_vector for the right sequence
        sd_vector<>finalSparseReverse(constructReverse);   //Construction of the final sd_vector for the reverse sequence
        sd_vector<>finalAll = merge(finalSparseRight, finalSparseReverse, finalSparseRight.size(), finalSparseReverse.size());  //merging
        file.close();
        return finalAll;
    }else{
        cout << "Error while opening" << endl;
    }
    return bit_vector{0};
}

/* Transform sequences which are contain in a file in a sd_vector
 * Merge the reverse complement, fastest ASCII version for reverse complement and encode/decode
 * @param path - a string which is the path to the file which contains the generated sequence
 * @return a sd_vector which contains an encoding merging version of the sequence of the file and of the reverse complement of the sequence
 */
sd_vector<>fromFileToSdVectorWithReverseEcoli(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(file){   // File is now open
        string word;
        string line("");
        file >> word;   //Take the first word to analyze size of one k-mer
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        uint64_t myWordLen(word.size()); //Size of k_mer, it is the 'k'
        uint64_t myOneLen(0);
        while(getline(file, line)){ //Counts the number of ones in the file
            myOneLen++;
        }
        file.clear();
        file.seekg(0, ios::beg);    //Return to the beginning of the file
        cout << "length of ones : " << myOneLen << endl;
        cout << "length of a seq : " << myWordLen << endl;
        uint64_t myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create sd_vector_builders
        int_vector<> reverser(myOneLen, 0); //to simplify the reverse complement build                                                              //SLOW DOWN
        cout << "Total length : " << myTotalLen << endl;
        int i = 0;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones, right version
        sd_vector_builder constructReverse(myTotalLen, myOneLen);   //same total size and ones size, reverse version
        while(file >> word){
                constructSparse.set(encodeEcoli(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
                reverser[i] = reverseComplementGATBLibEcoli(encodeEcoli(word, myWordLen), myWordLen);  //each reverses of each words, compressed version                                           //SLOW DOWN
                i++;
                file >> word;
        }
        sort(reverser.begin(), reverser.end()); //sorting of the reverse otherwise it will be impossible to use sd_vector_builder                               //SLOW DOWN
        for(int i = 0 ; i < myOneLen ; i++){
            constructReverse.set(reverser[i]);  //filled for the reverse version                                                                                //SLOW DOWN
        }
        vector<sd_vector<>> seq;
        sd_vector<>finalSparseRight(constructSparse);    //Construction of the final sd_vector for the right sequence
        sd_vector<>finalSparseReverse(constructReverse);   //Construction of the final sd_vector for the reverse sequence
        sd_vector<>finalAll = merge(finalSparseRight, finalSparseReverse, finalSparseRight.size(), finalSparseReverse.size());
        file.close();
        return finalAll;
    }else{
        cout << "Error while opening" << endl;
    }
    return bit_vector{0};
}

/**
 * Transform a number in base 10 in base 64.
 * @param valueInBase10
 * @return a string representing the number in base 64
 */
string convertToBase64(uint64_t valueInBase10){
    string res("");
    uint64_t biggest_divider = 1;
    while((uint64_t)(valueInBase10/biggest_divider) > 63) biggest_divider <<= 6;
    uint64_t remainder = valueInBase10;
    uint64_t quotient = valueInBase10;
    biggest_divider <<= 6; //*= 64 <-> 
    while(remainder != 0){
        biggest_divider >>= 6;
        quotient = remainder/biggest_divider;
        remainder = remainder % biggest_divider;
        res += (char)(';'+quotient);
    }
    if(biggest_divider > 1){
        string rest((int)(log(biggest_divider)/log(64)), ';');
        res += rest;
    }
    return res;
}

/**
 * Transform a number in base 10 in base 64.
 * @param valueInBase64
 * @return a number in base 10
 */
uint64_t convertFromBase64ToBase10(string valueInBase64){
    uint64_t res = 0;
    uint64_t pow = 1;
    for(int i = valueInBase64.size()-1; i >= 0; i--){
        res += (valueInBase64[i]-';') * pow;
        pow *= 64;
    }
    return res;
}

/**
 * Transform a file which contains k-mers in an sd_vector.
 * This version, generates a txt file which contains the values of reverse complement.
 * @param path of the file
 * @return an sd_vector of size 4^(P-1)
 */
sd_vector<> fromFileToSdVector_TXTversion(string path){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(!file) {
        cout << "Error while opening" << endl;
        return bit_vector{0};
    }

    //obtention of the k-mers' size
    string kmer("");
    file >> kmer;
    file.seekg(0, ios::beg);
    int K = kmer.size();
    
    //counting the number of k-mer in the file
    uint64_t nb_of_kmer = 0;
    string line("");
    while(getline(file, line)) nb_of_kmer++;
    cout << "nb_of_kmer : " << nb_of_kmer << endl;
    file.clear();
    file.seekg(0, ios::beg);
    sd_vector_builder sdv_builder_classic(pow(ALPHABET,K), nb_of_kmer);
    
    //first parsing to get the classic read
    while(file >> kmer){
        sdv_builder_classic.set(encode(kmer, K));//reverseComplementLexico(i_classic, K);
        file >> kmer;
    }
    file.clear();
    file.seekg(0, ios::beg);
    sd_vector<> classic(sdv_builder_classic);
    cout << "Length of classic : " << classic.size() << endl;
    cout << "size of classic (MB) : " << size_in_mega_bytes(classic) << endl;
    
    cout << "Completion of distinct_rev_comp.txt." << endl;
    auto b1 = std::chrono::high_resolution_clock::now();
    //second parsing to get the reverse complement
    ofstream output;
    output.open("distinct_rev_comp.txt");
    uint64_t nb_of_distinct_rev_comp = 0;
    while(file >> kmer){
        uint64_t i_classic = encode(kmer, K);
        uint64_t i_rc = reverseComplementLexico(i_classic, K);
        if(classic[i_rc]==0){ //if the reverse complement is not present in the classical reads then add it in the output
            output << i_rc;
            output << '\n';
            nb_of_distinct_rev_comp++;
        }
        file >> kmer;
    }
    output.close();
    auto e1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = e1 - b1;
    cout << "--> Time (s): " << elapsed1.count() << endl;
    
    //sort the output
    cout << "distinct_rev_comp.txt completed. Now we sort." << endl;
    auto b2 = std::chrono::high_resolution_clock::now();
    string command = "sort -n distinct_rev_comp.txt -o distinct_rev_comp.txt";
    system(command.c_str());
    cout << "Sort finished. " << endl;
    auto e2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = e2 - b2;
    cout << "--> Time (s): " << elapsed2.count() << endl;
    cout << "Now we build the final sd_vector." << endl;
    
    //now we build the sd_vector
    auto b3 = std::chrono::high_resolution_clock::now();
    ifstream input;
    input.open("distinct_rev_comp.txt");  // Reading of the file which contains k-mers sequences
    if(!input) {
        cout << "Error while opening distinct_rev_comp.txt" << endl;
        return bit_vector{0};
    };
    sd_vector_builder sdvB(pow(ALPHABET,K), nb_of_kmer + nb_of_distinct_rev_comp);
    uint64_t i = 0, i_classic = 1, i_rc = 0, last_rc = 0;
    
    sd_vector<>::select_1_type sdb_sel(&classic);
    uint64_t a = sdb_sel(i_classic);
    uint64_t b;
    input >> b;//get first value for b
    while(i_classic <= nb_of_kmer && i_rc < nb_of_distinct_rev_comp){
        if(a < b){
            sdvB.set(a);
            i_classic++;
            a = sdb_sel(i_classic);
        } else if (b < a){
            sdvB.set(b);
            i_rc++;
            if(i_rc < nb_of_distinct_rev_comp) input >> b;//update value
        } else { // a = b
            cout << "Should not happen because we have removed the duplicate precedently" << endl;
        }
    }
    cout << "i_classic after while loop: " << i_classic << endl;
    cout << "i_rc after while loop: " << i_rc << endl;
    //at this point, either i_classic = classic.size() or i_rc = rc.size()
    if(i_classic <= nb_of_kmer){
        while(i_classic <= nb_of_kmer){
            sdvB.set(sdb_sel(i_classic));
            i_classic++;
        }
    } else { //i_rc < rc.size()
        while(i_rc < nb_of_distinct_rev_comp){
            sdvB.set(b);
            i_rc++;
            input >> b;
        }
    }
    auto e3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed3 = e3 - b3;
    cout << "--> Time (s): " << elapsed3.count() << endl;
    
    cout << "i_classic at the end: " << i_classic << endl;
    cout << "i_rc at the end: " << i_rc << endl;
    input.close();
    sd_vector<> sdv(sdvB);
    return sdv;
}
