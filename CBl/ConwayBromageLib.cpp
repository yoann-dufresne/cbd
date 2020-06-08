//
// Created by Alexandra and Murat on 05/06/2020.
//

#include <iostream>
#include <fstream>
#include <vector>
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
 * @return a vector<string> representing the (at most 4) successors of the kmer.
 */
vector<string> next(string nonCompressedKMer, sd_vector<> const& currentCompressedSeq){
    vector<string> res;
    //We check if the nonCompressedKMer exists (i.e set to 1) in the currentCompressedSeq.
    //It exists if it has a proper size and if his index is set to 1 in the sd_vector.
    uint64_t K = (int)(log(currentCompressedSeq.size())/log(4));
    if(nonCompressedKMer.size() != K)
        return res;
    uint64_t index = encode(nonCompressedKMer, nonCompressedKMer.size());
    if(!currentCompressedSeq[index]) //if set to 0, there are no successors
        return res;
    
    int method = 1;
    if(method == 1){ //METHOD1 for k = 9 and n = 262144, the function returns the successors of all kmer in txt-file in 45993 ms
        uint64_t lastIndexInString = K-1;
        string kmer = nonCompressedKMer.substr(1, lastIndexInString) + "A";
        //next(n) = {4n%size, (4n+1)%size, (4n+2)%size, (4n+3)%size}
        uint64_t i_kmer1 = (index << 2)%currentCompressedSeq.size();
        uint64_t i_kmer2 = i_kmer1+1;
        uint64_t i_kmer3 = i_kmer1+2;
        uint64_t i_kmer4 = i_kmer1+3;
        if(index != i_kmer1 && currentCompressedSeq[i_kmer1]) res.push_back(kmer);
        kmer[lastIndexInString] = 'C';
        if(index != i_kmer2 && currentCompressedSeq[i_kmer2]) res.push_back(kmer);
        kmer[lastIndexInString] = 'G';
        if(index != i_kmer3 && currentCompressedSeq[i_kmer3]) res.push_back(kmer);
        kmer[lastIndexInString] = 'T';
        if(index != i_kmer4 && currentCompressedSeq[i_kmer4]) res.push_back(kmer);
    } else if(method == 2) { //METHOD2 for k = 9 and n = 262144, the function returns the successors of all kmer in txt-file in 43726 ms
        uint64_t i_kmer1 = (index << 2)%currentCompressedSeq.size();
        uint64_t i_kmer2 = i_kmer1+1;
        uint64_t i_kmer3 = i_kmer1+2;
        uint64_t i_kmer4 = i_kmer1+3;
        if(index != i_kmer1 && currentCompressedSeq[i_kmer1]) res.push_back(decode(i_kmer1, K));
        if(index != i_kmer2 && currentCompressedSeq[i_kmer2]) res.push_back(decode(i_kmer2, K));
        if(index != i_kmer3 && currentCompressedSeq[i_kmer3]) res.push_back(decode(i_kmer3, K));
        if(index != i_kmer4 && currentCompressedSeq[i_kmer4]) res.push_back(decode(i_kmer4, K));
    } else { //METHOD 3 for k = 9 and n = 262144, the function returns the successors of all kmer in txt-file in 45807 ms
        std::bitset<64> i = encode(nonCompressedKMer.substr(1, K-1), K-1);
        std::bitset<66> kmer1(i.to_string()+"00");
        uint64_t i_kmer1 = static_cast<uint64_t>(kmer1.to_ulong());
        uint64_t i_kmer2 = i_kmer1+1;
        uint64_t i_kmer3 = i_kmer1+2;
        uint64_t i_kmer4 = i_kmer1+3;
        if(index != i_kmer1 && currentCompressedSeq[i_kmer1]) res.push_back(decode(i_kmer1, K));
        if(index != i_kmer2 && currentCompressedSeq[i_kmer2]) res.push_back(decode(i_kmer2, K));
        if(index != i_kmer3 && currentCompressedSeq[i_kmer3]) res.push_back(decode(i_kmer3, K));
        if(index != i_kmer4 && currentCompressedSeq[i_kmer4]) res.push_back(decode(i_kmer4, K));
    }
    return res;
}

/* Verifiy the size of 2 k-mer
 * if size is equal, return true, else return false
 */
bool isTheSameSize(int KMerLen, int currentCompressedSeqLen){
    int currentKMerLen = log(currentCompressedSeqLen) / log(ALPHABET);  //Original size of k-mers of the compressed sequence
    if(KMerLen != currentKMerLen){
        //cout << "Invalidated : Your sequence is a " << KMerLen << "-mers, we need a " << currentKMerLen << "-mers" << endl;
        return false;
    }else{
        return true;
    }
}

/* Look for the previous k-mer for a given k-mer in the generated sequence
 * We can have at most 4 different previous k-mers
 * We can use previous only if the given k-mer (the one for we search for previous elements) is present in the sequence
 * (isThisKMerHere is true) and its size is equal to the size of sequence k-mers (isTheSameSize is true)
 * returns a vector which contains all the previous k-mers which are present in the generated sequence
 */
vector<string> previous(std::string nonCompressedKMer, sdsl::sd_vector<> const& currentCompressedSeq){
    string potentialPrevious [4] = {"", "", "", ""};    //At most 4 potential previous k-mers
    int potentialCompressed [4] = {0, 0, 0, 0};         //Compressed version of the above array
    vector<string> prev;
    //verify same size and existence in the sequence
    if(isTheSameSize(nonCompressedKMer.size(), currentCompressedSeq.size()) && isThisKMerHere(nonCompressedKMer, currentCompressedSeq)) {
        for(int i = 0 ; i < 4 ; i++){
            //building of the 4 potential previous k-mers
            potentialPrevious[i] = NUCLEOTIDES[i];
            for(int j = 0 ; j < nonCompressedKMer.size()-1 ; j++){
                potentialPrevious[i] = potentialPrevious[i] + nonCompressedKMer[j];
            }
        }
        for(int i = 0 ; i < 4 ; i++){
            potentialCompressed[i] = encode(potentialPrevious[i], nonCompressedKMer.size());    //compressed version of the potential previous
            //If it is set to one in the compressed sequence, push it in the final vector
            if(currentCompressedSeq[potentialCompressed[i]] == 1){
                prev.push_back(potentialPrevious[i]);
            }
        }
    }
    return prev;
}

/* Transform sequences which are contain in a file in a sd_vector
 * Need a string which is the path to the file
 * Return a sd_vector which contains encoding version of sequences of the file
 */
sd_vector<>fromFileToSdVector(string path){
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
        long int myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create the sd_vector_builder
        cout << "Total length : " << myTotalLen << endl;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones
        while(file >> word){
            if(word != "1"){
                //cout << "The seq is : " << word << endl;
                constructSparse.set(encode(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
            }
        }
        sd_vector<>finalSparse(constructSparse);    //Construction of the final sd_vector
        file.close();
        return finalSparse;
    }else{
        cout << "Error while opening" << endl;
    }
    return bit_vector{0};
}

/* Verify if a given k-mer is present in the generated sequence
 * If the given k-mer and generated sequence members have not the same size, comparison is impossible
 * it returns true if the given k-mer is present, false if it is absent
 */
bool isThisKMerHere(std::string nonCompressedKMer, sdsl::sd_vector<> const& currentCompressedSeq){
    if(isTheSameSize(nonCompressedKMer.size(), currentCompressedSeq.size())){   //call of isTheSameSize to verify the size
        uint64_t myEncodingKMer = encode(nonCompressedKMer, nonCompressedKMer.size());  //encoding version of the K-mer
        if(currentCompressedSeq[myEncodingKMer] == 1){  //verify if the case is set to one
            //cout << nonCompressedKMer << " is present" << endl;
            return true;
        }else{
            //cout << nonCompressedKMer << " is absent" << endl;
            return false;
        }
    }
    //cout << "Sequence size problem" << endl;
    return false;
}
