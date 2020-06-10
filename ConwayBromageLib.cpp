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

/* previous version for already compressed KMer
 * @param compressedKMer - a uint64_t which represent a compressed k-mer. We can find it in the generated sequence
 * @param currentCompressedSeq - the generated sequence which have to contains nonCompressedKMer and its potantial previous K-mers
 * @return a vector of uint64_t which contains all the previous k-mers which are present in the generated sequence, compressed form
 */
vector<uint64_t> previous(uint64_t compressedKMer, sdsl::sd_vector<> const& currentCompressedSeq){
    vector<uint64_t>prev;
    if(compressedKMer < currentCompressedSeq.size() && compressedKMer >= 0){
        uint64_t potentialPrevious[4];
        potentialPrevious[0] = (compressedKMer >> 2)%currentCompressedSeq.size();
        for(int i = 1 ; i < 4 ; i++){
            potentialPrevious[i] = potentialPrevious[i-1] + (currentCompressedSeq.size() / 4);
        }
        for(int i = 0 ; i < 4 ; i++){
            prev.push_back(potentialPrevious[i]);
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
        long int myTotalLen(pow(ALPHABET,myWordLen));   //Creation of the total length to create the sd_vector_builder
        cout << "Total length : " << myTotalLen << endl;
        sd_vector_builder constructSparse(myTotalLen, myOneLen);    //A size of myTotalLen, contains myOneLen ones
        while(file >> word){
            cout << "The seq is : " << word << endl;
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
 * @param nonCompressedKMer - a string which represents the k-mer we study
 * @param currentCompressedSeq - the generated sequence where we want to verify if nonCompressedKMer is in
 * @return true if the givenl-mer is present, else false
 */
bool isThisKMerHere(uint64_t compressedKMer, sd_vector<> const& currentCompressedSeq){
    if(compressedKMer < currentCompressedSeq.size()){   //call of isTheSameSize to verify the size
        int currentKMerLen = log(currentCompressedSeq.size()) / log(ALPHABET);
        if(currentCompressedSeq[compressedKMer]){  //verify if the case is set to one
            cout << decode(compressedKMer, currentKMerLen) << " is present" << endl;
            return true;
        }else{
            cout << decode(compressedKMer, currentKMerLen) << " is absent" << endl;
            return false;
        }
    }
    cout << "Sequence size problem" << endl;
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
        cout << "IF : " << endl;
        int position(len-complem);
        final = position-1;
        cout << final << endl;
    }else{
        cout << "ELSE : " << endl;
        int position(complem);
        final = len-position-1;
        cout << final << endl;
    }
    return final;
}

/* Transform sequences which are contain in a file in a sd_vector
 * Also give the reverse complement in an other sd_vector, for ecoli only
 * @param path - a string which is the path to the file which contains the generated sequence
 * @return a vector of sd_vector which contains encoding version of the sequence of the file and of the reverse complement of the sequence
 */
vector<sd_vector<>>fromFileToSdVectorWithReverseEcoli(string path){
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
                constructSparse.set(encodeEcoli(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
                reverser[i] = reverseComplement(word, myTotalLen);  //each reverses of each words, compressed version
                i++;
                file >> word;
        }
        sort(reverser.begin(), reverser.end()); //sorting of the reverse otherwise it will be impossible to use sd_vector_builder
        for(int i = 0 ; i < myOneLen ; i++){
            constructReverse.set(reverser[i]);  //filled for the reverse version
        }
        vector<sd_vector<>> seq;
        sd_vector<>finalSparseRight(constructSparse);    //Construction of the final sd_vector for the right sequence
        sd_vector<>finalSparseReverse(constructReverse);   //Construction of the final sd_vector for the reverse sequence
        seq.push_back(finalSparseRight);
        seq.push_back(finalSparseReverse);
        file.close();
        return seq;
    }else{
        cout << "Error while opening" << endl;
    }
    //return bit_vector{0}; //need an other return here
}
