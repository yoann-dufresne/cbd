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
//Switch version
uint64_t encodeEcoliSwitchVers(string word, uint64_t size){
    uint64_t hash = 0;
    char c;
    for(uint i = 0 ; i < size ; i++){   //We go through the sequence to encode each nucleotides
        hash <<= 2; // We shift 2 bits to the right
        c = word[i];    //Take each nucleotides to encode them
        uint64_t charval = 0; // 'A' = 0
        switch(c){
            case 'C': charval = 1; break;
            case 'T': charval = 2; break;
            case 'G': charval = 3; break;
        }
        hash += charval;    //creation of the hash for the given sequence
    }
    return hash;    //return the final hash of the sequence
}


//Switch version
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
//If version
string decodeEcoliIfVers(uint64_t seq, uint64_t size){
    string res(size, ' ');
    uint64_t lastIndex = res.size()-1;
    for(int i(0) ; i < size ; i++){
        res[lastIndex-i] = 'A';
        if((seq & 0x3) == 1){
            res[lastIndex-i] = 'C';
        }else if ((seq & 0x3) == 2){
            res[lastIndex-i] = 'T';
        }else if ((seq & 0x3) == 3){
            res[lastIndex-i] = 'G';
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
 * Return the canonical form of a kmer according to a specified encoding.
 * @param kmer
 * @param kmerSize - Size of the kmer
 * @param encodingIsACGT - true if the wanted encoding is ACGT. If it's false, then the encoding ACTG will be considered.
 * @return canonical version of the kmer
 */
uint64_t getCanonical (uint64_t kmer, uint64_t kmerSize, bool encodingIsACGT){
    if(encodingIsACGT){ //ACGT encoding
        uint64_t reverseComplement = reverseComplementLexico(kmer, kmerSize);
        return (kmer < reverseComplement)?kmer:reverseComplement;
    }
    //ACTG encoding
    uint64_t reverseComplement = reverseComplementGATBLibEcoli(kmer, kmerSize);
    return (kmer < reverseComplement)?kmer:reverseComplement;
}

/**
 * Returns the successors of a canonical kmer.
 * @param Kmer : the Kmer that we want to find its successors. Can either be canonical or not.
 * @param compressedSeq : the sd_vector which stores all the canonical kmers.
 * @param encodingIsACGT : true if the wanted encoding is ACGT. If it's false, then the encoding ACTG will be considered.
 * @return a vector<uint64_t> representing the (at most 8) successors of the kmer. The successors can be NON CANONICAL.
 */
vector<uint64_t> successors(uint64_t Kmer, sd_vector<> const& compressedSeq, bool encodingIsACGT){
    vector<uint64_t> res;
    int PmerSize = (int)(log(compressedSeq.size())/log(4));
    int KmerSize = PmerSize-1;
    uint64_t limit = compressedSeq.size() >> 2;
    if(Kmer >= limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << limit << endl;
        return res; //empty
    }  
    
    //we retrieve the reverse complement of the Kmer
    uint64_t KmerRevComp = (encodingIsACGT)?reverseComplementLexico(Kmer, KmerSize):reverseComplementGATBLibEcoli(Kmer, KmerSize);
    
    int numberOfBitsToShift = KmerSize << 1; //<-> 2*KmerSize
    uint64_t one = 1 << numberOfBitsToShift;
    uint64_t nextForward = Kmer << 2; //the forward version of the current pmer that we are interested in
    uint64_t nextRevComp = KmerRevComp + (3 << numberOfBitsToShift); //the reverse complement of the current pmer that we are interested in
    for(int i = 0; i < 4; i++){
        if(!encodingIsACGT) nextRevComp = KmerRevComp + (((i<2)?i+2:i%2) << numberOfBitsToShift);
        
        uint64_t canonical = (nextForward < nextRevComp)?nextForward:nextRevComp; 
        if(compressedSeq[canonical]){ //we look if the canonical version of the pmer is present
            uint64_t successorKmer = nextForward%limit;
            if(find(res.begin(), res.end(), successorKmer) == res.end()) //check if vector doesn't contain successorKmer
                res.push_back(successorKmer);
        }
        
        if(encodingIsACGT) nextRevComp -= one;
        nextForward++;
    }
    
    //predecessors
    uint64_t prevForward = Kmer; //the forward version of the current pmer that we are interested in
    uint64_t prevRevComp = (KmerRevComp << 2) + 3; //the reverse complement of the current pmer that we are interested in
    for(int i = 0; i < 4; i++){
        if(!encodingIsACGT) prevRevComp = (KmerRevComp << 2) + ((i<2)?i+2:i%2);

        uint64_t canonical = (prevForward < prevRevComp)?prevForward:prevRevComp;
        if(compressedSeq[canonical]){ //we look if the canonical version of the pmer is present
            uint64_t successorKmer = prevForward >> 2;
            if(find(res.begin(), res.end(), successorKmer) == res.end()) //check if vector doesn't contain successorKmer
                res.push_back(successorKmer);
        }
        prevForward += one; //equals to AX, CX, GX, then TX where X is the Kmer
        if(encodingIsACGT) prevRevComp--;
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
/* Transform sequences which are contain in a file in a sd_vector
 * @param path - The path to the file which contains info
 * @param format - The format we want for the transformation : ACGT or ACTG
 * @return a sd_vector which contains elements of the file transformed according to the format
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
        cout << "encoding in ACGT format... " << endl;
        while(file >> word){
            if(word != "1") {
                constructSparse.set(encode(word,myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
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

/* Transform sequences which are contain in a file in a sd_vector
 * @param path - The path to the file which contains info
 * @param format - The format we want for the transformation : ACGT or ACTG
 * @return a sd_vector which contains elements of the file transformed according to the format
 */
sd_vector<>fromFileToSdVectorChooser(string path, string format){
    ifstream file(path, ios::in);  // Reading of the file which contains k-mers sequences
    if(file){   // File is now open
        if(format == "ACGT" || format == "ACTG"){
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
            if(format == "ACGT"){
                cout << "encoding in ACGT format... " << endl;
                sd_vector_builder constructACGT(myTotalLen, myOneLen);
                while(file >> word){
                    uint64_t pmer = encode(word,myWordLen);
                    if(getCanonical(pmer, word.size(), true) != pmer){
                        cout << "The file is not completely canonical" << endl;
                        exit(1); //EXIT_FAILURE
                    }
                    constructACGT.set(pmer); //filled to one each element which is represent by the encoding version of the sequence
                    file >> word;
                }
                sd_vector<>finalACGT(constructACGT);
                file.close();
                return finalACGT;
            }else{
                cout << "encoding in ACTG format... " << endl;
                sd_vector_builder constructACTG(myTotalLen, myOneLen);
                while(file >> word){
                    uint64_t pmer = encodeEcoli(word,myWordLen);
                    if(getCanonical(pmer, word.size(), false) != pmer){
                        cout << "The file is not completely canonical" << endl;
                        exit(1); //EXIT_FAILURE
                    }
                    constructACTG.set(pmer);
                    file >> word;
                }
                sd_vector<>finalACTG(constructACTG);
                file.close();
                return finalACTG;
            }
        }else {
            cout << "invalide format : need ACGT or ACTG" << endl;
        }
    }else{
        cout << "Error while opening" << endl;
    }
    return bit_vector{0};
}

/**
 * Check if the given Kmer is present in the sequence.
 * @param Kmer - An integer representing the k-mer.
 * @param compressedSeq - The compressed p-mer sequence.
 * @param encodingIsACGT - true if the encoding is ACGT and false if it's ACTG.
 * @return true if the k-mer is present.
 */
bool isThisKMerHere(uint64_t Kmer, sd_vector<> const& compressedSeq, bool encodingIsACGT){
    int PmerSize = (int)(log(compressedSeq.size())/log(4));
    int KmerSize = PmerSize-1;
    uint64_t limit = compressedSeq.size() >> 2;
    if(Kmer >= limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << limit << endl;
        return false; 
    }  
    
    uint64_t KmerRevComp = (encodingIsACGT)?reverseComplementLexico(Kmer, KmerSize):reverseComplementGATBLibEcoli(Kmer, KmerSize);
    int numberOfBitsToShift = KmerSize << 1;
    uint64_t one = 1 << numberOfBitsToShift;
    uint64_t nextForward = Kmer << 2;
    uint64_t nextRevComp;
    
    nextRevComp = KmerRevComp + (((encodingIsACGT)?3:2) << numberOfBitsToShift);
    if(compressedSeq[(nextForward < nextRevComp)?nextForward:nextRevComp]) return true;
    
    nextForward++;
    nextRevComp = KmerRevComp + (((encodingIsACGT)?2:3) << numberOfBitsToShift);
    if(compressedSeq[(nextForward < nextRevComp)?nextForward:nextRevComp]) return true;
    
    nextForward++;
    nextRevComp = KmerRevComp + (((encodingIsACGT)?1:0) << numberOfBitsToShift);
    if(compressedSeq[(nextForward < nextRevComp)?nextForward:nextRevComp]) return true;
    
    nextForward++;
    nextRevComp = KmerRevComp + (((encodingIsACGT)?0:1) << numberOfBitsToShift);
    if(compressedSeq[(nextForward < nextRevComp)?nextForward:nextRevComp]) return true;
    
    uint64_t prevForward = Kmer;
    uint64_t prevRevComp = KmerRevComp << 2;
    
    prevRevComp = KmerRevComp + ((encodingIsACGT)?3:2);
    if(compressedSeq[(prevForward < prevRevComp)?prevForward:prevRevComp]) return true;
    
    prevForward += one; 
    prevRevComp = KmerRevComp + ((encodingIsACGT)?2:3);
    if(compressedSeq[(prevForward < prevRevComp)?prevForward:prevRevComp]) return true;
    
    prevForward += one;
    prevRevComp = KmerRevComp + ((encodingIsACGT)?1:0);
    if(compressedSeq[(prevForward < prevRevComp)?prevForward:prevRevComp]) return true;
    
    prevForward += one;
    prevRevComp = KmerRevComp + ((encodingIsACGT)?0:1);
    if(compressedSeq[(prevForward < prevRevComp)?prevForward:prevRevComp]) return true;
    
    return false;
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

/* Successor counter to verify veracity of successors function
 * It is a non optimized function (slow), we use it for tests exclusively
 * @param compressedKmer - The k-mer for which we want to know successors
 * @param currentCompressedSeq - The given which is supposed to contain successors of compressedKer
 * @param format - a boolean to konw if we use ACGT or ACTG encoding
 * @return a vector which contain successors of compressedKMer
 */
vector<uint64_t> successorCounter(uint64_t compressedKMer, sd_vector<>currentCompressedSeq, bool format){
    uint64_t sizeOfSeq = log(currentCompressedSeq.size()) / log(ALPHABET);
    uint64_t sizeOfKMer = sizeOfSeq - 1;
    uint64_t list[8];
    string istr;
    vector<uint64_t> reader;
    if(format){
        istr = decode(compressedKMer, sizeOfKMer);
        list[0] = encode(istr + "A", sizeOfSeq);
        list[1] = encode(istr + "C", sizeOfSeq);
        list[2] = encode(istr + "G", sizeOfSeq);
        list[3] = encode(istr + "T", sizeOfSeq);
        list[4] = encode("A" + istr, sizeOfSeq);
        list[5] = encode("C" + istr, sizeOfSeq);
        list[6] = encode("G" + istr, sizeOfSeq);
        list[7] = encode("T" + istr, sizeOfSeq);
    }else{
        istr = decodeEcoli(compressedKMer, sizeOfKMer);
        list[0] = encodeEcoli(istr + "A", sizeOfSeq);
        list[1] = encodeEcoli(istr + "C", sizeOfSeq);
        list[2] = encodeEcoli(istr + "G", sizeOfSeq);
        list[3] = encodeEcoli(istr + "T", sizeOfSeq);
        list[4] = encodeEcoli("A" + istr, sizeOfSeq);
        list[5] = encodeEcoli("C" + istr, sizeOfSeq);
        list[6] = encodeEcoli("G" + istr, sizeOfSeq);
        list[7] = encodeEcoli("T" + istr, sizeOfSeq);
    }
    for(int i = 0 ; i < 8 ; i++){
        uint64_t lecteur;
        if(!isCanonical(list[i], sizeOfSeq, format)){
            //cout << "not cano : " << decodeEcoli(list[i], 3) << endl;
            uint64_t canoVers = getCanonical(list[i], sizeOfSeq, format);
            //cout << "cano version : " << canoVers << endl;
            if(currentCompressedSeq[canoVers]){
                //cout << "PASSED" << endl;
                if(format){
                    if(i > 3){
                        lecteur = encode(decode(list[i], sizeOfSeq).erase(sizeOfSeq-1, 1), sizeOfKMer);
                    }else{
                        lecteur = encode(decode(list[i], sizeOfSeq).erase(0, 1), sizeOfKMer);
                    }
                }else{
                    if(i > 3){
                        lecteur = encodeEcoli(decodeEcoli(list[i], sizeOfSeq).erase(sizeOfSeq-1, 1), sizeOfKMer);
                    }else{
                        lecteur = encodeEcoli(decodeEcoli(list[i], sizeOfSeq).erase(0, 1), sizeOfKMer);
                    }
                }
                //cout << "lecteur : " << lecteur << endl;
                if(find(reader.begin(), reader.end(), lecteur) == reader.end()){
                   reader.push_back(lecteur);
               }
            }
        }else{
            if(currentCompressedSeq[list[i]]){
                if(format){
                    if(i > 3){
                        lecteur = encode(decode(list[i], sizeOfSeq).erase(sizeOfSeq-1, 1), sizeOfKMer);
                    }else{
                        lecteur = encode(decode(list[i], sizeOfSeq).erase(0, 1), sizeOfKMer);
                    }
                }else{
                    if(i > 3){
                        lecteur = encodeEcoli(decodeEcoli(list[i], sizeOfSeq).erase(sizeOfSeq-1, 1), sizeOfKMer);
                    }else{
                        lecteur = encodeEcoli(decodeEcoli(list[i], sizeOfSeq).erase(0, 1), sizeOfKMer);
                    }
                }
                //cout << "lecteur : "<< lecteur << endl;
                if(find(reader.begin(), reader.end(), lecteur) == reader.end()){
                   reader.push_back(lecteur);
               }
            }
        }
    }
    return reader;
}

/* Variation of getCanonical : just tell if it is canonical or not
 * Just use for tests
 * @param kmer - the k-mer we want to test
 * @param kmerSize - the size if this k-mer
 * @param encodingIsACGT - boolean to know if encoding is ACGT or ACTG
 * @return true if it is canonical, else false
 */
bool isCanonical (uint64_t kmer, uint64_t kmerSize, bool encodingIsACGT){
    if(encodingIsACGT){ //ACGT encoding
        uint64_t reverseComplement = reverseComplementLexico(kmer, kmerSize);
        return (kmer < reverseComplement)?true:false;
    }
    //ACTG encoding
    uint64_t reverseComplement = reverseComplementGATBLibEcoli(kmer, kmerSize);
    return (kmer < reverseComplement)?true:false;
}

