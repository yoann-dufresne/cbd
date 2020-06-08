//
// Created by Alexandra and Murat on 05/06/2020.
//

#include <iostream>
#include <fstream>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include "ConwayBromageLib.h"

using namespace std;
using namespace sdsl;

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
    //cout << hash << endl;   // print of verification : hash building
    return hash;    //return the final hash of the sequence
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
                //cout << "The seq is : " << word << endl;  //test
                constructSparse.set(encode(word, myWordLen)); //filled to one each element which is represent by the encoding version of the sequence
            }
        }
        sd_vector<>finalSparse(constructSparse);    //Construction of the final sd_vector
        //tests
        /*sd_vector<>::iterator itSparse;
        cout << "value of the sd_vector : " << endl;
        for(itSparse = finalSparse.begin() ; itSparse != finalSparse.end() ; itSparse++){
          cout << *itSparse << " ";
        }*/
        cout << endl;
        file.close();
        return finalSparse;
    }else{
        cout << "Error while opening" << endl;
    }
}

/* Verifiy the size of 2 k-mer
 * if size is equal, return true, else return false
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

/* Verify if a given k-mer is present in the generated sequence
 * If the given k-mer and generated sequence members have not the same size, comparison is impossible
 * it returns true if the given k-mer is present, false if it is absent
 */
bool isThisKMerHere(std::string nonCompressedKMer, sdsl::sd_vector<> const& currentCompressedSeq){
    if(isTheSameSize(nonCompressedKMer.size(), currentCompressedSeq.size())){   //call of isTheSameSize to verify the size
        uint64_t myEncodingKMer = encode(nonCompressedKMer, nonCompressedKMer.size());  //encoding version of the K-mer
        if(currentCompressedSeq[myEncodingKMer] == 1){  //verify if the case is set to one
            cout << nonCompressedKMer << " is present" << endl;
            return true;
        }else{
            cout << nonCompressedKMer << " is absent" << endl;
            return false;
        }
    }
    cout << "Sequence size problem" << endl;
    return false;
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

