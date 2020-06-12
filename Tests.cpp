/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   Tests.cpp
 * Author: muratokutucu
 *
 * Created on 8 juin 2020, 11:07
 */

#include <cstdlib>
#include <string>
//#include "Tests.h"
#include "ConwayBromageLib.h"
#include <lest/lest_basic.hpp>        //Use of lest for unit tests, see here : https://github.com/martinmoene/lest.git
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>


using namespace std;
using namespace sdsl;
using namespace lest;


static const uint64_t  totalLen = 4611686018427387904;
static const sd_vector<>littleTestPrev = bit_vector{0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1};

/*
 * Compilation line : g++ -Wall -Wextra -std=c++11 -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -O3 -DNDEBUG -I ~/include -L ~/lib -o Tests.exe Tests.cpp ConwayBromageLib.cpp -lsdsl -ldivsufsort -ldivsufsort64
 */

//Appears to other tests, need to be verify in first
const test critical[] = {
        CASE("merge : "){
            bit_vector ba = {0, 1, 1, 1, 0, 1};
            bit_vector bb = {1, 0, 0, 0};
            sd_vector<>a(ba);
            sd_vector<>b(bb);
            sd_vector<> res = merge(a, b, 4, 1);
            uint64_t finalSize(0);
            if(a.size() > b.size()){
                finalSize = a.size();
            }else{
                finalSize = b.size();
            }
            EXPECT(res.size() == finalSize);
            for(int i = 0 ; i < finalSize ; i++){
                if(a[i] || b[i]){
                    EXPECT(res[i] == 1);
                }
            }
        },
        //Tests of decode :
        CASE("decode case A = 0 ; C = 1 ; G = 2 ; T = 3 : original "){ //decode the 1000000 firsts elements of the generated sequence res (with fromFileToSdVector)
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(uint64_t i = 0 ; i < 1000000 ; i++){
                EXPECT(i == encode(decode(i, len), len));   //expect to find i with a encode-decode combination
                // cout << decode(i, len) << " : " << i << endl;

            }
        },
        CASE("decode case A = 0 ; C = 1 ; T = 2 ; G = 3 (ecoli file) : original for ecoli"){ //decode the 1000000 firsts elements of the generated sequence res (with fromFileToSdVector)
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(uint64_t i = 0 ; i < 1000000 ; i++){
                EXPECT(i == encodeEcoli(decodeEcoli(i, len), len));   //expect to find i with a encode-decode combination
            }
        },
        //Tests of previous :
        CASE("previous : perfect use expected for the 1000000 firsts : "){  //Test with a perfect use of previous + a perfect generated sequence
            vector<vector<uint64_t>> prevOfLittlePrev(16);
            prevOfLittlePrev[1] = {4, 8};
            prevOfLittlePrev[3] = {4, 8};
            prevOfLittlePrev[4] = {1, 5, 13};
            prevOfLittlePrev[5] = {1, 5, 13};
            prevOfLittlePrev[13] = {3, 15};
            prevOfLittlePrev[15] = {3, 15};
            uint64_t currentKMerLen = log(littleTestPrev.size()) / log(ALPHABET);
            for(int i = 0 ; i < littleTestPrev.size() ; i++){ //for the firsts 1000000
                vector<uint64_t> prev = previous(i, littleTestPrev);
                EXPECT(prevOfLittlePrev[i] == prev);
            }
        },
        CASE("previous : unexpected use : a non existant k-mer : "){    //Test with an example of bad use : try to find a non existing k-mer
            vector<uint64_t> prev = previous(-404, littleTestPrev);    //a compressed k-mer can't have a negative number
            EXPECT(prev.size() == 0);   //In this case, they is no possible previous
        },
};

const test lessCritical[] = {
        CASE("reverse complement : original"){
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(int i = 0 ; i < 1000000 ; i++){
                EXPECT(i == reverseComplement(decode(reverseComplement(decode(i, len), totalLen), len), totalLen));
            }
        },
        CASE("reverse complement : faster lexicographical order"){
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(int i = 0 ; i < 1000000 ; i++){
                EXPECT(i == reverseComplementLexico(reverseComplementLexico(i, len), len));
            }
        },
        CASE("reverse complement : fastest for ASCII version"){
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(int i = 0 ; i < 1000000 ; i++){
                EXPECT(i == reverseComplementGATBLibEcoli(reverseComplementGATBLibEcoli(i, len), len));
            }
        },
};

const test atTheEnd[] = {
        //Tests of fromFileToSdVector :
        CASE("fromFileToSdVector : reading from a file : original "){ //verify os fromFileToSdVector can read a file correctly
            sd_vector<>res = fromFileToSdVector("../sorted_kmers.txt");
            EXPECT_NOT((res.size() == 1 && res[0] == 0));   //if res.size() == 1 and res[0] == 0, the file read has failed
        },
        CASE("fromFileToSdVector : reading from a file : original with reverse "){ //verify os fromFileToSdVector can read a file correctly
            sd_vector<>res = fromFileToSdVectorWithReverse("../sorted_kmers.txt");
            EXPECT_NOT((res.size() == 1 && res[0] == 0));   //if res.size() == 1 and res[0] == 0, the file read has failed
        },
        CASE("fromFileToSdVector : reading from a file : Fastest ASCII reading version"){
            sd_vector<>res = fromFileToSdVectorEcoli("./ecoli_count.txt");
            EXPECT_NOT((res.size() == 1 && res[0] == 0));
        },
        CASE("fromFileToSdVector : reading from a file : Fastest ASCII reading version reverse complement included"){
            sd_vector<>res = fromFileToSdVectorWithReverseEcoli("./ecoli_count.txt");
            EXPECT_NOT((res.size() == 1 && res[0] == 0));
        },
};

int main(){
    int failures(0);
    cout << "*** First tests : merge, decode and previous ***" << endl;
    if(failures = run(critical)){  //run all the test cases
        return failures;
    }
    cout << "*** First tests passed ! ***" << endl;
    cout << "*** Second tests : reverse complement ***" << endl;
    if (failures = run(lessCritical)){
        return failures;
    }
    cout << "*** Second tests passed ! ***" << endl;
    cout << "*** Third tests : fromFileToSdVector ***" << endl;
    if (failures = run(atTheEnd)){
        return failures;
    }
    cout << "*** Third tests passed ! ***" << endl;
    return cout << "ALL TESTS PASSED !!!\n", EXIT_SUCCESS;
}

//Ancient tests next and decode : need to implement next with lest library

/*
void test_next(){
    //lambda to print a vector<string>
    auto print_vec = [&](vector<string> successors){
        cout << "{ ";
        for(int i = 0; i < successors.size(); i++)
            cout << successors[i] << " ";
        cout << "}" << endl;
    };
    printf("Test of the function next : \n");
    //test for k=4 => size = 256 (4^4)
    sd_vector<> res = fromFileToSdVector("./sorted_kmers.txt");
    int k = (int)(log(res.size())/log(4));
    int sd_vector_size = 256;
    for(int i = 0; i < sd_vector_size; i++){
        vector<string> successors = next(decode(i, k),res);
        cout << "next of " << decode(i, k) << " -> "; 
        print_vec(successors);
    }
    //test of nexts : what happens if the user give a non existing kmer in the sd_vector
    // for example, what happens if we give "AAAAA" or "AAA" in a sd_vector of 4-mer
    vector<string> successors;
    successors = next("AATTA",res);
    cout << "next of AATTA -> "; 
    print_vec(successors); //has to return {} because AATTA has a size superior to 4

    successors = next("AAT",res);
    cout << "next of AAT -> "; //has to return {} because AATTA has a size inferior to 4
    print_vec(successors);
    
}*/


/*void test_decode(){
    printf("Test of the function decode : \n");
    //test for k=4 => size = 256 (4^4)
    int k = 4;
    int sd_vector_size = 256;
    int cpt_success = 0;
    for(int i = 0; i < sd_vector_size; i++){
        string dec = decode(i, k);
        if(encode(dec, k) == i){
            cpt_success++;
            printf("%d <-> %s SUCCESS \n",i, dec.c_str());
        } else {
            printf("%d <-> %s FAIL    \n",i, dec.c_str());
        }
    }
    printf("Number of successful test : %d/%d\n", cpt_success, sd_vector_size);
}*/

