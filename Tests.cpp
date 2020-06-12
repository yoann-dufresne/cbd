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
        //Tests of previous :
        CASE("previous : perfect use expected for a manual sd_vector : "){  //Test with a perfect use of previous + a perfect generated sequence
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
        //Test of next upon sorted_kmers.txt
        CASE("next : test upon sorted_kmers.txt "){
            sd_vector<> sdv = fromFileToSdVector("./sorted_kmers.txt");
            int currentKMerLen = log(sdv.size()) / log(ALPHABET);
            for(int i = 0 ; i < sdv.size() ; i++){
                //generation of the true results
                string kmer_initial = decode(i, currentKMerLen);
                string sub = kmer_initial.substr(1, currentKMerLen-1);
                string kmer1 = sub + "A";
                uint64_t next_kmer1 = encode(sub + "A", currentKMerLen);
                uint64_t next_kmer2 = encode(sub + "C", currentKMerLen);
                uint64_t next_kmer3 = encode(sub + "G", currentKMerLen);
                uint64_t next_kmer4 = encode(sub + "T", currentKMerLen);
                //tests if exists
                vector<uint64_t> TrueAnswer;
                if(sdv[next_kmer1]) TrueAnswer.push_back(next_kmer1);
                if(sdv[next_kmer2]) TrueAnswer.push_back(next_kmer2);
                if(sdv[next_kmer3]) TrueAnswer.push_back(next_kmer3);
                if(sdv[next_kmer4]) TrueAnswer.push_back(next_kmer4);

                //comparison with expected results
                vector<uint64_t> successors = next(i, sdv);
                bool result = true;
                if(TrueAnswer.size() == successors.size()){
                    for(int j = 0; j < successors.size(); j++){
                        if(successors[j] != TrueAnswer[j]){
                            result = false;
                        }
                    }
                } else {
                    result = false;
                }
                EXPECT(result);
            }
        },

        //Test of next on a small sd_vector
        CASE("next : small sd_vector"){
            sd_vector<> sdv = bit_vector{0,1,0,1,1,1,0,0,1,0,0,0,0,1,0,1};
            vector<vector<uint64_t>> TrueNext(16);
            TrueNext[1]  = {4,   5};
            TrueNext[3]  = {13, 15};
            TrueNext[4]  = {1,   3};
            TrueNext[5]  = {4,   5};
            TrueNext[8]  = {1,   3};
            TrueNext[13] = {4,   5};
            TrueNext[15] = {13, 15};
            bool result = true;
            for(int i(0); i < sdv.size(); i++){
                vector<uint64_t> successors = next(i, sdv);
                //check if same
                for(int j(0); j < successors.size(); j++){
                    if(TrueNext[i].size() == successors.size()){
                        for(int j = 0; j < successors.size(); j++){
                            if(successors[j] != TrueNext[i][j]){
                                result = false;
                            }
                        }
                    } else {
                        result = false;
                    }
                }
            }
            EXPECT(result);
        }
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
        }
};

const test atTheEnd[] = {
        //Tests of fromFileToSdVector :
        CASE("fromFileToSdVector : reading from a file : original "){ //verify os fromFileToSdVector can read a file correctly
            sd_vector<>res = fromFileToSdVector("./sorted_kmers.txt");
            EXPECT_NOT((res.size() == 1 && res[0] == 0));   //if res.size() == 1 and res[0] == 0, the file read has failed
        }
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

