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
static const sd_vector<> ret = fromFileToSdVector("./sorted_kmers.txt", "ACGT");
/*
 * Compilation line : g++ -Wall -Wextra -std=c++11 -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -O3 -DNDEBUG -I ~/include -L ~/lib -o Tests.exe Tests.cpp ConwayBromageLib.cpp -lsdsl -ldivsufsort -ldivsufsort64
 */

//Appears to other tests, need to be verify in first
const test critical[] = {
        CASE("merge : "){
            cout << "\t--> merge" << endl;
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
            cout << "\t--> original decode" << endl;
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(uint64_t i = 0 ; i < 1000000 ; i++){
                        EXPECT(i == encode(decode(i, len), len));   //expect to find i with a encode-decode combination
            }
        },
        //Tests of previous :
        CASE("previous : perfect use expected without fromFileToSdVector call : "){  //Test with a perfect use of previous + a perfect generated sequence
            cout << "\t--> previous with little sd_vector" << endl;
            vector<vector<uint64_t>> prevOfLittlePrev(4);
            prevOfLittlePrev[0] = { 1, 2 };
            prevOfLittlePrev[1] = { 0, 1, 3 };
            prevOfLittlePrev[3] = { 0, 3 };
            for(int i = 0 ; i < 4 ; i++){
                vector<uint64_t> prev = previous(i, littleTestPrev);
                EXPECT(prevOfLittlePrev[i] == prev);
            }
        },
        CASE("previous : unexpected use : a non existant k-mer : "){    //Test with an example of bad use : try to find a non existing k-mer
            cout << "\t--> previous for a non existant k-mer" << endl;
            vector<uint64_t> prev = previous(-404, littleTestPrev);    //a compressed k-mer can't have a negative number
                    EXPECT(prev.size() == 0);   //In this case, they is no possible previous
        },
};

const test lessCritical[] = {
        CASE("reverse complement : faster lexicographical order"){
            cout << "\t--> faster lexicographical order" << endl;
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(int i = 0 ; i < 1000000 ; i++){
                        EXPECT(i == reverseComplementLexico(reverseComplementLexico(i, len), len));
            }
        },
};

const test lessLessCrit[] = {
        //Tests of fromFileToSdVector :
        CASE("fromFileToSdVector : reading from a file : original with normal use "){ //verify os fromFileToSdVector can read a file correctly
            cout << "\t--> fromFileToSdVector with normal use" << endl;
            sd_vector<> res = fromFileToSdVector("./sorted_kmers.txt", "ACGT");
            EXPECT_NOT((res.size() == 1 && res[0] == 0));   //if res.size() == 1 and res[0] == 0, the file read has failed
        },
        CASE("fromFileToSdVector : unexpected format "){ //verify reaction to an unexpected format
            cout << "\t--> fromFileToSdVector with unexpected format" << endl;
            sd_vector<> res = fromFileToSdVector("./sorted_kmers.txt", "RRRR");
                    EXPECT((res.size() == 1 && res[0] == 0));
        },
        CASE("fromFileToSdVector : ACGT format : last element "){
            cout << "\t--> Last element of ACGT" << endl;
            sd_vector<> ret = fromFileToSdVector("./sorted_kmers.txt", "ACGT");
            int currentKMerLen = log(ret.size()) / log(ALPHABET);
            string lastRet = decode(ret.size()-1, currentKMerLen);
            for(int i = 0 ; i < lastRet.size() ; i++){
                EXPECT(lastRet[i] == 'T');
            }
        },
        CASE("fromFileToSdVector : ACTG format : last element "){
            cout << "\t--> Last element of ACTG format" << endl;
            sd_vector<> ret = fromFileToSdVector("./ecoli_count.txt", "ACTG");
            int currentKMerLen = log(ret.size()) / log(ALPHABET);
            string lastRet = decodeEcoli(ret.size()-1, currentKMerLen);
            for(int i = 0 ; i < lastRet.size() ; i++){
                        EXPECT(lastRet[i] == 'G');
            }
        },
};

const test atTheEnd[] = {
        CASE("previous : perfect use expected with file fromFileToSdVectorCall : "){  //Test with a perfect use of previous + a perfect generated sequence
            cout << "\t--> previous with fromFileToSdVector call" << endl;
            uint64_t currentKMerLen = log(ret.size()) / log(ALPHABET);
            for (int i = 0; i < ret.size() / 4; i++) { //for the firsts 1000000
                string current = decode(i, currentKMerLen-1);
                vector<uint64_t> prev = previous(i, ret);
                for(int j = 0 ; j < prev.size() ; j++){
                    vector<uint64_t> proof;
                    proof.push_back(encode((decode(prev[j], currentKMerLen-1).erase(0, 1)) + "A", currentKMerLen-1));
                    proof.push_back(encode((decode(prev[j], currentKMerLen-1).erase(0, 1)) + "C", currentKMerLen-1));
                    proof.push_back(encode((decode(prev[j], currentKMerLen-1).erase(0, 1)) + "G", currentKMerLen-1));
                    proof.push_back(encode((decode(prev[j], currentKMerLen-1).erase(0, 1)) + "T", currentKMerLen-1));
                    EXPECT(i == proof[0] || i == proof[1] || i == proof[2] || i == proof[3]);
                }
            }
        },
        CASE("isThisKMerHere"){
            cout << "\t--> isThisKMerHere : existant K-mer : example without fromFileToSdVector calling :  " << endl;
            bool val = isThisKMerHere(0, littleTestPrev);
            EXPECT(val == true);
        },
        CASE("isThisKMerHere"){
            cout << "\t--> isThisKMerHere : non-existant K-mer : example without fromFileToSdVector calling :  " << endl;
            bool val = isThisKMerHere(-1, littleTestPrev);
                    EXPECT(val == false);
        }
};

int main(){
    int failures(0);
    cout << "*** First tests : merge, decode and previous without file call ***" << endl;
    if(failures = run(critical)){  //run all the test cases
        return failures;
    }
    cout << "*** First tests passed ! ***\n" << endl;
    cout << "*** Second tests : reverse complement ***" << endl;
    if (failures = run(lessCritical)){
        return failures;
    }
    cout << "*** Second tests passed ! ***\n" << endl;
    cout << "*** Third tests : fromFileToSdVector ***" << endl;
    if (failures = run(lessLessCrit)){
        return failures;
    }
    cout << "*** Third tests passed ! ***\n" << endl;
    cout << "*** Fourth tests : previous/next with fromFileToSdVector call and isThisKMerHere ***" << endl;
    if (failures = run(atTheEnd)){
        return failures;
    }
    cout << "*** Fourth tests passed ! ***\n" << endl;
    return cout << "ALL TESTS PASSED !!!\n", EXIT_SUCCESS;
}
