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
//static const sd_vector<> ret = fromFileToSdVectorChooser("./sorted_kmers.txt","ACGT");
/*
 * Compilation line : g++ -Wall -Wextra -std=c++11 -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -O3 -DNDEBUG -I ~/include -L ~/lib -o Tests.exe Tests.cpp ConwayBromageLib.cpp -lsdsl -ldivsufsort -ldivsufsort64
 */

//Appears to other tests, need to be verify in first
const test critical[] = {
        //Tests of decode :
        CASE("encode lexicographical order : "){
            cout << "\t--> encode lexicographical order" << endl;
            EXPECT(encode("GGTA", 4) == 172);
        },
        CASE("encode ASCII order "){
            cout << "\t--> encode ASCII order" << endl;
            EXPECT(encodeEcoli("GGTA", 4) == 248);
        },
        CASE("decode lexicographical order : "){
            cout << "\t--> decode lexicographical order" << endl;
            EXPECT(decode(172, 4) == "GGTA");
        },
        CASE("decode lexicographical order : "){
            cout << "\t--> decode ASCII order" << endl;
            EXPECT(decodeEcoli(248, 4) == "GGTA");
        }
};

const test lessCritical[] = {
        CASE("reverse complement : faster lexicographical order"){
            cout << "\t--> faster lexicographical order" << endl;
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(int i = 0 ; i < 1000000 ; i++){
                EXPECT(i == reverseComplementLexico(reverseComplementLexico(i, len), len));
            }
        },
        CASE("reverse complement : faster ASCII order"){
            cout << "\t--> faster ASCII order" << endl;
            uint64_t len = log(totalLen) / log(ALPHABET);
            for(int i = 0 ; i < 1000000 ; i++){
                EXPECT(i == reverseComplementGATBLibEcoli(reverseComplementGATBLibEcoli(i, len), len));
            }
        },
};

const test lessLessCrit[] = {
        CASE("getCanonical"){
            cout << "\t--> canonical k-mer ACGT format : " << endl;
            EXPECT(getCanonical(172, 4, 1) == 172);
        },
        CASE("getCanonical"){
            cout << "\t--> uncanonical k-mer ACGT format : " << endl;
            EXPECT(getCanonical(138, 4, 1) == 93);
        },
        CASE("getCanonical"){
            cout << "\t--> canonical k-mer ACTG format : " << endl;
            EXPECT(getCanonical(8, 4, 1) == 8);
        },
        CASE("getCanonical"){
            cout << "\t--> uncanonical k-mer ACTG format : " << endl;
            EXPECT(getCanonical(172, 4, 0) == 144);
        },
};

const test atTheEnd[] = {
        CASE("isThisKMerHere"){
            cout << "\t--> isThisKMerHere : existant K-mer : example without fromFileToSdVector calling :  " << endl;
            bool val = isThisKMerHere(0, littleTestPrev, 1);
            EXPECT(val == true);
        },
        CASE("isThisKMerHere"){
            cout << "\t--> isThisKMerHere : non-existant K-mer : example without fromFileToSdVector calling :  " << endl;
            bool val = isThisKMerHere(-1, littleTestPrev, 1);
            EXPECT(val == false);
        },
        CASE("isThisKMerHere"){
            cout << "\t--> isThisKMerHere : 100 4-mers ACTG encode : " << endl;
            sd_vector<>ret = fromFileToSdVectorChooser("./sortACTG.txt", "ACTG");
            ifstream file("./sortACTG.txt", ios::in);
            string line;
            string unCompressedKMer;
            bool isHere;
            for(int i = 0 ; i < 64 ; i++){
                isHere = false;
                unCompressedKMer = decodeEcoli(i, 3);
                while(getline(file,line)){
                    if(line.find(unCompressedKMer) != string::npos){
                        isHere = true;
                        break;
                    }
                }
                EXPECT(isThisKMerHere(i, ret, 0) == isHere);
                file.clear();
                file.seekg(0, ios::beg);
            }
        },
        CASE("isThisKMerHere"){
            cout << "\t--> isThisKMerHere : 100 4-mers ACGT encode : " << endl;
            sd_vector<>ret = fromFileToSdVectorChooser("./sortACGT.txt", "ACGT");
            ifstream file("./sortACGT.txt", ios::in);
            string line;
            string unCompressedKMer;
            bool isHere;
            for(int i = 0 ; i < 64 ; i++){
                isHere = false;
                unCompressedKMer = decode(i, 3);
                while(getline(file,line)){
                    if(line.find(unCompressedKMer) != string::npos){
                        isHere = true;
                        break;
                    }
                }
                EXPECT(isThisKMerHere(i, ret, 0) == isHere);
                file.clear();
                file.seekg(0, ios::beg);
            }
        },
        //Test of successors on a small sd_vectors
        CASE("successors"){
            cout << "\t--> successors with ACGT encoding : " << endl;
            sd_vector<> sdv = bit_vector{0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0};
            vector<vector<uint64_t>> TrueNext(16);
            TrueNext[0]  = {1, 3, 2};
            TrueNext[1]  = {0, 1, 3};
            TrueNext[2]  = {0, 2, 3};
            TrueNext[3]  = {1, 2, 0};
            bool result = true;
            for(int i(0); i < sdv.size()/4; i++){
                vector<uint64_t> succ = successors(i, sdv, true);
                //check if same
                for(int j(0); j < succ.size(); j++){
                    if(TrueNext[i].size() == succ.size()){
                        for(int j = 0; j < succ.size(); j++){
                            if(succ[j] != TrueNext[i][j]){
                                result = false;
                            }
                        }
                    } else {
                        result = false;
                    }
                }
            }
            EXPECT(result);
        },
        CASE("successors"){
            cout << "\t--> successors with ACTG encoding : " << endl;
            sd_vector<> sdv = bit_vector{0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0};
            vector<vector<uint64_t>> TrueNext(16);
            TrueNext[0]  = {1, 3, 2};
            TrueNext[1]  = {0, 1, 2};
            TrueNext[2]  = {0, 3, 1};
            TrueNext[3]  = {2, 3, 0};
            bool result = true;
            for(int i(0); i < sdv.size()/4; i++){
                vector<uint64_t> succ = successors(i, sdv, false);//alternativeSuccessor(i, sdv, false);//successors(i, sdv, false);
                //check if same
                for(int j(0); j < succ.size(); j++){
                    if(TrueNext[i].size() == succ.size()){
                        for(int j = 0; j < succ.size(); j++){
                            if(succ[j] != TrueNext[i][j]){
                                result = false;
                            }
                        }
                    } else {
                        result = false;
                    }
                }
            }
            EXPECT(result);
        },
        CASE("successors"){
            cout << "\t--> successors with the successor counter : " << endl;
            sd_vector<>ret = fromFileToSdVectorChooser("./sortACTG.txt", "ACTG");
            for(int i = 0 ; i < 64 ; i++){
                vector<uint64_t> counter = successorCounter(i, ret, 0);
                vector<uint64_t> succ = successors(i, ret, 0);
                sort(counter.begin(), counter.end());
                sort(succ.begin(), succ..end());
                EXPECT(succ == counter);
            }
        },
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
