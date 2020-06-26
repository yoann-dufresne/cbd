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
        CASE("encode : ACGT encoding : "){    //Up to KmerManipulatorACGT class
            cout << "\t--> encode ACGT encoding" << endl;
                    KmerManipulatorACGT encoder(4);
                    EXPECT(encoder.encode("GGTA") == 172);
        },
        CASE("encode : ACTG encoding : "){           //Up to KmerManipulatorACTG class
            cout << "\t--> encode ACTG encoding" << endl;
                    KmerManipulatorACTG encoder(4);
                    EXPECT(encoder.encode("GGTA") == 248);
        },
        CASE("decode : ACGT encoding : "){    //Up to KmerManipulatorACGT class
            cout << "\t--> decode ACGT encoding" << endl;
                    KmerManipulatorACGT decoder(4);
                    EXPECT(decoder.decode(172) == "GGTA");
        },
        CASE("decode : ACTG encoding : "){        //Up to KmerManipulatorACTG class
            cout << "\t--> decode ACTG encoding" << endl;
                    KmerManipulatorACTG decoder(4);
                    EXPECT(decoder.decode(248) == "GGTA");
        },
};

const test lessCritical[] = {
        CASE("reverse complement : ACGT encoding : "){  //Up to KmerManipulatorACGT class
            cout << "\t--> ACGT encoding" << endl;
            KmerManipulatorACGT reverser(31);
            for(int i = 0 ; i < 1000000 ; i++){
                        EXPECT(i == reverser.reverseComplement(reverser.reverseComplement(i)));
            }
        },
        CASE("reverse complement : ACTG encoding : "){    //Up to KmerManipulatorACTG class
            cout << "\t--> ACTG encoding" << endl;
            KmerManipulatorACTG reverser(31);
            for(int i = 0 ; i < 1000000 ; i++){
                EXPECT(i == reverser.reverseComplement(reverser.reverseComplement(i)));
            }
        },
};

const test lessLessCrit[] = {
        CASE("getCanonical : canonical k-mer ACGT format : "){                           //Up to KmerManipulatorACGT class
            cout << "\t--> canonical k-mer ACGT format" << endl;
                    KmerManipulatorACGT canonicaler(4);
                    EXPECT(canonicaler.getCanonical(172) == 172);
        },
        CASE("getCanonical : uncanonical k-mer ACGT format : "){                           //Up to KmerManipulatorACGT class
            cout << "\t--> uncanonical k-mer ACGT format" << endl;
                    KmerManipulatorACGT canonicaler(4);
                    EXPECT(canonicaler.getCanonical(138) == 93);
        },
        CASE("getCanonical : canonical k-mer ACTG format : "){                           //Up to KmerManipulationACTG class
            cout << "\t--> canonical k-mer ACTG format" << endl;
                    KmerManipulatorACTG canonicaler(4);
                    EXPECT(canonicaler.getCanonical(8) == 8);
        },
        CASE("getCanonical : uncanonical k-mer ACTG format :"){                           //Up to KmerManipulatorACTG class
            cout << "\t--> uncanonical k-mer ACTG format" << endl;
                    KmerManipulatorACTG canonicaler(4);
                    EXPECT(canonicaler.getCanonical(172) == 144);
        },
};

const test atTheEnd[] = {
        CASE("isThisKMerHere : existant K-mer : example without fromFileToSdVector calling : "){
            cout << "\t--> isThisKMerHere : existant K-mer : example without fromFileToSdVector calling" << endl;
            bool val = isThisKMerHere(0, littleTestPrev, 1);
            EXPECT(val == true);
        },
        CASE("isThisKMerHere : non-existant K-mer : example without fromFileToSdVector calling : "){
            cout << "\t--> isThisKMerHere : non-existant K-mer : example without fromFileToSdVector calling" << endl;
            bool val = isThisKMerHere(-1, littleTestPrev, 1);
            EXPECT(val == false);
        },
        CASE("isThisKMerHere : 100 4-mers ACTG encode : "){
            cout << "\t--> isThisKMerHere : 100 4-mers ACTG encode" << endl;
            sd_vector<>ret = fromFileToSdVectorChooser("./sortACTG.txt", "ACTG");
            ifstream file("./sortACTG.txt", ios::in);
            string line;
            string unCompressedKMer;
            string reverseComplement;
            bool isHere;
            for(int i = 0 ; i < 64 ; i++){
                isHere = false;
                unCompressedKMer = decodeEcoli(i, 3);
                reverseComplement = decodeEcoli(reverseComplementGATBLibEcoli(i, 3), 3);
                while(getline(file,line)){
                    if(line.find(unCompressedKMer) != string::npos){
                        isHere = true;
                        break;
                    }
                    if(line.find(reverseComplement) != string::npos){ //if reverse complement is present
                        isHere = true;
                        break;
                    }
                }
                EXPECT(isThisKMerHere(i, ret, false) == isHere);
                file.clear();
                file.seekg(0, ios::beg);
            }
        },
        CASE("isThisKMerHere : 100 4-mers ACGT encode : "){
            cout << "\t--> isThisKMerHere : 100 4-mers ACGT encode" << endl;
            sd_vector<>ret = fromFileToSdVectorChooser("./sortACGT.txt", "ACGT");
            ifstream file("./sortACGT.txt", ios::in);
            string line;
            string unCompressedKMer;
            string reverseComplement;
            bool isHere;
            for(int i = 0 ; i < 64 ; i++){
                isHere = false;
                unCompressedKMer = decode(i, 3);
                reverseComplement = decode(reverseComplementLexico(i, 3), 3);
                while(getline(file,line)){
                    if(line.find(unCompressedKMer) != string::npos){ //if forward is present
                        isHere = true;
                        break;
                    }
                    if(line.find(reverseComplement) != string::npos){ //if reverse complement is present
                        isHere = true;
                        break;
                    }
                }
                EXPECT(isThisKMerHere(i, ret, true) == isHere);
                file.clear();
                file.seekg(0, ios::beg);
            }
        },
        //Test of successors on a small sd_vectors
        CASE("successors with ACGT encoding : "){
            cout << "\t--> successors with ACGT encoding" << endl;
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
        CASE("successors with ACTG encoding"){
            cout << "\t--> successors with ACTG encoding" << endl;
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
        CASE("successors with the successors counter : "){
            cout << "\t--> successors with the successor counter" << endl;
            sd_vector<>ret = fromFileToSdVectorChooser("./sortACTG.txt", "ACTG");
            for(int i = 0 ; i < 64 ; i++){
                vector<uint64_t> counter = successorCounter(i, ret, 0);
                vector<uint64_t> succ = successors(i, ret, 0);
                sort(counter.begin(), counter.end());
                sort(succ.begin(), succ.end());
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
