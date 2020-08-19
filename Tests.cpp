#include <cstdlib>
#include <string>
#include <bitset>
#include "ConwayBromageLib.h"
#include <lest/lest_basic.hpp>      
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>

using namespace std;
using namespace sdsl;
using namespace lest;

//static const uint64_t  totalLen = 4611686018427387904;
//static const sd_vector<>littleTestPrev = bit_vector{0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1};
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
        //Test of successors on a small sd_vectors
        CASE("successors with ACGT encoding : "){
            cout << "\t--> successors with ACGT encoding" << endl;
            sd_vector<> sdv = bit_vector{0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0};
            KmerManipulatorACGT km(2);
            ConwayBromage cb(sdv, &km);
            
            vector<uint8_t> TrueNext(4);
            TrueNext[0]  = 86;  //01010110
            TrueNext[1]  = 205; //11001101
            TrueNext[2]  = 179; //10110011
            TrueNext[3]  = 106; //01101010
            bool result = true;
            for(int i(0); i < sdv.size()/4; i++){
                uint8_t succ = cb.successors(i);
                //check if they are the same
                if(TrueNext[i] != succ){
                    result = false;
                }

            }
            EXPECT(result);
        },
        CASE("successors with ACTG encoding"){
            cout << "\t--> successors with ACTG encoding" << endl;
            sd_vector<> sdv = bit_vector{0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0};
            KmerManipulatorACTG km(2);
            ConwayBromage cb(sdv, &km);
            
            vector<uint8_t> TrueNext(4);
            TrueNext[0]  = 101;  //01100101
            TrueNext[1]  = 220;  //11011100
            TrueNext[2]  = 166;  //10100110
            TrueNext[3]  = 59;   //00111011
            bool result = true;
            for(int i(0); i < sdv.size()/4; i++){
                uint8_t succ = cb.successors(i);
                //check if they are the same
                if(TrueNext[i] != succ){
                    result = false;
                }

            }
            EXPECT(result);
        },
        CASE("contains : 100 4-mers with ACGT encoding : "){
            cout << "\t--> contains : 100 4-mers with ACGT encoding" << endl;
            ifstream f("./sortACGT.txt", ios::in);
            KmerManipulatorACGT km(4);
            ConwayBromage cb(f, &km);
            f.close();
            ifstream file("./sortACGT.txt", ios::in);
            string line;
            string unCompressedKMer;
            string reverseComplement;
            bool isHere;
            KmerManipulatorACGT kmSize3(3);
            for(int i = 0 ; i < 64 ; i++){
                isHere = false;
                unCompressedKMer = kmSize3.decode(i);
                reverseComplement = kmSize3.decode(kmSize3.reverseComplement(i));
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
                EXPECT(cb.contains(i) == isHere);
                file.clear();
                file.seekg(0, ios::beg);
            }
            file.close();
        },
        CASE("contains : 100 4-mers with ACTG encoding : "){
            cout << "\t--> contains : 100 4-mers with ACTG encoding" << endl;
            ifstream f("./sortACTG.txt", ios::in);
            KmerManipulatorACTG km(4);
            ConwayBromage cb(f, &km);
            f.close();
            ifstream file("./sortACTG.txt", ios::in);
            string line;
            string unCompressedKMer;
            string reverseComplement;
            bool isHere;
            KmerManipulatorACGT kmSize3(3);
            for(int i = 0 ; i < 64 ; i++){
                isHere = false;
                unCompressedKMer = kmSize3.decode(i);
                reverseComplement = kmSize3.decode(kmSize3.reverseComplement(i));
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
                EXPECT(cb.contains(i) == isHere);
                file.clear();
                file.seekg(0, ios::beg);
            }
            file.close();
        },
        CASE("successors with the successors counter : "){  //Modify to fit with the OOP, abandonment of successorCounter
            cout << "\t--> successors : 100 4-mers (no differences for encoding)" << endl;
            ifstream f("./sortACGT.txt", ios::in);
            string call[4]{"A", "C", "G", "T"}; //build comrades manually
            KmerManipulatorACGT k(4);
            KmerManipulatorACGT k3(3);
            ConwayBromage cb(f, &k);
            for(int i = 0 ; i < 64 ; i++){
                bitset<8> bitForm((unsigned)cb.successors(i));  //uint8_t of successors, bit version
                string uncompressedNex(k3.decode(i).erase(0,1));
                string uncompressedPre(k3.decode(i).erase(k3.decode(i).size()-1,1));
                vector<uint64_t> potSucc;
                for(int j = 0 ; j < 8 ; j++){
                    if(bitForm[7-j] == 1){
                        if(j < 4){
                            potSucc.push_back(k.getCanonical(k3.encode(uncompressedNex + call[j])));
                        }else{
                            potSucc.push_back(k.getCanonical(k3.encode(call[j-4] + uncompressedPre)));
                        }
                    }
                }
                for(int j = 0 ; j < potSucc.size() ; j++){
                    EXPECT(cb.contains(potSucc[j]));
                }
            }
        }
};

int main(){
    int failures(0);
    cout << "*** First tests : decode and encode ***" << endl;
    if(failures = run(critical)){  //run all the test cases
        return failures;
    }
    cout << "*** First tests passed ! ***\n" << endl;
    cout << "*** Second tests : reverse complement ***" << endl;
    if (failures = run(lessCritical)){
        return failures;
    }
    cout << "*** Second tests passed ! ***\n" << endl;
    cout << "*** Third tests : getCanonical ***" << endl;
    if (failures = run(lessLessCrit)){
        return failures;
    }
    cout << "*** Third tests passed ! ***\n" << endl;
    cout << "*** Fourth tests : successors and contains ***" << endl;
    if (failures = run(atTheEnd)){
        return failures;
    }
    cout << "*** Fourth tests passed ! ***\n" << endl;
    return cout << "ALL TESTS PASSED !!!\n", EXIT_SUCCESS;
}
