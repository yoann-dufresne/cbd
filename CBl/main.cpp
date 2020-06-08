#include <iostream>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include "ConwayBromageLib.h"

using namespace std;
using namespace sdsl;

/*
 * Compilation line : g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib main.cpp ConwayBromageLib.cpp -o CBl -lsdsl -ldivsufsort -ldivsufsort64
 */

int main() {
    sd_vector<>res = fromFileToSdVector("../sorted_kmers.txt");
    sd_vector<>::iterator itRes;
    //test
    /*cout << "Print of the sd_vector : " << endl;
    for(itRes = res.begin() ; itRes != res.end() ; itRes++){
        cout << *itRes << " ";
    }
    cout << endl;*/
    cout << "Size of res : " << size_in_mega_bytes(res) << " MB" << endl;
    cout << "Size of res in bytes : " << size_in_bytes(res) << " B" << endl;

    isThisKMerHere("TGTTGATAAC", res);  //test of isThisKmerHere
    vector<string> prevRes = previous("TTATATAACC", res);   //test of previous
    vector<string>::iterator prevIt;
    cout << "Previous elements which are in the sequence : " << endl;
    for(prevIt = prevRes.begin() ; prevIt != prevRes.end() ; prevIt++){
        cout << *prevIt << endl;
    }
    return 0;
}
