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
    cout << "Print of the sd_vector : " << endl;
    for(itRes = res.begin() ; itRes != res.end() ; itRes++){
        cout << *itRes << " ";
    }
    cout << endl;
    cout << "Size of res : " << size_in_mega_bytes(res) << " MB" << endl;
    cout << "Size of res in bytes : " << size_in_bytes(res) << " B" << endl;
    return 0;
}
