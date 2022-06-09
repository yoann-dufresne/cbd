#include <cstdlib>
#include <string>
#include <bitset>
//#include "ConwayBromageLib.h"
#include "lest.hpp"
//#include <sdsl/sd_vector.hpp>
//#include <sdsl/vectors.hpp>

using namespace std;
//using namespace sdsl;

/*Tests.cpp performs little unit tests for the library
you don't need it to run*/

//static const uint64_t  totalLen = 4611686018427387904;
//static const sd_vector<>littleTestPrev = bit_vector{0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1};
//static const sd_vector<> ret = fromFileToSdVectorChooser("./sorted_kmers.txt","ACGT");
lest::tests & specification()
{
    static lest::tests test;
    return test;
}

int main( int argc, char * argv[] )
{
    return lest::run( specification(), argc, argv /*, std::cout */ );
}