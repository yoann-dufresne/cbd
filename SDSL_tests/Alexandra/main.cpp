#include <sdsl/bit_vectors.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>
#include <fstream>
#include <ctime>
#include <cstdlib>

/*
 * Compilation line (after installing sdsl libraries):
 * g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib main.cpp -o myOut -lsdsl -ldivsufsort -ldivsufsort64
 */

using namespace sdsl;
using namespace std;

//build of a bit vector in function of its size
bit_vector bitVectorCreator(int size){
    bit_vector myCustomB(size);
    int myRandomNumber = 0;
    //Use of random number to simulate sparsing
    for(int i = 0 ; i < size ; i++){
        myRandomNumber = rand() % 10;
        if(myRandomNumber < 1){
            myCustomB[i] = 1;
        }else{
            myCustomB[i] = 0;
        }
    }
    return myCustomB;
}
/*try to create concatenation operator for sd_vector
 * We can't add an element at the end of a sd_vector directly, ex : sd_b[3] = 0
 * To add, we use a temporary bit_vector, it can be a problem because of size
 * It needs improvements
*/
 sd_vector<> operator+(sd_vector<> const& a, sd_vector<> const& b){
    int myCustomLen = a.size() + b.size();
    bit_vector myCustomBitVec(myCustomLen,0);   //temporary bit_vector
    int j(0);
    for(int i = 0 ; i < myCustomLen ; i++){
        if(i >= a.size()){
            myCustomBitVec[i] = b[j];
            j++;
        }else{
            myCustomBitVec[i] = a[i];
        }
    }
    sd_vector<>myCustomSdVec(myCustomBitVec);   //final sd_vector we will return
    return myCustomSdVec;
}

int main() {
    srand(time(0));
    //Creation of a bit vector
    bit_vector b = bitVectorCreator(10000);
    int count = 0;
    for(int i = 0 ; i < b.size() ; i++){    //counting zeros
        if(b[i] == 0){
            count++;
        }
    }
    //We will make size comparison of bit_vector and sd_vector
    cout << "Number of 0 : " << count << " for a size of " << b.size() << " elements" << endl;
    cout << "Size of b (bytes) : " << size_in_bytes(b) << " B" << endl;
    cout << "Size of b (MB) : " << size_in_mega_bytes(b) << " MB" << endl;
    sd_vector<>sd_b(b); //creation of an sd_vector, based on the bit_vector b
    cout << "Size of sd_b (bytes) : " << size_in_bytes(sd_b) << " B" << endl;
    cout << "Size of sd_b (MB) : " << size_in_mega_bytes(sd_b) << " MB" << endl;
    // test of rank on a sd_vector
    sd_vector<>::rank_1_type sd_b_rank(&sd_b); //rank for 1, can work for 0 or others
    //Print results of rank for each i
    /*for(int i = 0 ; i < sd_b.size() ; i++){
        cout << i << " : " << sd_b_rank(i) << endl;
    }*/
    // test of select on a sd_vector
    size_t len = sd_vector<>::rank_1_type (&sd_b)(sd_b.size()); //count the number of ones in the sd_vector
    cout << "There are " << len << " ones" << endl;
    sd_vector<>::select_1_type  sd_b_select(&sd_b); //select for 1, can work for 0 or others
    //print results of select for each i
    /*for(size_t i = 1 ; i <= len ; i++){
        cout << "One " << i << " : " << sd_b_select(i) << endl;
    }*/
    // try to build a concatenator operation : tests
    bit_vector a = {0, 1, 0, 0, 1}; //smalls bit_vector for tests
    bit_vector c = {1, 1, 1, 0};
    sd_vector<>sd_a(a);
    sd_vector<>sd_c(c);
    sd_vector<>conc = sd_a + sd_c;  //concatenation operator
    for(int i = 0 ; i < conc.size() ; i++){
        cout << conc[i] << " " << endl;
    }
}

