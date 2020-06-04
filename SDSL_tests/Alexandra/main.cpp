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
/* Update of the function : 04/06
 * New version of the concatenation operator
 * Use of sd_vector_builder to create an sd_vector without a bit_vector
 * No needs to count ones of sd_vector a and b because of rank operations
 * Higher performances than previous version
*/
 sd_vector<> operator+(sd_vector<> const& a, sd_vector<> const& b){
    int myCustomLen = a.size() + b.size();  //length of the future concatenate sd_vector
    int j(0);
    int aLenForOnes = sd_vector<>::rank_1_type (&a)(a.size()); //count the number of ones in a
    int bLenForOnes = sd_vector<>::rank_1_type (&b)(b.size()); //same with b
    /* sd_vector_builder allows us to create a sd_vector without depending on a bit_vector
     * It needs the len of the final sd_vector, and the number of ones in it
     */
    sd_vector_builder svb(myCustomLen, (aLenForOnes + bLenForOnes));
    for(int i = 0 ; i < myCustomLen ; i++){
        if(i >= a.size()){
            if(b[j] == 1){  //If the element is ones, we set the ith case to one
                svb.set(i);
            }
            j++;
        }else{
            if(a[i] == 1){
                svb.set(i);
            }
        }
    }
    cout << "svb study : " << endl;
    sd_vector<>myCustomSdVec(svb);  //Creation of the final sd_vector thanks to the sd_vector_builder svb
    for(int i = 0 ; i < myCustomSdVec.size() ; i++){
        cout << myCustomSdVec[i] << " ";
    }
    cout << endl << "End of svb study" << endl;
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

