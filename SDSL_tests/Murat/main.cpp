#include <sdsl/suffix_arrays.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/select_support.hpp>
#include <fstream>
#include <ctime>
#include <cstdlib>

using namespace sdsl;
using namespace std;

/*
 * Author : Murat Okutucu
 */

/*
 * Param len the wanted size of the bit vector
 * Returns a bit_vector with the size len. It contains ones and zeros. Ones have a probability of about 1% to appear.
 */
bit_vector init(size_t len){
    bit_vector res(len);
    for(int i = 0; i < len; i++){
        res[i] = (rand()%100 < 1)? 1:0;
    }
    return res;
}

// How to compile the file ?
// answer -> g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib main.cpp -o main -lsdsl -ldivsufsort -ldivsufsort64
int main() {
    srand(time(0));
    //default code present in the README
    /* csa_wt<> fm_index;
    construct_im(fm_index, "mississippi!", 1);
    std::cout << "'si' occurs " << count(fm_index,"si") << " times.\n";
    store_to_file(fm_index,"fm_index-file.sdsl");
    std::ofstream out("fm_index-file.sdsl.html");
    write_structure<HTML_FORMAT>(fm_index,out);*/

    //study of the size of bit_vector and sd_vector for different length
    for(int i = 1; i < 10000000; i*=10) {
        bit_vector bv = init(i);
        sd_vector<> sdv = bv;
        cout << "i : " << i << " taille de bit_vector : " << size_in_bytes(bv) << endl;
        cout << "i : " << i << " taille de sdv : " << size_in_bytes(sdv) << endl;
    }

    //I have tried to concatenate two sd_vector but unfortunately I did'nt find how to
    bit_vector b1 = init(10);
    sd_vector<> s1 = b1;

    bit_vector b2 = init(10);
    sd_vector<> s2 = b2;

    for(int i = 0; i < s1.size(); i++) cout << s1[i];
    cout << endl;

    //tests on bit_vector and int vector
    bit_vector b = {1,1,0,1,0,0,1};
    bit_vector bvec(10000000, 0); // taille, contenu
    b.size(); b[2];
    int_vector<> v = {3,2,1,0,2,1,3,4,1,1,1,3,2,3};
    v.size();
    v[1]=0;
    cout << size_in_bytes(v) << endl; //size before compression
    util::bit_compress(v); //permits to compress the vector
    cout << v << endl;
    cout << size_in_bytes(v) << endl; //size after compression

    //tests on the rank operator
    bit_vector bv2 = bit_vector(8000, 0);
    for (size_t i=0; i < bv2.size(); i+=100)
        bv2[i] = 1;
    rank_support_v<1> b_rank(&b); // <- pointer to b
    for (size_t i=0; i<=b.size(); i+= b.size()/4)
        cout << "(" << i << ", " << b_rank(i) << ") ";
    cout << endl;

    return 0;
}