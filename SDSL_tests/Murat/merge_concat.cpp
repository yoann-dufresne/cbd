#include <iostream>
#include <sdsl/sd_vector.hpp>   //in order to manipulate sd_vector

/* Author : Murat Okutucu */
using namespace std;
using namespace sdsl;

/*
 * Return a tuple of two elements :
 * - the first element is the number of one in the resulting sd_vector
 * - the second element is the sd_vector which is the concatenation of a and b
 * @param a sd_vector to whom we want to concatenate
 * @param b sd_vector to concatenate
 * @param nb_of_1_in_A the number of 1-bit in the sd_vector a
 * @param nb_of_1_in_A the number of 1-bit in the sd_vector b
 */
std::tuple<int, sd_vector<>> concatenate (const sd_vector<> &a, const sd_vector<> &b,
                                            int nb_of_1_in_A, int nb_of_1_in_B){
    int len = a.size()+b.size();
    int nb_of_1 = nb_of_1_in_A + nb_of_1_in_B;
    sd_vector_builder s(len, nb_of_1);
    //filling builder
    int i = 0;
    for(; i < a.size(); i++)
        if(a[i]) //==1
            s.set(i); //affect 1 to the indice in argument
    for(; i < len; i++)
        if(b[i-a.size()])
            s.set(i);
    //creation of the sd_vector
    sd_vector<> concat(s);
    return std::make_tuple(nb_of_1, concat);
}

/*
 * Return a tuple of two elements :
 * - the first element is the number of one in the resulting sd_vector
 * - the second element is the sd_vector which is the merge of a and b
 * @param a an sd_vector of the same size as b
 * @param b an sd_vector
 * @param nb_of_1_in_A the number of 1-bit in the sd_vector a
 * @param nb_of_1_in_A the number of 1-bit in the sd_vector b
 */
std::tuple<int, sd_vector<>> merge (const sd_vector<> &a, const sd_vector<> &b, int nb_of_1_in_A, int nb_of_1_in_B){
    int len = a.size();
    int nb_of_1 = nb_of_1_in_A + nb_of_1_in_B;
    int nb_of_1_in_merged = 0;
    sd_vector_builder s(len, nb_of_1);
    //filling builder
    for(int i = 0; i < len; i++) {
        if (a[i] == 1 && b[i] == 1)
            s.set(i);
        if (a[i] == 1 || b[i] == 1) {
            s.set(i);
            nb_of_1_in_merged++;
        }
    }
    //creation of the sd_vector
    sd_vector<> merge(s);
    return std::make_tuple(nb_of_1_in_merged, merge);
}
/*
 * Create a bit_vector of the desired size.
 * @param len the size
 * @return the created bit_vector
 */
bit_vector create_bit_vector(size_t len){
    bit_vector res(len);
    for(int i = 0; i < len; i++){
        res[i] = (rand()%2 < 1)? 1:0;
    }
    return res;
}
/*
 * Print an sd_vector
 * @param sdv the sd_vector
 */
void print_sd_vector(const sd_vector<> &sdv){
    for(int i = 0; i < sdv.size(); i++)
        cout << sdv[i];
    cout << endl;
}

//g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib merge_concat.cpp -o mc -lsdsl -ldivsufsort -ldivsufsort64
int main() {
    srand(time(0));
    //generation of the sd_vectors that we will manipulate in the tests
    int size = 5;
    sd_vector<> A = create_bit_vector(size); //create a sd_vector of a specific size
    sd_vector<> B = create_bit_vector(size);
    int number_of_1_in_A = sd_vector<>::rank_1_type(&A)(A.size());
    int number_of_1_in_B = sd_vector<>::rank_1_type(&B)(B.size());

    //print of A and B
    cout << "A = ";
    print_sd_vector(A);
    cout << "B = ";
    print_sd_vector(B);

    //test of the function "concatenate"
    int nb_of_1_after_concatenation;
    sd_vector<> concat_v;
    tie(nb_of_1_after_concatenation, concat_v) = concatenate(A,B,number_of_1_in_A,number_of_1_in_B);
    cout << "After concatenation, the resulting sd_vector has " << nb_of_1_after_concatenation << " 1-bit in it." << endl;
    //print concat_v
    print_sd_vector(concat_v);

    //test of the function "merge"
    int nb_of_1_after_merging;
    sd_vector<> merge_v;
    tie(nb_of_1_after_merging, merge_v) = merge(A,B,number_of_1_in_A,number_of_1_in_B);
    cout << "After merging, the resulting sd_vector has " << nb_of_1_after_merging << " 1-bit in it." << endl;
    //print merge
    print_sd_vector(merge_v);
    return 0;
}


