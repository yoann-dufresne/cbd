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
#include "Tests.h"
#include "ConwayBromageLib.h"

using namespace std;
using namespace sdsl;

void test_decode(){
    printf("Test of the function decode : \n");
    //test for k=4 => size = 256 (4^4)
    int k = 4;
    int sd_vector_size = 256;
    int cpt_success = 0;
    for(int i = 0; i < sd_vector_size; i++){
        string dec = decode(i, k);
        if(encode(dec, k) == i){
            cpt_success++;
            printf("%d <-> %s SUCCESS \n",i, dec.c_str());
        } else {
            printf("%d <-> %s FAIL    \n",i, dec.c_str());
        }
    }
    printf("Number of successful test : %d/%d\n", cpt_success, sd_vector_size);
} 

void test_next(){
    //lambda to print a vector<string>
    auto print_vec = [&](vector<string> successors){
        cout << "{ ";
        for(int i = 0; i < successors.size(); i++)
            cout << successors[i] << " ";
        cout << "}" << endl;
    };
    printf("Test of the function next : \n");
    //test for k=4 => size = 256 (4^4)
    sd_vector<> res = fromFileToSdVector("./sorted_kmers.txt");
    int k = (int)(log(res.size())/log(4));
    int sd_vector_size = 256;
    for(int i = 0; i < sd_vector_size; i++){
        vector<string> successors = next(decode(i, k),res);
        cout << "next of " << decode(i, k) << " -> "; 
        print_vec(successors);
    }
    //test of nexts : what happens if the user give a non existing kmer in the sd_vector
    // for example, what happens if we give "AAAAA" or "AAA" in a sd_vector of 4-mer
    vector<string> successors;
    successors = next("AATTA",res);
    cout << "next of AATTA -> "; 
    print_vec(successors); //has to return {} because AATTA has a size superior to 4

    successors = next("AAT",res);
    cout << "next of AAT -> "; //has to return {} because AATTA has a size inferior to 4
    print_vec(successors);
    
}

//Test for previous function
void test_previous(){
    sd_vector<>res = fromFileToSdVector("../sorted_kmers.txt")[0];
    int currentKMerLen = log(res.size()) / log(ALPHABET);
    for(int i = 0 ; i < res.size() ; i++){
        vector<uint64_t> prevDecode = previous(i, res);
        cout << "previous of " << decode(i, currentKMerLen) << " is : { ";
        for(int j = 0 ; j < prevDecode.size() ; j++){
            cout << decode(prevDecode[j], currentKMerLen) << " ";
        }
        cout << "} " << endl;
    }
}
