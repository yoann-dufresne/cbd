/*
 * File:   Tests.cpp
 * Author: muratokutucu
 *
 * Created on 29 juin 2020, 11:07
 */

#include <cstdlib>
#include <string>
#include <iostream>
#include <istream>
#include <random>
#include <chrono>
#include "Functions.hpp"
#include "ConwayBromage.hpp"
#include <lest/lest_basic.hpp> 
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <cmath>
#include <future>


using namespace std;
using namespace sdsl;
using namespace lest;
                              
/*
 * Compilation line : g++ -std=c++11 -Dlest_FEATURE_AUTO_REGISTER=1 -Dlest_FEATURE_COLOURISE=1 -O3 -DNDEBUG -I ~/include -L ~/lib -o exec Performances.cpp ConwayBromage.cpp -lsdsl -ldivsufsort -ldivsufsort64
 */

/**
 * Return an sd_vector with X percent of 1-bit and an vector<int> of all existing k-mers.
 * It's not always possible to generate the exact X percent. For example, for a size of 11 we can't produce an sd_vector with exactly 10% of 1.
 * This is why the parameter producedX will store the produced percentage. In this example, it will be of 9.1 percent.
 * @param len - size of the sd_vector
 * @param X - the percentage
 * @param producedX - the obtained percentage
 * @return the sd_vector and an int_vector
 */
tuple<sd_vector<>, vector<int>> buildSDVWithXPercentOfPresence(uint64_t len, int KmerSize, int X, double &producedX){
    if(X == 0){ //if X=0 then it must be no 1-bit in the sd_vector
        sd_vector_builder builder(len, 0);
        vector<int> iv(0);
        sd_vector<> sdv(builder);
        return make_tuple(sdv, iv);
    }
    //construction of the variable which permits to take the left-side kmer of a pmer
    uint64_t takeLeftSideKmer = 0;
    for(int i = 0; i < 2*KmerSize; i++){
        takeLeftSideKmer <<= 1;
        takeLeftSideKmer += 1;
    }
    takeLeftSideKmer <<= 2; //for KmerSize = 9, we must have 11 11 11 11 11 11 11 11 11 00
    cout << "takeLeftSideKmer = " << takeLeftSideKmer << " Must be equal to (4^P)-4. " << endl;

    //variables
    uint64_t nbOf1 = len * X/100.0;
    uint64_t nbOf1Placed = 0;
    uint64_t pas = 100.0/X;
    //builder of the final sd_vector
    sd_vector_builder builder(len, nbOf1);    
    vector<int> existingKmers(nbOf1*2, 0);
    
    uint64_t intVecIndex = 0;
    for(uint64_t pmer = 0; nbOf1Placed < nbOf1 && pmer < len; pmer+=pas){
        //build of the sd_vector
        builder.set(pmer);
        nbOf1Placed++;
        //build of the int_vector
        existingKmers[intVecIndex] = pmer >> 2; //we take the right-side kmer of the pmer
        intVecIndex++;
        existingKmers[intVecIndex] = (pmer & takeLeftSideKmer) >> 2; //we take the left-side kmer of the pmer
        intVecIndex++;
    }
    
    producedX = nbOf1Placed*100.0/len;
    cout << "SD_VECTOR : Wanted percentage is " << X << "%. Obtained percentage is " << producedX << "%." << endl;
    return make_tuple(sd_vector<>(builder), existingKmers);
}

/**
 * Return a vector<int> composed of X percent of existing k-mer.
 * @param size - desired size of the list.
 * @param percent - desired percentage of presence of the k-mers.
 * @param allExistingKmers - a vector<int> composed of all the existing k-mers in the sd_vector.
 * @param sdvSize - the size of the sd_vector.
 * @param producedPercent - the percentage of presence in the produced vector<int>.
 * @return a vector<int> of k-mer with some percentage existing.
 */
vector<int> generateTestList(int size, int percent, const vector<int> &allExistingKmers, int sdvSize, double &producedPercent){
    vector<int> res(size, 0);
    
    //creation of objects in order to generate a random non existing kmer
    random_device rd1; 
    mt19937 genKmer(rd1());
    uniform_int_distribution<uint64_t> randomKmer(0, sdvSize/4-2);
    
    //permits to obtain the desired percentage
    random_device rd2;
    mt19937 genForChoice(rd2());
    uniform_int_distribution<uint64_t> randomForChoice(0, allExistingKmers.size() * 100.0/percent);
    
    int cptOfPresent = 0, i = 0;
    uint64_t randomNb = 0;
    while(i < size){ //we must fill all the final int_vector
        uint64_t choice = randomForChoice(genForChoice); //generate a number between 0 and present.size() * 100.0/percent
        if(choice < allExistingKmers.size()){ //put in the result a number of the int_vector named present with a proba of 'percent' percent
            randomNb = allExistingKmers[choice]; //we take a k-mer and not a pmer
            cptOfPresent++; //we count it in order to obtain the real percentage of presence that we have produced
        } else {//generation of a non-existing k-mer
            randomNb = randomKmer(genKmer);//generate a number between 0 and sdvSize-1 which coresponds to a pmer
            while(find(allExistingKmers.begin(), allExistingKmers.end(), randomNb) == allExistingKmers.end()){ //check if the k-mer is present
                randomNb = randomKmer(genKmer);
            }    
        }
        res[i] = randomNb;
        i++;
    }
    producedPercent = cptOfPresent*100.0/res.size();
    cout << "Wanted percentage for the list(i.e vector<int>) of kmers : "<< percent << "%. ";
    cout << "Obtained percentage : " << producedPercent << "%." << endl;
    return res;
}
/**
 * Produce the file perf.txt with the percentage of presence in the generated test list. 
 * This file will have k-mers apparition probability in each test associated to time (in microsec) taken to accomplish
 * calls to successors function on all k-mers.
 * @param sizeOfRequestLists - size of the generated test list
 */
void launchPerformanceTests(int sizeOfTestList){
    int PmerSize = 10;
    double producedX = 0;
    ofstream file; 
    auto result = buildSDVWithXPercentOfPresence(pow(4, PmerSize), PmerSize-1, 10, producedX);
    sd_vector<> sdv = get<0>(result);
    vector<int> allExistingKmers = get<1>(result);
    
    KmerManipulatorACGT km(PmerSize);
    ConwayBromage cb(sdv, &km);
    file.open ("perf.txt");
    
    double producedPercent = 0;
    for(int percent = 10; percent <= 90; percent+=10){
        //generation of the list of test with X percent of existing Kmers
        vector<int> testList = generateTestList(sizeOfTestList, percent, allExistingKmers, sdv.size(), producedPercent);
        auto b = chrono::high_resolution_clock::now();
        for(int i = 0; i < testList.size(); i++) cb.successors(testList[i]);
        auto e = chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(e-b).count(); 
        file << producedPercent;
        file << "\t";
        file << elapsed;
        file << "\tmicrosecondes\n";
    }
    file.close();   
}

int main(){
    launchPerformanceTests(10000); //size of the test list
    string command = "/usr/local/bin/python3 chart.py";
    system(command.c_str());
    return 0;
}