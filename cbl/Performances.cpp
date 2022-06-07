#include <cstdlib>
#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <random>
#include <chrono>
#include "Functions.hpp"
#include "ConwayBromage.hpp"
#include <lest/lest_basic.hpp>
#include <vector>
#include <bitset>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <cmath>
#include <future>

using namespace std;
using namespace sdsl;
using namespace lest;
using namespace std::chrono;

/*Performances.cpp was use to study execution time and memory-space occupation performances
It is not used in the library (ConwayBromageLib.cpp and COnwayBromageLib.h) and you don't need it to run*/
                              
/**
 * Return the next canonical kmer starting from an index and incrementing it.
 * @param start - beginning index
 * @param km - A kmer manipulator
 * @return an index
 */
uint64_t nextCanonical(uint64_t start, KmerManipulator* km){
    uint64_t index = start;
    while(km->getCanonical(index) != index){ //index is canonical if getCanonical returns the same
        index++;
    }
    return index;
}

/**
 * Return an sd_vector with X percent of 1-bit and an vector<int> of all existing k-mers.
 * It's not always possible to generate the exact X percent. For example, for a size of 11 we can't produce an sd_vector with exactly 10% of 1.
 * This is why the parameter producedX will store the produced percentage. In this example, it will be of 9.1 percent.
 * @param len - size of the sd_vector
 * @param X - the percentage
 * @param producedX - the obtained percentage
 * @return the sd_vector and an int_vector
 */
tuple<sd_vector<>, vector<int>> buildSDVWithXPercentOfPresence(uint64_t len, int PmerSize, int X, KmerManipulator* km){
    if(X == 0){ //if X=0 then it must be no 1-bit in the sd_vector
        sd_vector_builder builder(len, 0);
        vector<int> iv(0);
        sd_vector<> sdv(builder);
        return make_tuple(sdv, iv);
    }
    //construction of the variable which permits to take the left-side kmer of a pmer
    uint64_t takeLeftSideKmer = 0;
    int KmerSize = PmerSize-1;
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
    for(uint64_t pmer = 0; nbOf1Placed < nbOf1 && pmer < len; pmer = nextCanonical(pmer+1, km)){
        //build of the sd_vector
        builder.set(pmer);
        nbOf1Placed++;
        //build of the int_vector
        existingKmers[intVecIndex] = pmer >> 2; //we take the right-side kmer of the pmer
        intVecIndex++;
        existingKmers[intVecIndex] = (pmer & takeLeftSideKmer) >> 2; //we take the left-side kmer of the pmer
        intVecIndex++;
    }
    
    cout << "SD_VECTOR : Wanted percentage is " << X << "%. Obtained percentage is " << (nbOf1Placed*100.0/len)<< "%." << endl;
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
vector<int> generateTestList(int size, int percent, const vector<int> &allExistingKmers, const ConwayBromage &cb){
    vector<int> res(size, 0);
    int sdvSize = cb.size();
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
            //while(find(allExistingKmers.begin(), allExistingKmers.end(), randomNb) == allExistingKmers.end()){ //check if the k-mer is present
            while(cb.contains(randomNb)){
                randomNb = randomKmer(genKmer);
            }    
        }
        res[i] = randomNb;
        i++;
    }
    cout << "Wanted percentage for the list of kmers : "<< percent << "%. ";
    cout << "Obtained percentage : " << (cptOfPresent*100.0/res.size()) << "%." << endl;
    return res;
}
/**
 * Produce the file perf.txt with the percentage of presence in the generated test list. 
 * This file will have k-mers apparition probability in each test associated to time (in microsec) taken to accomplish
 * calls to successors function on all k-mers.
 * @param sizeOfRequestLists - size of the generated test list
 */
void launchPerformanceTests(int sizeOfTestList, int PmerSize, KmerManipulator* km){
    ofstream file;
    int probabilityForSDV = 10;
    //creation of the sd_vector with 10% of 1-bit and retrieve of a vector<int> of these k-mers
    auto result = buildSDVWithXPercentOfPresence(pow(4, PmerSize), PmerSize, probabilityForSDV, km);
    sd_vector<> sdv = get<0>(result);
    vector<int> allExistingKmers = get<1>(result);
    
    //KmerManipulatorACTG km(PmerSize);
    ConwayBromage cb(sdv, km);
    cout << "Pmer size : " << cb.getKmerSize() << endl;
    file.open ("perf.txt");
    //calls to successors on lists with different percent of existing k-mers
    for(int percentOfPresence = 10; percentOfPresence <= 90; percentOfPresence+=10){
        //generation of the list of test with X percent of existing Kmers
        vector<int> testList = generateTestList(sizeOfTestList, percentOfPresence, allExistingKmers, cb);
        cout << "testList's size : " << testList.size();
        int cpt = 0;
        for(int i = 0; i < testList.size(); i++){ 
            if(cb.contains(testList[i])==true)
                cpt++;
        }
        double realPercentage = cpt*100.0/testList.size();
        cout << " | Percentage of isPresent " << realPercentage << endl;
        //time measurement of the calls to successors
        auto b = chrono::high_resolution_clock::now();
        for(int i = 0; i < testList.size(); i++) cb.successors(testList[i]);
        auto e = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration_cast<chrono::milliseconds>(e-b).count(); 
        file << realPercentage;
        file << "\t";
        file << elapsed;
        file << "\tms\n";
    }
    file.close();   
}

/**
 * Build a row of element to test timing performance of the function isPresent depending on percentage of elements
 * user wants to be in the original sd_vector
 * @param ratioIn - represent the percentage of elements which are really present in the sd_vector
 * @param nbOfOnes - the number of elements we want in the test row, typically we need 1000000
 * @param cb - ConwayBromage object which contains the sd_vector and shows the encoding format
 * @return an int_vector which contains a percentage of numbers which are present in the sd_vector and which are not contain in
 */
int_vector<> ratioForIsPresent(int ratioIn, int nbOfOnes, ConwayBromage cb){
    srand (time(NULL));
    uint64_t kMerSize = cb.m_kmerSize - 1;  //Size of the k-mer we will build, if sd_vector element is 11-mers, we will build 10-mers
    uint64_t kMerSizeSeq = cb.size() / 4;
    cout << "kmer size : " << kMerSize << endl;
    cout << "kMerSizeSeq : " << kMerSizeSeq << endl;
    cout << "pMerSize : " << cb.size() << endl;
    uint64_t ratioInNumber = ((double)ratioIn/(double)100) * nbOfOnes;  //number which represent the percentage of elements
    int_vector<> al(nbOfOnes, 0);
    vector<uint64_t> isIn;  //Will contain all the elements which are present in the sd_vector
    vector<uint64_t> isOut; //Elements which are NOT present
    uint64_t myMer;
    for(int i = 0 ; i < kMerSizeSeq ; i++){ //We test all k-mers elements. For example, all the 10-mers which can be contain in 11-mers sd_vector
        if(cb.contains(i)){
            isIn.push_back(i);  //is in
        }else{
            isOut.push_back(i); //is not in
        }
    }
    cout << "ratio in number : " << ratioInNumber << endl;
    int i(0);
    while(ratioInNumber > 0){   //We start to complete with elements which are present in the sd_vector
        myMer = rand() % isIn.size();   //Take a random number in the vector which contains present elements
        al[i] = isIn[myMer];
        i++;
        ratioInNumber--;    //count how many present elements we have
    }
    while(i < nbOfOnes){    //Same with elements which are NOT present
        myMer = rand() % isOut.size();  //We complete the int_vector until we have the number of elements we want, typically 1000000 elements
        al[i] = isOut[myMer];           //OBSERVATION : we don't care about duplication, it does not affect isPresent working
        i++;
    }
    return al;
}
/**
 * Call ratioForIsPresent and do time tests on 11 differents percentages (from 0 to 100)
 * Send time and percentage data in a file for graphs building
 */
void metricForIsPresent(){
    ifstream f("./plageTest.txt", ios::in); //File which contains elements to build tests : 1000000 11-mers
    ofstream result("./perfIsPresent.txt"); //File where we will send time and percentage data
    if(result){
        KmerManipulatorACGT k(11);  //Creation of objects
        ConwayBromage cb(f, &k);
        for(int j = 0 ; j <= 100 ; j = j+10){   //We will test isPresent for percentage from 0 to 100
            cout << "cas j = " << j << endl;
            int_vector<> ratio = ratioForIsPresent(j, 1000000, cb); // row test build
            high_resolution_clock::time_point beg = high_resolution_clock::now();   //time measurement starting
            for(int i = 0 ; i < ratio.size() ; i++){
                cb.contains(ratio[i]);
            }
            high_resolution_clock::time_point en = high_resolution_clock::now();    //time measurement ending
            auto durationmicro = duration_cast<microseconds>( en - beg ).count();   //Calculation of time in microseconds
            //auto durationsec = duration_cast<seconds>( en - beg ).count();        //Same in seconds
            result << j;    //Send elements to perfIsPresent.txt, format : presentage   time needed
            result << "\t";
            result << durationmicro << endl;
        }
    }else{
        cout << " Fail while opening file" << endl;
    }
}

int main(){
    int PmerSize = 11;
    KmerManipulatorACTG km(PmerSize);
    launchPerformanceTests(10000000, PmerSize, &km); //launch the test for the function successors and save the results in perf.txt
    string command = "/usr/local/bin/python3 chart.py";
    system(command.c_str()); //plots the results of perf.txt

    metricForIsPresent();   //call of time tests for isPresent
    string command = "/usr/local/bin/python3 graphsForIsPresent.py";    //call of Python script to build graphs
    system(command.c_str());
    return 0;
}
