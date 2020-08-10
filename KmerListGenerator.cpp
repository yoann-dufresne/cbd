#include <cstdlib>
#include <string>
#include <iostream>
#include <istream>
#include <random>
#include <chrono>
#include "ConwayBromageLib.h"
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <cmath>

using namespace std;
using namespace sdsl;
using namespace chrono;

/**
* Returns a list of k-mers contained in a query file.
* @param fastqFilePath - query file
* @param KmerSize - size of cut sequences
*/
vector<uint64_t> getKmersFromQueryFile(string fastqFilePath, KmerManipulator *km){
    vector<uint64_t> res;
    int KmerSize = km->getSize();
    ifstream fi(fastqFilePath, ios::in);

    string kmer_str;
    //sequences are in lineNumber = 2*k with k an integer    
    while(getline(fi, kmer_str)){ //we ignore lines containing ">kmer"
        getline(fi, kmer_str);    //takes the kmer
        res.push_back(km->encode(kmer_str));
    }

    fi.close();
    cout << "K-mers obtained in vector from : " << fastqFilePath << endl;
    return res;
}

/**
* Write in a (fastq) file genomic k-mers from a file containing sequences. The produced file will be used for query.
* @param fastqFilePath - fastq file
* @param queryFile - fastq file containing genomic k-mers to query 
* @param KmerSize - size of cut sequences
* @param numberOfKmer - number of kmer that should be present in the query file
*/
void generateGenomicKmers(string fastqFilePath, string queryFile, int KmerSize, int numberOfKmer){
    KmerManipulatorACTG km(KmerSize);
    ifstream fi(fastqFilePath, ios::in);
    ofstream fo(queryFile);

    string sequence;
    uint64_t lineNumber = 1;
    uint64_t cpt = 1;
    //sequences are in lineNumber = 2*k with k an integer    
    while(getline(fi, sequence)){
        if(lineNumber==2){
            uint64_t numberOfKmerInSequence = 1 + (sequence.size()-KmerSize);
            for(uint64_t i = 0; i < numberOfKmerInSequence; i++){
                string kmer = sequence.substr(i, KmerSize);
                fo << ">kmer" << cpt << "\n" << kmer << endl;
                if(cpt == numberOfKmer){
                    fi.close();
                    fo.close();
                    cout << "Genomic k-mers file generated : " << queryFile << endl;
                    return;
                }
                cpt++;
            }
        }
        lineNumber++;
        if(lineNumber==5) lineNumber = 1;
    }

    fi.close();
    fo.close();
    cout << "Genomic k-mers file generated : " << queryFile << endl;
}

/**
* Generate a fastq file containing random k-mers (doesn't depend on encoding).
* @param outputFile - fastq file that will contain k-mers
* @param numberOfKmer - number of kmer that should be present in the output file
* @kmerSize - size of the k-mers
*/
void generateRandomKmers(string outputFile, uint64_t numberOfKmer, int kmerSize){
    ofstream fo(outputFile);
    KmerManipulatorACTG km(kmerSize);
    //generator of random kmers
    random_device rd; 
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> random(0, pow(4, kmerSize)-1);

    //fill the list
    for(uint64_t i = 0; i < numberOfKmer; i++){
        fo << ">kmer" << (i+1) << "\n" << km.decode(random(gen)) << endl;
    }

    cout << "Random k-mers file generated : " << outputFile << endl;
    fo.close();
}    

/**
* Generate a fastq file containing random existing k-mers.
* @param inputFile - counts file (txt)
* @param outputFile - fastq file containing k-mers
* @kmForKmer - a kmer manipulator for encoding
* @param numberOfKmer - number of kmer that should be present in the output file
*/
void generateRandomExistingKmers(string inputFile, string outputFile, KmerManipulator* km, uint64_t numberOfKmer){
    ifstream fi(inputFile, ios::in);
    ofstream fo(outputFile);

    //retrieve some (not all) kmers
    string line;
    vector<uint64_t> existingKmers;
    int kmerSize = km->getSize();

    while(getline(fi, line))
        existingKmers.push_back(km->encode(line));

    //write random existing kmers in file
    srand((unsigned)time(NULL)); 
    uint64_t limit = existingKmers.size()-1;

    for(uint64_t i = 0; i < numberOfKmer; i++)
        fo << ">kmer" << (i+1) << "\n" << km->decode(existingKmers[rand()%limit]) << endl; //pick a random kmer in existing ones

    fi.close();
    fo.close();
    cout << "Random existing k-mers file generated : " << outputFile << endl;
}    

//EXAMPLE OF USE
int main(){
    KmerManipulatorACTG km30(30);
    KmerManipulatorACTG km31(31);

    /*---------------------------
      Generation of genomic Kmers
      ---------------------------*/
    int numberOfKmer = 100;
    generateGenomicKmers("./human.fastq", "./query/queryGenomicKmers_human.fastq", 30, numberOfKmer);
    //vector<uint64_t> genomicKmers = getKmersFromQueryFile("./query/queryGenomicKmers_human.fastq", &km30);

    /*---------------------------
      Generation of random Kmers
      ---------------------------*/
    generateRandomKmers("./query/queryRandomKmers.fastq", numberOfKmer, 30);
    //vector<uint64_t> randomKmers = getKmersFromQueryFile("./query/queryRandomKmers.fastq", &km30);

    /*---------------------------
      Generation of random existing Kmers
      ---------------------------*/
    generateRandomExistingKmers("./ecoli_sorted.txt", "./query/queryExistingKmers_ecoli.fastq", &km30, numberOfKmer);
    
    //check if all are present (just as a precaution, not mandatory)
    ifstream f("./ecoli_sorted.txt", ios::in);
    auto a = chrono::high_resolution_clock::now();
    ConwayBromage cb(f, &km31);
    auto b = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration_cast<chrono::milliseconds>(b-a).count(); 
    cout << "CB built in " << elapsed << " ms" << endl;
    f.close();

    vector<uint64_t> randomExistingKmers = getKmersFromQueryFile("./query/queryExistingKmers_ecoli.fastq", &km30);
    cout << "randomExistingKmers size : " << randomExistingKmers.size() << endl;
    bool allArePresent = true;
    for(int i = 0; i < randomExistingKmers.size(); i++){
        if(!cb.contains(randomExistingKmers[i])){
            allArePresent = false;
            break;
        }
    }

    if(allArePresent) cout << "All are present !" << endl;
    else cout << "Not all present." << endl;

    return 0;
}
