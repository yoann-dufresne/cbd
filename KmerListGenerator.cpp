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
* Returns a list of random k-mers (doesn't depend on encoding).
* @param listSize - size of the list
* @kmerSize - size of the k-mers
*/
vector<uint64_t> getRandomKmerList(int listSize, int kmerSize){
    vector<uint64_t> res(listSize);

    //generator of random kmers
    random_device rd; 
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> random(0, pow(4, kmerSize)-1);

    //fill the list
    for(int i = 0; i < res.size(); i++)
        res[i] = random(gen);

    cout << "Random k-mers list size : " << res.size() << endl;
    return res;
}    

/**
* Returns a list of random k-mers among existing ones in a file.
* @param listSize - size of the list
* @km - a kmer manipulator for encoding
* @file_path - path of the file containing p-mers and their counts.
* @cb - ConwayBromage built on file_path
*/
vector<uint64_t> getRandomKmerListAmongExisting(int listSize, KmerManipulator* km, string file_path, ConwayBromage const& cb){
    vector<uint64_t> res(listSize);
    ifstream f(file_path, ios::in);

    //retrieve some (not all) kmers
    string line;
    int numberOfPmer = 0;
    vector<uint64_t> existingKmers;
    int kmerSize = (km->getSize()-1);
    while(getline(f, line)){
        numberOfPmer++;
        uint64_t pmer = km->encode(line);
        uint64_t kmer1 = pmer >> 2;
        if(cb.contains(kmer1))
            existingKmers.push_back(kmer1);
    }
    f.close();

    //build random kmer list
    srand((unsigned)time(NULL)); 
    int limit = existingKmers.size()-1;;
    for(int i = 0; i < res.size(); i++){
        res[i] = existingKmers[rand()%limit]; //pick a random kmer in existing ones
    }

    //check if all are present
    int presenceCount = 0;
    for(int i = 0; i < res.size(); i++){
        if(cb.contains(res[i])){
            presenceCount++;
        }
    }
    if(presenceCount != res.size()){
        cout << "getRandomKmerListAmongExisting exit 1 because k-mers are not all present";
        exit(1);
    }

    cout << "Random existing k-mers list size : " << res.size() << endl;
    return res;
}    

/**
* Returns a list of(genomic) k-mers contained in the file.
* @param fastqFilePath - query file
* @param KmerSize - size of cut sequences
*/
vector<uint64_t> getKmersFromQueryFile(string fastqFilePath, KmerManipulator *km){
    vector<uint64_t> res;
    int KmerSize = km->getSize();
    ifstream fi(fastqFilePath, ios::in);

    string sequence;
    //sequences are in lineNumber = 2*k with k an integer    
    while(getline(fi, sequence)){ //we ignore the lines containing ">kmer"
        getline(fi, sequence);    //take the kmer
        int numberOfKmer = 1 + (sequence.size()-KmerSize);
        for(int i = 0; i < numberOfKmer; i++){
            string kmer = sequence.substr(i, KmerSize);
            res.push_back(km->encode(kmer));
        }
    }

    fi.close();
    cout << "Genomic k-mers obtained from : " << fastqFilePath << ". List size : " << res.size() << endl;
    return res;
}

/**
* Write in a (fastq) file genomic k-mers from a file containing sequences. The produced file will be used for query.
* @param fastqFilePath - fastq file
* @param queryFile - fastq file containing genomic k-mers to query 
* @param KmerSize - size of cut sequences
* @param maximumNumberOfKmer - maximum number of kmer that will be present in the query file
*/
void fastqToKmerFile(string fastqFilePath, string queryFile, int KmerSize, int maximumNumberOfKmer){
    KmerManipulatorACTG km(KmerSize);
    ifstream fi(fastqFilePath, ios::in);
    ofstream fo(queryFile);

    string sequence;
    int lineNumber = 1;
    int cpt = 1;
    //sequences are in lineNumber = 2*k with k an integer    
    while(getline(fi, sequence)){
        if(lineNumber==2){
            int numberOfKmer = 1 + (sequence.size()-KmerSize);
            for(int i = 0; i < numberOfKmer; i++){
                string kmer = sequence.substr(i, KmerSize);
                fo << ">kmer" << cpt << "\n" << kmer << endl;
                if(cpt == maximumNumberOfKmer){
                    fi.close();
                    fo.close();
                    cout << "Query file has been generated under the name : " << queryFile << endl;
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
    cout << "Query file has been generated under the name : " << queryFile << endl;
}

//EXAMPLE OF USE
int main(){
    //produce the file query_file_human.fastq which contains 300 k-mers (k=30) from human.fastq
    fastqToKmerFile("./human.fastq", "query_file_human.fastq", 30, 10);

    //to get the genomic k-mers, we must precise the query file and the kmer manipulator (k=30) in which we want the sequences to be cut
    KmerManipulatorACTG km30(30);
    vector<uint64_t> genomicKmers = getKmersFromQueryFile("./query_file_human.fastq", &km30);

    //random kmers
    vector<uint64_t> randomKmers = getRandomKmerList(1000000, 30); //list size of 1000000

    //random kmers among existing ones
    ifstream f("./ecoli_count.txt", ios::in);
    KmerManipulatorACTG km31(31);
    ConwayBromage cb(f, &km31);
    f.close();
    vector<uint64_t> randomExistingKmers = getRandomKmerListAmongExisting(1000000, &km31, "./ecoli_count.txt", cb);
    
    return 0;
}
