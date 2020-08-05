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
    cout << "K-mers obtained in vector from : " << fastqFilePath << ". List size : " << res.size() << endl;
    return res;
}

/**
* Generate a fastq file containing random k-mers (doesn't depend on encoding).
* @param outputFile - fastq file that will contain k-mers
* @param numberOfKmer - number of kmer that should be present in the output file
* @kmerSize - size of the k-mers
*/
void generateRandomKmers(string outputFile, int numberOfKmer, int kmerSize){
    ofstream fo(outputFile);
    KmerManipulatorACTG km(kmerSize);
    //generator of random kmers
    random_device rd; 
    mt19937 gen(rd());
    uniform_int_distribution<uint64_t> random(0, pow(4, kmerSize)-1);

    //fill the list
    for(int i = 0; i < numberOfKmer; i++){
        fo << ">kmer" << (i+1) << "\n" << km.decode(random(gen)) << endl;
    }

    cout << "Random k-mers query file generated under the name : " << outputFile << endl;
    fo.close();
}    

/**
* Generate a fastq file containing random existing k-mers.
* @param inputFile - counts file (txt)
* @param outputFile - fastq file containing k-mers
* @cb - ConwayBromage built on file_path
* @km - a kmer manipulator for encoding
* @param numberOfKmer - number of kmer that should be present in the output file
*/
void generateRandomExistingKmers(string inputFile, string outputFile, ConwayBromage const& cb, KmerManipulator* kmForPmer, KmerManipulator* kmForKmer, int numberOfKmer){
    ifstream f(inputFile, ios::in);
    vector<uint64_t> res(numberOfKmer);

    //retrieve some (not all) kmers
    string line;
    int numberOfPmer = 0;
    vector<uint64_t> existingKmers;
    int kmerSize = kmForKmer->getSize();
    while(getline(f, line)){
        numberOfPmer++;
        uint64_t pmer = kmForPmer->encode(line);
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

    ofstream fo(outputFile);
    for(int i = 0; i < res.size(); i++){
        fo << ">kmer" << (i+1) << "\n" << kmForKmer->decode(res[i]) << endl;
    }
    fo.close();
    cout << "Random existing k-mers query file generated under the name : " << outputFile << endl;
}    

/**
* Write in a (fastq) file genomic k-mers from a file containing sequences. The produced file will be used for query.
* @param fastqFilePath - fastq file
* @param queryFile - fastq file containing genomic k-mers to query 
* @param KmerSize - size of cut sequences
* @param maximumNumberOfKmer - maximum number of kmer that will be present in the query file
*/
void generateGenomicKmers(string fastqFilePath, string queryFile, int KmerSize, int maximumNumberOfKmer){
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
    cout << "Genomic k-mers query file generated under the name : " << queryFile << endl;
}

//EXAMPLE OF USE
int main(){
    /*---------------------------
      Generation of genomic Kmers
      ---------------------------*/

        //produce the file query_file_human.fastq which contains 10 k-mers (k=30) from human.fastq
    string output1 = "./query/queryGenomicKmers_human.fastq";
    generateGenomicKmers("./human.fastq", output1, 30, 10);

        //to get the genomic k-mers, we must precise the query file and the kmer manipulator (k=30) in which we want the sequences to be cut
    KmerManipulatorACTG km30(30);
    vector<uint64_t> genomicKmers = getKmersFromQueryFile(output1, &km30);

    /*---------------------------
      Generation of random Kmers
      ---------------------------*/

    string output2 = "./query/queryRandomKmers.fastq";
    generateRandomKmers(output2, 10, 30); //1000 30-mers
    vector<uint64_t> randomKmers = getKmersFromQueryFile(output2, &km30);

    /*---------------------------
      Generation of random existing Kmers
      ---------------------------*/

    ifstream f("./ecoli_count.txt", ios::in);
    KmerManipulatorACTG km31(31);
    ConwayBromage cb(f, &km31);
    f.close();
    string output3 = "./query/queryExistingKmers_ecoli.fastq";
    generateRandomExistingKmers("./ecoli_count.txt", output3, cb, &km31, &km30, 10);
    vector<uint64_t> randomExistingKmers = getKmersFromQueryFile(output3, &km30);
    
    /*
    int numberOfKmer = 100000000;
    //  Generation of genomic Kmers
        //produce the file query_file_human.fastq which contains 10 k-mers (k=30) from human.fastq
    string output1 = "./query/queryGenomicKmers_human.fastq";
    generateGenomicKmers("/pasteur/projets/policy02/seqbio/yo/genomes/human/SRR2052337.1/SRR2052337.1.fastq", output1, 30, numberOfKmer);

        //to get the genomic k-mers, we must precise the query file and the kmer manipulator (k=30) in which we want the sequences to be cut
    KmerManipulatorACTG km30(30);
    vector<uint64_t> genomicKmers = getKmersFromQueryFile(output1, &km30);

    //  Generation of random Kmers
    string output2 = "./query/queryRandomKmers.fastq";
    generateRandomKmers(output2, numberOfKmer, 30); //1000 30-mers
    vector<uint64_t> randomKmers = getKmersFromQueryFile(output2, &km30);

    //  Generation of random existing Kmers
    ifstream f("/pasteur/projets/policy02/seqbio/murat/metagenome.txt", ios::in);
    KmerManipulatorACTG km31(31);
    auto a = chrono::high_resolution_clock::now();
    ConwayBromage cb(f, &km31);
    auto b = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration_cast<chrono::milliseconds>(b-a).count(); 
    cout << "CB construit pour metagenome en " << elapsed << "ms" << endl;
    f.close();
    string output3 = "./query/queryExistingKmers_metagenome.fastq";
    generateRandomExistingKmers("/pasteur/projets/policy02/seqbio/murat/metagenome.txt", output3, cb, &km31, &km30, numberOfKmer);
    vector<uint64_t> randomExistingKmers = getKmersFromQueryFile(output3, &km30);
    */
    return 0;
}
