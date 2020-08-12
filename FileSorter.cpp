#include <iostream>
#include <istream>
#include <vector>
#include <string>
#include <chrono>
#include "ConwayBromageLib.h"

using namespace std;

/**
* Returns a sorted vector with the content of two vectors. 
* @param v1 - first vector
* @param v2 - second vector 
*/
vector<uint64_t> merge(const vector<uint64_t> &v1, const vector<uint64_t> &v2){
    vector<uint64_t> res;
    uint64_t i1 = 0, i2 = 0;
    //uint64_t v1_len = v1.size(), v2_len = v2.size();
    uint64_t v1MaxIndice = v1.size()-1, v2MaxIndice = v2.size()-1;
    //while(i1 < v1_len && i2 < v2_len){
    while(true){
        uint64_t a = v1[i1];
        uint64_t b = v2[i2];
        if(a < b){
            res.push_back(a);
            i1++;
            if(i1 > v1MaxIndice) break;
        } else if(b < a){
            res.push_back(b);
            i2++;
            if(i2 > v2MaxIndice) break;
        } else { //both equal
            res.push_back(a);
            res.push_back(b);
            i1++;
            i2++;
            if(i1 > v1MaxIndice || i2 > v2MaxIndice) break;
        }
    }
    //here i1=v1.size() or i2=v2.size()
    if(i1 == v1.size()){
        while(i2 < v2.size()){
            res.push_back(v2[i2]);
            i2++;
        }
    } else {
        while(i1 < v1.size()){
            res.push_back(v1[i1]);
            i1++;
        }
    }
    return res;
}

/**
* Returns a sorted vector with the contents of many vectors
* @param groups - a vector of vector
*/
vector<uint64_t> mergeGroups(vector<vector<uint64_t>> groups){
    while(groups.size() != 1){
        vector<vector<uint64_t>> final;
        int groups_len = groups.size();
        if(groups_len%2 == 0){
            #pragma omp parallel for
            for(int i = 0; i < groups_len; i+=2){
                final.push_back(merge(groups[i], groups[i+1]));
            }
        } else { //for impair nb of groups, the last one will not be merged during this tour
            #pragma omp parallel for
            for(int i = 0; i < groups_len-1; i+=2){
                final.push_back(merge(groups[i], groups[i+1]));
            }
            final.push_back(groups[groups_len-1]);   
        }
        groups = final;
    }
    cout << "vector max size : " << groups[0].max_size() << endl;
    cout << "final vector size : " << groups[0].size() << endl;
    return groups[0];
}

/**
* Write in a file the k-mers of a vector.
* @param v - a vector
* @param km - a kmer manipulator
* @param path - where to write
*/
void writeInFile(const vector<uint64_t> &v, KmerManipulator* km, string path){
    ofstream fo(path);
    for(uint64_t i = 0; i < v.size(); i++){
        fo << km->decode(v[i]) << "\t1" << endl;
    }
    fo.close();
}

/**
* Sort a file and write in another file
* @param fileToSort - path of the file to sort
* @param outputFile - sorted output
* @param km - a kmer manipulator for the k-mers in the file to sort
*/
void sortFile(string fileToSort, string outputFile, KmerManipulator *km){
    auto X = chrono::high_resolution_clock::now();

    //open file
    ifstream input(fileToSort, ios::in);

    cout << "Creation of the sorted groups. " << endl;;
    //creation of the sorted groups
    auto a = chrono::high_resolution_clock::now(); 
    vector<uint64_t> endLines;
    vector<vector<uint64_t>> groups;
    cout << "vector max size : " << groups.max_size() << endl;
    
    string kmer_str;
    getline(input, kmer_str);
    uint64_t kmer_previous = km->encode(kmer_str);

    uint64_t groupIndex = 0, lineNumber = 2;

    vector<uint64_t> pack;
    groups.push_back(pack);
    groups[groupIndex].push_back(kmer_previous);
    
    while(getline(input, kmer_str)){
        uint64_t kmer = km->encode(kmer_str);

        if(kmer < kmer_previous){
            vector<uint64_t> pack;
            groups.push_back(pack);
            groupIndex++;
            endLines.push_back(lineNumber);
        }
        groups[groupIndex].push_back(kmer);

        lineNumber++;
        kmer_previous = kmer;
    }

    if(endLines.back() < lineNumber) endLines.push_back(lineNumber);

    auto b = chrono::high_resolution_clock::now();
    double elapsed1 = chrono::duration_cast<chrono::seconds>(b-a).count();
    cout << "   >Done in " << elapsed1 << " s." << endl;

    //print info
    cout << "Number of groups : " << groups.size() << endl;
    uint64_t beginning_line = 0;
    for(uint64_t i = 0; i < groups.size(); i++){
        cout << "\t > Group " << (i+1) << " size : " << groups[i].size() << " | in lines [ "<< beginning_line << " ; " << (endLines[i]-1) <<" ]" << endl;
        beginning_line = endLines[i];
    }

    //sort
    cout << "Sorting. " << endl;
    auto c = chrono::high_resolution_clock::now();
    vector<uint64_t> sorted = mergeGroups(groups);
    auto d = chrono::high_resolution_clock::now();
    double elapsed2 = chrono::duration_cast<chrono::seconds>(d-c).count();
    cout << "   >Done in " << elapsed2 << " s." << endl;

    //write in file
    cout << "Writing in file. " << endl;
    auto e = chrono::high_resolution_clock::now();
    writeInFile(sorted, km, outputFile);
    auto f = chrono::high_resolution_clock::now();
    double elapsed3 = chrono::duration_cast<chrono::seconds>(f-e).count();
    cout << "   >Done in " << elapsed3 << " s." << endl;

    //close file
    input.close();

    auto Y = chrono::high_resolution_clock::now();
    double total_time = chrono::duration_cast<chrono::seconds>(Y-X).count();
    cout << "Sorted file produced. TOTAL TIME : " << total_time << " s." << endl;

}

/**
* Print if the file is sorted or not. Generally used after a call to sortFile() in order to check that the output is really sorted.
* @param path - the file
* @param km - a k-mer manipulator to analyze k-mers in the file to sort 
*/
void checkIfSorted(string path, KmerManipulator* km){
    auto X = chrono::high_resolution_clock::now();

    ifstream fi(path, ios::in);
    string line;

    getline(fi, line);
    uint64_t lineNumber = 1, k_previous = km->encode(line);

    bool isSorted = true;
    while(getline(fi, line)){
        lineNumber++;
        uint64_t k = km->encode(line);

        if(k < k_previous){
            isSorted = false;
            cout << path << " : NOT SORTED AT LINE : " << lineNumber << endl;
            cout << "\t" << lineNumber << " : " << line << endl;
        }
        k_previous = k;
    }

    if(isSorted) cout << path << " : FILE IS SORTED" << endl;

    auto Y = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration_cast<chrono::seconds>(Y-X).count();
    cout << "   >Done in " << elapsed << " s." << endl;
}


int main(){
    KmerManipulatorACTG km(31);

    //produce ecoli_sorted.txt
    sortFile("ecoli_count.txt", "ecoli_sorted.txt", &km);

    //check if there is no mistake 
    checkIfSorted("ecoli_sorted.txt", &km);

    return 0;
}
