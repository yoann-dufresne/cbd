#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include "ConwayBromageLib.h"
#include <immintrin.h>  //need -mavx2

using namespace std;
using namespace sdsl;

/**
 * @brief need to be used only for creating version of this object with the vector in it so don't use it directly
 * 
 * @param km 
 */
ConwayBromage::ConwayBromage(KmerManipulator* km){
    m_kmerManipulator = km;

    string line("");
    uint64_t one = 1;
    uint64_t sdvSize = one << ((2*m_kmerManipulator->getSize())-1); //Creation of the total length to create the sd_vector_builder
    //sd_vector_builder builder(sdvSize, numberOfKmer);
    m_limit = (sdvSize >> 1);//changed from the original code source, otherwise it block contains and successor for some case

    //initialization of the cache
    int nbOfBitsToShift = 2 * (m_kmerManipulator->getSize()-1);
    m_numberOfBitsToShift = nbOfBitsToShift;

    m_RC[0] = km->reverseComplementOfNucleotide(0);
    m_RC[1] = km->reverseComplementOfNucleotide(1);
    m_RC[2] = km->reverseComplementOfNucleotide(2);
    m_RC[3] = km->reverseComplementOfNucleotide(3);

    m_RC_shifted[0] = m_RC[0] << nbOfBitsToShift;
    m_RC_shifted[1] = m_RC[1] << nbOfBitsToShift;
    m_RC_shifted[2] = m_RC[2] << nbOfBitsToShift;
    m_RC_shifted[3] = m_RC[3] << nbOfBitsToShift;
    

    m_nucleotides_shifted[0] = (uint64_t)0 << nbOfBitsToShift;
    m_nucleotides_shifted[1] = (uint64_t)1 << nbOfBitsToShift;
    m_nucleotides_shifted[2] = (uint64_t)2 << nbOfBitsToShift;
    m_nucleotides_shifted[3] = (uint64_t)3 << nbOfBitsToShift;

    for(int i = 0; i < 4; i++){
        //for next and previous pmers
        char nucleotideOfI = m_kmerManipulator->decodeNucleotide(i);
        if(nucleotideOfI == 'A'){ 
            m_correspondingBitValueForNextPmers[i] = 128;
            m_correspondingBitValueForPrevPmers[i] = 8;  
        } else if(nucleotideOfI == 'C') {     
            m_correspondingBitValueForNextPmers[i] = 64; 
            m_correspondingBitValueForPrevPmers[i] = 4; 
        } else if(nucleotideOfI == 'G') {
            m_correspondingBitValueForNextPmers[i] = 32; 
            m_correspondingBitValueForPrevPmers[i] = 2; 
        } else {//'T'
            m_correspondingBitValueForNextPmers[i] = 16; 
            m_correspondingBitValueForPrevPmers[i] = 1; 
        }
    }
}
/**
 * Returns the kmer manipulator.
 * @return a KmerManipulator object.
 */
KmerManipulator* ConwayBromage::getKmerManipulator(){
    return m_kmerManipulator;
}

/**
 * Returns the size k of the k-mer that are stored.
 * @return an int
 */
int ConwayBromage::getKmerSize(){
    return m_kmerManipulator->getSize();
}

//
//  ConwayBromage object with a sdsl succint bit-vector, first implementation used in this project
//


/**
 * First ConwayBromageSD constructor : Stores memory-efficiently k-mers coming from an istream in a specific format.  
 * @param kmerFlux - An istream of k-mers. For example a file represented by an ifstream.
 * @param km - Object which contains all the information about encoding and decoding in a specific format.
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * int kmerSize = 31;
 * ifstream f("./k-mers.txt", ios::in);
 * KmerManipulatorACGT km(kmerSize);
 * ConwayBromage cb(f, &km);        
 * f.close();
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * @warning The k-mers in the istream must me canonical and have a size k <= 32.
 * @warning The istream must have at each line (if it's a file for example) a unique k-mers.
 * @warning The k-mers in the istream must be sorted either in lexicographical order (A<C<G<T) or in A<C<T<G. It depends on the encoding format.
 */
ConwayBromageSD::ConwayBromageSD(istream& kmerFlux, KmerManipulator* km): ConwayBromage(km){

    string line("");
    uint64_t numberOfKmer = 0;
    while(getline(kmerFlux, line)) //Counts the number of k-mer in the file
        numberOfKmer++;
 
    kmerFlux.clear();
    kmerFlux.seekg(0, ios::beg);    //Return to the beginning of the file
    uint64_t one = 1;
    uint64_t sdvSize = one << ((2*m_kmerManipulator->getSize())-1); //Creation of the total length to create the sd_vector_builder
    //sdsl need a builder to manipulate the bit-vector at creation
    sd_vector_builder builder(sdvSize, numberOfKmer);
    uint64_t previousKmer(0);
    while(getline(kmerFlux, line)){
        uint64_t k = m_kmerManipulator->encode(line);
        //first chack : canonical elements expected
        if(k >  m_kmerManipulator->reverseComplement(k)){
            
            cout << line<<"The file is not completely canonical" << endl;
            exit(1);
        }
        //second check : ascending order expected
        if(k < previousKmer){
            cout << "The file is not sort in the ascending order" << endl;
            exit(1);
        }
        builder.set(k);   
        previousKmer = k;
    }
    
    m_sequence = builder;
}

/**
 * Second ConwayBromageSD constructor : Build the ConwayBromage object based on the sd_vector in parameter and the KmerManipulator.
 * @param sdv - An sd_vector which represents the sequence.
 * @param km - A KmerManipulator.
 */
ConwayBromageSD::ConwayBromageSD(sdsl::sd_vector<> const& sdv, KmerManipulator* km) : ConwayBromage(km){
    m_sequence=sd_vector(sdv); //copy of the sd_vector
}

/**
 * Check if the given (k-1)-mer is present. The (k-1)-mer can either be canonical or not.
 * @param Kmer - An uint64_t representing the (k-1)-mer.
 * @return True if the (k-1)-mer is present. False otherwise. 
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * bool kmer_exists = cb.contains(230);
 *
 * KmerManipulatorACGT km(3);
 * uint64_t intGTT = km.encode("GTT");
 * bool GTT_exists = cb.contains(intGTT);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
bool ConwayBromageSD::contains(uint64_t Kmer) const{
    if(Kmer > m_limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << m_limit << endl;
        return false; 
    }  
    
    uint64_t PmerPrev, RC_PmerNext, RC_PmerPrev;
    uint64_t PmerNext = Kmer << 2; //<-> XA where X is the Kmer
    uint64_t RC_Kmer = m_kmerManipulator->reverseComplement(Kmer) >> 2;
    uint64_t RC_Kmer_ShiftedOf2Bits = RC_Kmer << 2;

    for(int i = 0; i < 4; i++){ 
        //next
        RC_PmerNext = RC_Kmer + m_RC_shifted[i]; //YX where X is the reverse complement of the kmer and Y is the reverse complement of the last nucleotide of the pmer
        if(m_sequence[(PmerNext < RC_PmerNext)?PmerNext:RC_PmerNext]) return true;
        //previous
        PmerPrev = Kmer + m_nucleotides_shifted[i]; //equals to AX then CX then GX then TX where X is the Kmer
        RC_PmerPrev = RC_Kmer_ShiftedOf2Bits + m_RC[i];
        if(m_sequence[(PmerPrev < RC_PmerPrev)?PmerPrev:RC_PmerPrev]) return true;
        
        PmerNext++; //equals to XA then XC then XG then XT
    }
    return false;
}

/**
 * Returns the successors of a (k-1)-mer. The (k-1)-mer can either be canonical or not. A (k-1)-mer has at minimum 0 successors and maximum 8 successors.
 * ### How it works ? 
 * Let's say the function takes as a parameter the integer which represents the (k-1)-mer GTT and we assume we are in ACGT encoding.
 * First, the method will generate the following k-mers : GTTA, GTTC, GTTG, GTTT, AGTT, CGTT, GGTT, TGTT.
 * Then, it will check if the canonical version of these k-mers are present.
 * Let's suppose that only the canonical version of these k-mers are present : GTTA, GTTT, CGTT, TGTT.
 * GTTA and GTTT are, what we call next k-mers so we will store the information in the 4 left bits of the result (uint8_t).
 * CGTT and TGTT are previous k-mers so, this time, the information will be stored in the 4 right bits of the result.
 * Thus, the function will return 1001 0101 which corresponds to 149 in base 10.
 * It means that, in this example, the successors of the (k-1)-mer GTT are TTA, TTT, CGT and TGT.
 *
 * @param Kmer : a (k-1)-mer.
 * @return a uint8_t which carry information about the presence/absence of the 8 potential successors of the (k-1)-mer.
 * Read of the return : 
 * read the uint8_t under bit form : the first four elements are "next" k-mers, the last four are "previous"
 * Always the same reading direction : A then C then G then T regardless of the encoding
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * uint8_t successors = cb.successors(78);
 * KmerManipulatorACGT km(3);
 * uint64_t intGTT = km.encode("GTT");
 * uint8_t successorsOfGTT = cb.successors(intGTT);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * @warning The method doesn't check the presence of the (k-1)-mer (so you have to do it on your own with the method contains for example).
 */
uint8_t ConwayBromageSD::successors(uint64_t Kmer) const{
    if(Kmer > m_limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be equal or inferior to 4^(P-1) i.e " << m_limit << endl;
        return 0; //empty
    }
    uint8_t res = 0;
    //we build the eight possible successors    
    uint64_t PmerPrev, RC_PmerNext, RC_PmerPrev;
    uint64_t PmerNext = Kmer << 2; //<-> XA where X is the Kmer
    uint64_t RC_Kmer = m_kmerManipulator->reverseComplement(Kmer) >> 2;
    uint64_t RC_Kmer_ShiftedOf2Bits = RC_Kmer << 2;

    for(int i = 0; i < 4; i++){ 
        RC_PmerNext = RC_Kmer + m_RC_shifted[i]; //YX where X is the reverse complement of the kmer and Y is the reverse complement of the last nucleotide of the pmer
        if(m_sequence[(PmerNext < RC_PmerNext)?PmerNext:RC_PmerNext]){ //if the pmer is present 
            res ^= m_correspondingBitValueForNextPmers[i];
        }
        PmerPrev = Kmer + m_nucleotides_shifted[i]; //equals to AX then CX then GX then TX where X is the Kmer
        RC_PmerPrev = RC_Kmer_ShiftedOf2Bits + m_RC[i];
        if(m_sequence[(PmerPrev < RC_PmerPrev)?PmerPrev:RC_PmerPrev]){
            res ^= m_correspondingBitValueForPrevPmers[i];
        }
        PmerNext++; //equals to XA then XC then XG then XT
    }
    
    return res;
}

/**
 * Returns the compressed sequence.
 * @return an sd_vector
 */
sdsl::sd_vector<> ConwayBromageSD::getSequence(){
    return m_sequence;
}

int ConwayBromageSD::serialize(std::string filename){
    store_to_file(m_sequence,filename);
    return 0;
}
ConwayBromageSD ConwayBromageSD::deserialize(std::string filename,KmerManipulator* km){  
    sdsl::sd_vector<> tmp;
    load_from_file(tmp,filename);
    ConwayBromageSD ret(tmp,km);
    return ret;
}


//
//  an Intermediate object to contain a vector of bitvector to circumvent the fact that the max of bm is 2^48
//

/**
 * @brief 
 * 
 * @param nb the number of bvector we want in our vector
 */
Intermediate::Intermediate(int nb){
    std::vector<bm::bvector<>> tmp2;
    nbv=nb;
    vbv=tmp2;
    vbv.resize(nb);
    for(int i=0;i<nb;i++){
        bm::bvector<> tmp;
        tmp.set(0,false);
        tmp.set_new_blocks_strat(bm::BM_GAP);
        vbv.push_back(tmp);
    }
    
}
Intermediate::Intermediate(){
    Intermediate(1);
}
Intermediate::Intermediate(std::vector<bm::bvector<>>& tmp2,int a){
    vbv=tmp2;
    nbv=a;
}
void Intermediate::set(uint64_t id){
    uint64_t i=id & mask;
    i=i>>48;
    id=id&mask2;
    vbv.at(i).set(id);
}

int Intermediate::present(uint64_t id)const{
    uint64_t i=id & mask;
    i=i>>48;
    id=id&mask2;
    return vbv.at(i)[id];
}
void Intermediate::optimize(){
    for(int i=0;i<vbv.size();i++){
        vbv[i].optimize();
    }
}
//
//  ConwayBromage object with a bit-magic succint bit-vector, to test the perfomance difference between the 2 implementation
//


/**
 * First ConwayBromage constructor : Stores memory-efficiently k-mers coming from an istream in a specific format.  
 * @param kmerFlux - An istream of k-mers. For example a file represented by an ifstream.
 * @param km - Object which contains all the information about encoding and decoding in a specific format.
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * int kmerSize = 31;
 * ifstream f("./k-mers.txt", ios::in);
 * KmerManipulatorACGT km(kmerSize);
 * ConwayBromage cb(f, &km);        
 * f.close();
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * @warning The k-mers in the istream must me canonical and have a size k <= 32.
 * @warning The istream must have at each line (if it's a file for example) a unique k-mers.
 * @warning The k-mers in the istream must be sorted either in lexicographical order (A<C<G<T) or in A<C<T<G. It depends on the encoding format.
 */
ConwayBromageBM::ConwayBromageBM(istream& kmerFlux, KmerManipulator* km) : ConwayBromage(km){
    Intermediate builder(1);
    m_sequence=builder;
    if(km->getSize()>24){
        uint64_t number=0;
        number=(uint64_t)1<<((km->getSize()-24)*2);// the number of bvector needed to have the equivalent of one bvector with the  range we need
        Intermediate builder(number);
        m_sequence=builder;
    }

    string line("");
    uint64_t numberOfKmer = 0;
    while(getline(kmerFlux, line)) //Counts the number of k-mer in the file
        numberOfKmer++;
 
    kmerFlux.clear();
    kmerFlux.seekg(0, ios::beg);    //Return to the beginning of the file
    uint64_t one = 1;
    uint64_t sdvSize = one << ((2*m_kmerManipulator->getSize())-1); //Creation of the total length to create the sd_vector_builder
    uint64_t previousKmer(0);
    while(getline(kmerFlux, line)){
        uint64_t k = m_kmerManipulator->encode(line);
        //first chack : canonical elements expected
        if(k >  m_kmerManipulator->reverseComplement(k)){
            cout << "The file is not completely canonical" << endl;
            exit(1);
        }
        //second check : ascending order expected
        if(k < previousKmer){
            cout << "The file is not sort in the ascending order" << endl;
            exit(1);
        }
        m_sequence.set(k);   
        previousKmer = k;
    }
    m_sequence.optimize();
    m_limit = (sdvSize >> 1);//changed from the original code source, otherwise it block contains and successor for some case

}

/**
 * Second ConwayBromage constructor : Build the ConwayBromage object based on the bvector in parameter and the KmerManipulator.
 * @param sdv - An sd_vector which represents the sequence.
 * @param km - A KmerManipulator.
 */
ConwayBromageBM::ConwayBromageBM(Intermediate&  sdv, KmerManipulator* km) : ConwayBromage(km){
    m_sequence = sdv; //copy of the sd_vector
}

/**
 * Check if the given (k-1)-mer is present. The (k-1)-mer can either be canonical or not.
 * @param Kmer - An uint64_t representing the (k-1)-mer.
 * @return True if the (k-1)-mer is present. False otherwise. 
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * bool kmer_exists = cb.contains(230);
 *
 * KmerManipulatorACGT km(3);
 * uint64_t intGTT = km.encode("GTT");
 * bool GTT_exists = cb.contains(intGTT);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
bool ConwayBromageBM::contains(uint64_t Kmer) const{
    if(Kmer > m_limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be strictly inferior to 4^(P-1) i.e " << m_limit << endl;
        return false; 
    }  
    
    uint64_t PmerPrev, RC_PmerNext, RC_PmerPrev;
    uint64_t PmerNext = Kmer << 2; //<-> XA where X is the Kmer
    uint64_t RC_Kmer = m_kmerManipulator->reverseComplement(Kmer) >> 2;
    uint64_t RC_Kmer_ShiftedOf2Bits = RC_Kmer << 2;

    for(int i = 0; i < 4; i++){ 
        //next
        RC_PmerNext = RC_Kmer + m_RC_shifted[i]; //YX where X is the reverse complement of the kmer and Y is the reverse complement of the last nucleotide of the pmer
        
        if(m_sequence.present((PmerNext < RC_PmerNext)?PmerNext:RC_PmerNext)){
            return true;
        } 

        //previous
        PmerPrev = Kmer + m_nucleotides_shifted[i]; //equals to AX then CX then GX then TX where X is the Kmer
        RC_PmerPrev = RC_Kmer_ShiftedOf2Bits + m_RC[i];
        if(m_sequence.present((PmerPrev < RC_PmerPrev)?PmerPrev:RC_PmerPrev)){ 
            return true;
        }
        
        PmerNext++; //equals to XA then XC then XG then XT
    }
    return false;
}

/**
 * Returns the successors of a (k-1)-mer. The (k-1)-mer can either be canonical or not. A (k-1)-mer has at minimum 0 successors and maximum 8 successors.
 * ### How it works ? 
 * Let's say the function takes as a parameter the integer which represents the (k-1)-mer GTT and we assume we are in ACGT encoding.
 * First, the method will generate the following k-mers : GTTA, GTTC, GTTG, GTTT, AGTT, CGTT, GGTT, TGTT.
 * Then, it will check if the canonical version of these k-mers are present.
 * Let's suppose that only the canonical version of these k-mers are present : GTTA, GTTT, CGTT, TGTT.
 * GTTA and GTTT are, what we call next k-mers so we will store the information in the 4 left bits of the result (uint8_t).
 * CGTT and TGTT are previous k-mers so, this time, the information will be stored in the 4 right bits of the result.
 * Thus, the function will return 1001 0101 which corresponds to 149 in base 10.
 * It means that, in this example, the successors of the (k-1)-mer GTT are TTA, TTT, CGT and TGT.
 *
 * @param Kmer : a (k-1)-mer.
 * @return a uint8_t which carry information about the presence/absence of the 8 potential successors of the (k-1)-mer.
 * Read of the return : 
 * read the uint8_t under bit form : the first four elements are "next" k-mers, the last four are "previous"
 * Always the same reading direction : A then C then G then T regardless of the encoding
 * ### Example
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.cpp
 * uint8_t successors = cb.successors(78);
 * KmerManipulatorACGT km(3);
 * uint64_t intGTT = km.encode("GTT");
 * uint8_t successorsOfGTT = cb.successors(intGTT);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * @warning The method doesn't check the presence of the (k-1)-mer (so you have to do it on your own with the method contains for example).
 */
uint8_t ConwayBromageBM::successors(uint64_t Kmer) const{
    if(Kmer > m_limit) { //we must have nonCompressedKmer < 4^(P-1)
        cout << "The value of the kmer must be equal or inferior to 4^(P-1) i.e " << m_limit << endl;
        return 0; //empty
    }
    uint8_t res = 0;
    //we build the eight possible successors    
    uint64_t PmerPrev, RC_PmerNext, RC_PmerPrev;
    uint64_t PmerNext = Kmer << 2; //<-> XA where X is the Kmer
    uint64_t RC_Kmer = m_kmerManipulator->reverseComplement(Kmer) >> 2;
    uint64_t RC_Kmer_ShiftedOf2Bits = RC_Kmer << 2;

    for(int i = 0; i < 4; i++){ 
        RC_PmerNext = RC_Kmer + m_RC_shifted[i]; //YX where X is the reverse complement of the kmer and Y is the reverse complement of the last nucleotide of the pmer
        if(m_sequence.present((PmerNext < RC_PmerNext)?PmerNext:RC_PmerNext)){ //if the pmer is present 
            res ^= m_correspondingBitValueForNextPmers[i];
        }
        PmerPrev = Kmer + m_nucleotides_shifted[i]; //equals to AX then CX then GX then TX where X is the Kmer
        RC_PmerPrev = RC_Kmer_ShiftedOf2Bits + m_RC[i];
        if(m_sequence.present((PmerPrev < RC_PmerPrev)?PmerPrev:RC_PmerPrev)){
            res ^= m_correspondingBitValueForPrevPmers[i];
        }
        PmerNext++; //equals to XA then XC then XG then XT
    }
    
    return res;
}



/**
 * Returns the compressed sequence.
 * @return an sd_vector
 */
Intermediate ConwayBromageBM::getSequence(){
    return m_sequence;
}




/**
 * @brief serialize the CBBM object 
 * 
 * @param output 
 * @return int 
 */
//TODO:serialization need to be redo after change
int ConwayBromageBM::serialize(std::string dirpath){
    std::filesystem::create_directory(dirpath);
    for(int i=0;i<m_sequence.nbv;i++){
        std::ofstream tmp(dirpath+"/"+std::to_string(i));
        serializeaux(tmp,m_sequence.vbv.at(i));
        tmp.close();
    }
}
int ConwayBromageBM::serializeaux(std::ostream& output, bm::bvector<>& sequence){
    unsigned char* buf = 0;
    try{
        bm::serializer<bm::bvector<> > bvs;
        bvs.byte_order_serialization(true);
        bvs.gap_length_serialization(true);
        bm::serializer<bm::bvector<> >::buffer sbuf;
        {
            bvs.serialize(sequence, sbuf);
            buf = sbuf.data();
            auto sz = sbuf.size();
            std::cout<<sz<<std::endl;
            output.write((const char*)sbuf.data(),sz);

        }
    }catch(std::exception& ex){
        std::cerr << ex.what() << std::endl;
        delete [] buf;
        return 1;
    }


    return 0;
}
ConwayBromageBM ConwayBromageBM::deserialize(std::string dirpath,KmerManipulator* km){
    int tmp=std::pow(2,(km->getSize()-24)*2);
    std::vector<bm::bvector<>> tmp2;
    for(int i=0;i<tmp;i++){
        std::ifstream ftmp(dirpath+"/"+std::to_string(i),std::ios::binary);

        bm::bvector bmtmp;//=deserializeaux(ftmp);
        bm::serializer<bm::bvector<> >::buffer sbuf;
        read_bvector(ftmp,bmtmp,sbuf);
        tmp2.push_back(bmtmp);
        ftmp.close();
    }
    Intermediate i(tmp2,tmp);
    return ConwayBromageBM(i,km);
    


}

int ConwayBromageBM::read_bvector(std::ifstream& bv_file,
                  bm::bvector<>& bv,
                  bm::serializer<bm::bvector<> >::buffer& sbuf)
{
    if (!bv_file.good())
        return -1;
    bv_file.seekg(0, ios::end);
    unsigned int len = bv_file.tellg();
    bv_file.seekg(ios::beg);
    if (!bv_file.good())
        return -3;
    if (!len)
        return -2; // 0-len detected (broken file)
    std::cout<<len<<std::endl;
    sbuf.resize(len, false); // resize without content preservation
    bv_file.read((char*) sbuf.data(), std::streamsize(len));
    if (!bv_file.good())
        return -4;
    
    bm::deserialize(bv, sbuf.data());
 
    return 0;
}







