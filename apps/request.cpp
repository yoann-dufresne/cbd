#include <ConwayBromageLib.h>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
void randomcontainsrequest(int nb,ConwayBromageBM cb){
    for(int i=0;i<nb;i++){
        cb.contains(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}
void randomcontainsrequest(int nb,ConwayBromageSD cb){
    for(int i=0;i<nb;i++){
        cb.contains(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}

void linearcontainsrequest(int nb,ConwayBromageBM cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}
void linearcontainsrequest(int nb,ConwayBromageSD cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}

void randomsuccessorsrequest(int nb,ConwayBromageBM cb){
    for(int i=0;i<nb;i++){
        cb.successors(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}
void randomsuccessorsrequest(int nb,ConwayBromageSD cb){
    for(int i=0;i<nb;i++){
        cb.successors(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}

void linearsuccessorsrequest(int nb,ConwayBromageBM cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}

void linearsuccessorsrequest(int nb,ConwayBromageSD cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}



//./request k file nb           with k the length of the kmer in the file, nb the number of request
int main(int argc, char* argv[]){
    std::cout<<atoi(argv[1])<<std::endl;
    KmerManipulatorACGT tmp(31);
    std::ifstream f(argv[2], std::ios::in);
    KmerManipulatorACGT km(30);
    ConwayBromageBM test(f,&tmp);
    f.close();
    std::ifstream f2(argv[2], std::ios::in);
    ConwayBromageSD test2(f2,&tmp);
    f2.close();
    for(int i=0;i<60;i++){
        auto start = std::chrono::steady_clock::now();
        linearsuccessorsrequest(atoi(argv[3]),test,806,9);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        
        
        auto start2 = std::chrono::steady_clock::now();
        linearsuccessorsrequest(atoi(argv[3]),test2,806,9);
        auto end2 = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds2 = end2-start2;
        if(elapsed_seconds.count()<elapsed_seconds2.count()){
            std::cout<<"BM better"<<std::endl;
        }else{
            std::cout<<"SD better"<<std::endl;
        }
        std::cout<<"SD : "<<elapsed_seconds2.count()<<std::endl;
        std::cout<<"BM : "<<elapsed_seconds.count()<<std::endl;

    }

}
