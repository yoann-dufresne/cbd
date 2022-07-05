#include <ConwayBromageLib.h>
#include <sdsl/sd_vector.hpp>
#include <sdsl/vectors.hpp>
#include <bitset>

void randomcontainsrequest(int nb,ConwayBromageBM& cb){
    for(int i=0;i<nb;i++){
        cb.contains(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}
void randomcontainsrequest(int nb,ConwayBromageSD& cb){
    for(int i=0;i<nb;i++){
        cb.contains(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}

void linearcontainsrequest(int nb,ConwayBromageBM& cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}
void linearcontainsrequest(int nb,ConwayBromageSD& cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}

void randomsuccessorsrequest(int nb,ConwayBromageBM& cb){
    for(int i=0;i<nb;i++){
        cb.successors(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}
void randomsuccessorsrequest(int nb,ConwayBromageSD& cb){
    for(int i=0;i<nb;i++){
        cb.successors(rand()%((uint64_t)(std::pow(4,(cb.getKmerManipulator()->getSize()-1) ) ) ) );
    }
}

void linearsuccessorsrequest(int nb,ConwayBromageBM& cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}

void linearsuccessorsrequest(int nb,ConwayBromageSD& cb,int start=0,int jump=1){
    for(int i=start;i<nb+start*jump;i+=jump){
        cb.contains(i);
    }
}
void successivecontainsrequest(int nb,ConwayBromageSD& cb,uint64_t kmer){
    std::cout<<"test"<<std::endl;
    uint64_t mask=0;
    for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
    }

    //choose randomly one of the next successor to 
    for(int i=0;i<nb;i++){
        cb.contains(kmer);
        kmer=(kmer<<2)&mask;
        int r=rand()%4;
        kmer+=r;
        bitset<64> tmp(kmer);
        std::cout<<tmp<<std::endl;
           
    }
    
}
void successivecontainsrequest(int nb,ConwayBromageBM& cb,uint64_t kmer){
    uint64_t mask=0;
    for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
    }

    //choose randomly one of the next successor to 
    for(int i=0;i<nb;i++){
        cb.contains(kmer);
        kmer=(kmer<<2)&mask;
        int r=rand()%4;
        kmer+=r;
    }
}

void successivesuccessorrequest(int nb,ConwayBromageBM& cb,uint64_t kmer){
    uint64_t mask=0;
    for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
    }

    //choose randomly one of the next successor to 
    for(int i=0;i<nb;i++){
        cb.successors(kmer);
        kmer=(kmer<<2)&mask;
        int r=rand()%4;
        kmer+=r;
    }
}
void successivesuccessorrequest(int nb,ConwayBromageSD& cb,uint64_t kmer){
    uint64_t mask=0;
    for(int i=0;i<cb.getKmerManipulator()->getSize()-1;i++){
        mask=mask|(uint64_t)1<<i*2|(uint64_t)1<<i*2+1;
    }

    //choose randomly one of the next successor to 
    for(int i=0;i<nb;i++){
        cb.successors(kmer);
        kmer=(kmer<<2)&mask;
        int r=rand()%4;
        kmer+=r;
    }
}

std::string toBinary(int n)
{
    std::string r;
    while(n!=0) {r=(n%2==0 ?"0":"1")+r; n/=2;}
    return r;
}

/*
./ request file k nb succesive type firstkmer
file the file with the kmer(for now only the one needed to make the object, add the serialized version later)
k the length of the kmer(max31)
nb the number of the request to test
successive, a boolean that indicate the type of request tested(0 random 1 successive)
type contains or successor to test
library BM or SDSL
firstkmer the first kmer to test and for successive its the one that is used as a base for the next one 
*/
int main(int argc, char* argv[]){
    std::string first("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    if((argc!=8)&&(argc!=7)){
        std::cout<<argc<<std::endl;
        std::cout<<"./main file k nb random type library first\n";
        std::cout<<"file the file wherer the kmer are\n";
        std::cout<<"k the length of the kmer\n";
        std::cout<<"nb the number of the request to test\n";
        std::cout<<"random or successive kmer to test\n";
        std::cout<<"type, contain or successor\n";
        std::cout<<"library BM or SDSL\n";
        std::cout<<"firstkmer optional the first kmer to test and for successive its the one that is used as a base for the next one ";
        std::cout<<std::endl;
        exit(1);
    }
    if(argc!=7){
        first=argv[8];
    }
    std::ifstream f(argv[1], std::ios::in);
    KmerManipulatorACGT km(atoi(argv[2]));    
    KmerManipulatorACGT k1(atoi(argv[2])-1);
    uint64_t first1=k1.encode(first);
    if(argv[7]=="BM"){
        ConwayBromageBM cbd(f,&km);
        auto start = std::chrono::steady_clock::now();
        if(argv[4]=="successive"){
            if(argv[5]=="contains"){
                successivecontainsrequest(atoi(argv[3]),cbd,first1);
            }else{
                successivesuccessorrequest(atoi(argv[3]),cbd,first1);
            }
        }else{
            if(argv[5]=="contains"){
                randomcontainsrequest(atoi(argv[3]),cbd);
            }else{
                randomsuccessorsrequest(atoi(argv[3]),cbd);
            }
        }
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout<<elapsed_seconds.count()<<std::endl;
    }else{
        ConwayBromageSD cbd(f,&km);
        auto start = std::chrono::steady_clock::now();
        if(argv[4]=="successive"){
            if(argv[5]=="contains"){
                successivecontainsrequest(atoi(argv[3]),cbd,first1);
            }else{
                successivesuccessorrequest(atoi(argv[3]),cbd,first1);
            }
        }else{
            if(argv[5]=="contains"){
                randomcontainsrequest(atoi(argv[3]),cbd);
            }else{
                randomsuccessorsrequest(atoi(argv[3]),cbd);
            }
        }
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout<<elapsed_seconds.count()<<std::endl;

    }
    

}
