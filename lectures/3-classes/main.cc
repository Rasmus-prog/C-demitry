#include <iostream>
#include <vector>
#include "hello.h"
#include "vec3.h"
int main() {
    hello();
    double a = 1.0;
    double b = a;
    if(a==b){std::cout<<"a==b"<<std::endl;} 
    else {std::cout<<"a!=b"<<std::endl;}
    std::vector<double> v {1,2,3};

    for(size_t i =0;i<v.size();i++){
        std::cout<<v[i] << " ";
    }
    std::cout<<std::endl;
    for(auto vi : v){
        std::cout<<vi << " ";
    }
    std::cout<<std::endl;

    for (double vi : v){
        std::cout<<vi << " ";
    }
    std::cout<<std::endl;
    
    for(auto vi : v){vi =6;}
    for(auto vi : v){
        std::cout<<vi << " ";
    }

    std::cout<<std::endl;
    for(auto& vi : v){vi =6;}
    for(auto& vi : v){
        std::cout<<vi << " ";
    }
    std::cout<<std::endl<<"now comes while"<<std::endl;

    size_t i =0;
    while(i< v.size()){
        std::cout<<"v["<<i<<"] = "<<v[i]<<' ';
        i++;
    }
    std::cout<<std::endl;
    std::cout<<"now comes do-while"<<std::endl;
    i=0;
    do{
        std::cout<<"v["<<i<<"] = "<<v[i]<<' ';
        i++;
    }while(i< v.size());
    std::cout<<std::endl;

    std::cout<<"Now comes vec3 class"<<std::endl;
    pp::vec3 p1(1,2,3);
    p1.x = 6;
    std::cout<<"p1 = ("<<p1.x<<","<<p1.y<<","<<p1.z<<")"<<std::endl;


return 0;
}