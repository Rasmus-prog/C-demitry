#include<iostream>
#include<string>
#include<vector>
int main(int argc, char** argv){
    std::string arg;
    double dr = 0.1, rmax = 10;
    std::vector<double> xvec;
    for(int i = 0; i < argc; ++i){
        arg = argv[i];
        // std::cout << arg << "\n";
        if(arg == "-dr" && i+1 < argc) dr= std::stod(argv[i+1]);
        if(arg == "-rmax" && i+1 < argc) rmax= std::stod(argv[i+1]);
        if(arg == "-n" && i+1 < argc){
            double number=std::stod(argv[i+1]);
            xvec.push_back(number);
        }
    }
    std::cout << "dr = " << dr << "\n";
    std::cout << "rmax = " << rmax << "\n";
    std::cout << "xvec = ";
    for(size_t i = 0; i < xvec.size(); ++i){
        std::cout << xvec[i] << " ";
    }
    std::cout << "\n";
    return 0;

}