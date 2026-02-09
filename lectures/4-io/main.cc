#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
int main(){
    std::vector<double> x,y;
    double number1, number2;
    std::cout << "stdout stream" << "\n";
    std::cerr << "stderr stream" << "\n"; 

    while(std::cin >> number1 >> number2){
        x.push_back(number1);
        y.push_back(number2);
    }
    std::cout << "x y atan2(x,y)" << "\n";
    for(size_t i = 0; i < x.size(); ++i){
        std::cout << x[i] << " " << y[i] << " " << std::atan2(x[i],y[i]) << "\n";
    }
    std::ifstream myinput("data.txt");
    std::ofstream myoutput("out.txt");
    while(myinput >> number1 >> number2){
        myoutput << number1 << " " << number2 << " " << std::atan2(number1,number2) << "\n";
    }
}