#include <iostream>
#include "gradient_descent_lib/Derivator.h"

// Users need to write this function, but they have to use the exactly same prototype
double myMultivariateFunction(std::vector<double> args, std::vector<double> params) {
    return params[0]*args[0]*args[0] + params[1]*args[1]*args[1];
}

int main() {
    std::vector<double> args;
    args.push_back(1.0);
    args.push_back(3.0);
    std::vector<double> params;
    params.push_back(0.2);
    params.push_back(-5.1);
    double stepSize = 1E-1;

    auto dev = new Derivator(&myMultivariateFunction, stepSize);
    auto grad = dev->GetGradient(args, params);
    std::vector<std::vector<double>> hessian;
    hessian = dev->GetHessian(args, params);

    for (auto& elem : grad) std::cout << elem << std::endl;
    for (unsigned int i = 0; i < hessian.size(); ++i) {
        for (unsigned int j = 0; j < hessian.size(); j++) {
            std::cout << hessian.at(i).at(j) << "\t";
        }
        std::cout << std::endl;
    }

    return 0;
}