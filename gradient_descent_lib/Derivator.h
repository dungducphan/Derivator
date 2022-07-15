#ifndef GRADIENT_DESCENT_LIB_DERIVATOR_H
#define GRADIENT_DESCENT_LIB_DERIVATOR_H

#include <vector>

class Derivator {
public:
    Derivator();
    Derivator(double (*) (std::vector<double>, std::vector<double>));
    Derivator(double (*) (std::vector<double>, std::vector<double>), double);
    ~Derivator();

public:
    void SetFunction(double (*) (std::vector<double>, std::vector<double>));
    void SetStepSize(double);

    std::vector<double> GetGradient(const std::vector<double>&, const std::vector<double>&);
    std::vector<std::vector<double>> GetHessian(const std::vector<double>&,
                                                const std::vector<double>&);

    std::vector<double> StepArgs(const std::vector<double>&, unsigned int, double);
    double (*mFunc) (std::vector<double>, std::vector<double>);
    double mStepSize;
};

#endif //GRADIENT_DESCENT_LIB_DERIVATOR_H
