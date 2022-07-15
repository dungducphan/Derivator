#include "Derivator.h"

Derivator::Derivator() : mFunc(nullptr),
                         mStepSize(1E-4) {}

Derivator::Derivator(double (*func)(std::vector<double>, std::vector<double>)) : mStepSize(1E-4) {
    mFunc = func;
}

Derivator::~Derivator() = default;

void Derivator::SetFunction(double (*func)(std::vector<double>, std::vector<double>)) {
    mFunc = func;
}

Derivator::Derivator(double (*func)(std::vector<double>, std::vector<double>), double stepSize) {
    mFunc = func;
    mStepSize = stepSize;
}

void Derivator::SetStepSize(double stepSize) {
    mStepSize = stepSize;
}

std::vector<double> Derivator::StepArgs(const std::vector<double>& args,
                                        unsigned int varIndex,
                                        double stepSize) {
    std::vector<double> ret;
    for (unsigned int i = 0; i < args.size(); i++) {
        i == varIndex ? ret.push_back(args.at(i) + stepSize) : ret.push_back(args.at(i));
    }

    return ret;
}

std::vector<double> Derivator::GetGradient(const std::vector<double>& args,
                                           const std::vector<double>& params) {

    std::vector<double> ret;
    for (unsigned int i = 0; i < args.size(); i++) {
        double gradient   = - 1 * mFunc(StepArgs(args, i, + 2 * mStepSize), params)
                            + 8 * mFunc(StepArgs(args, i, + 1 * mStepSize), params)
                            - 8 * mFunc(StepArgs(args, i, - 1 * mStepSize), params)
                            + 1 * mFunc(StepArgs(args, i, - 2 * mStepSize), params);
        ret.push_back(gradient / (12 * mStepSize));
    }

    return ret;
}

std::vector<std::vector<double>> Derivator::GetHessian(const std::vector<double>& args,
                                                       const std::vector<double>& params) {
    std::vector<std::vector<double>> ret;

    for (unsigned int i = 0; i < args.size(); ++i) {
        std::vector<double> row;
        for (unsigned int j = 0; j < args.size(); ++j) {
            double hessian =
                    +  1 * mFunc(StepArgs(StepArgs(args, i, + 2 * mStepSize), j, + 2 * mStepSize), params)
                    -  8 * mFunc(StepArgs(StepArgs(args, i, + 2 * mStepSize), j, + 1 * mStepSize), params)
                    +  8 * mFunc(StepArgs(StepArgs(args, i, + 2 * mStepSize), j, - 1 * mStepSize), params)
                    -  1 * mFunc(StepArgs(StepArgs(args, i, + 2 * mStepSize), j, - 2 * mStepSize), params)

                    -  8 * mFunc(StepArgs(StepArgs(args, i, + 1 * mStepSize), j, + 2 * mStepSize), params)
                    + 64 * mFunc(StepArgs(StepArgs(args, i, + 1 * mStepSize), j, + 1 * mStepSize), params)
                    - 64 * mFunc(StepArgs(StepArgs(args, i, + 1 * mStepSize), j, - 1 * mStepSize), params)
                    +  8 * mFunc(StepArgs(StepArgs(args, i, + 1 * mStepSize), j, - 2 * mStepSize), params)

                    +  8 * mFunc(StepArgs(StepArgs(args, i, - 1 * mStepSize), j, + 2 * mStepSize), params)
                    - 64 * mFunc(StepArgs(StepArgs(args, i, - 1 * mStepSize), j, + 1 * mStepSize), params)
                    + 64 * mFunc(StepArgs(StepArgs(args, i, - 1 * mStepSize), j, - 1 * mStepSize), params)
                    -  8 * mFunc(StepArgs(StepArgs(args, i, - 1 * mStepSize), j, - 2 * mStepSize), params)

                    -  1 * mFunc(StepArgs(StepArgs(args, i, - 2 * mStepSize), j, + 2 * mStepSize), params)
                    +  8 * mFunc(StepArgs(StepArgs(args, i, - 2 * mStepSize), j, + 1 * mStepSize), params)
                    -  8 * mFunc(StepArgs(StepArgs(args, i, - 2 * mStepSize), j, - 1 * mStepSize), params)
                    +  1 * mFunc(StepArgs(StepArgs(args, i, - 2 * mStepSize), j, - 2 * mStepSize), params);

            row.push_back(hessian / (12 * 12 * mStepSize * mStepSize));
        }
        ret.push_back(row);
    }

    return ret;
}


