#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "types.hpp"

namespace OptimLight
{

class Optimizer
{
    public:
        Optimizer();
        ~Optimizer();


    private:
        int max_iter;
        double tol;
        bool verbose;
        OptimizerMethod * method;
};




} // namespace OptimLight

#endif // OPTIMIZER_H