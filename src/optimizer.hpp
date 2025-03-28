#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

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

#endif // OPTIMIZER_HPP