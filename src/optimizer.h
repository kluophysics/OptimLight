#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "types.h"

namespace OptimLight
{


class OptimizerMethod
{

};


class Optimizer
{
    public:
        Optimizer();
        ~Optimizer();

        // void setMaxIter(int max_iter);
        // void setTol(double tol);
        // void setVerbose(bool verbose);
        void setMethod(OptimizerMethod *method);

        void optimize();

    private:
        int max_iter;
        double tol;
        bool verbose;
        OptimizerMethod * method;
};




} // namespace OptimLight

#endif // OPTIMIZER_H