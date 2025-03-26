#ifndef TYPES_H
#define TYPES_H

#include <map>
#include <string>

namespace OptimLight
{

    typedef enum {
       ALGORITHM_SD, // Steepest descent 
       ALGORITHM_CG, // Conjugate gradient
       ALGORITHM_BFGS, // Broyden-Fletcher-Goldfarb-Shanno
       ALGORITHM_LBFGS, // Limited-memory BFGS
       ALGORITHM_TNEWTON, // Truncated Newton
       ALGORITHM_NEWTON, // Newton
       ALGORITHM_TR, // Trust region
       ALGORITHM_LENGTH
    } algorithm;

    typedef enum {
        RESULT_DIDNOTRUN = 0, // did not run at all
        RESULT_MAXITER_REACHED = -1, //Reached maximum number of allowed iterations 
        RESULT_MAXTIME_REACHED = -2, //Expected to reach maximum allowed time in next iteration
        RESULT_EXCEEDED_BOUNDARY = -3, // Exceeded specified boundaries 
        RESULT_INFINITE = -4, // Encountered non-finite fval/grad/hess
       
        RESULT_SUCCESS = 1, // Success
        RESULT_FTOL_REACHED = 2, // Converged according to fval difference
        RESULT_FTOLREL_REACHED = 3, // Converged according to  |f_k - f_{k+1}| / max(|f_k|, 1)
        RESULT_GTOL_REACHED = 4, // Converged according to gradient norm 
        RESULT_GTOLREL_REACHED = 5, // Converged according to  \|gf_k\| / \|gf_0\|*/

        RESULT_RESULT_LENGTH, // Success
    } result;

    typedef enum {
        VERBOSE_LOW,
        VERBOSE_HIGH,
        VERBOSE_COMPLETE,
        VERBOSE_LENGTH
    } verbose;

    typedef enum {
        LINESEARCH_ARMIJO,
        LINESEARCH_WOLFE,
        LINESEARCH_STRONG_WOLFE,
        LINESEARCH_EXACT,
        LINESEARCH_LENGTH
    } linesearch; // Line search types


    std::string  algorithm_to_string(algorithm algorithm);
    algorithm algorithm_from_string(std::string & name);


};  // namespace OptimLight



#endif // TYPES_H