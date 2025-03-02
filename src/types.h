#ifndef TYPES_H
#define TYPES_H

#include <map>
#include <string>

namespace OptimLight
{

    /*Specify what information will be output in the algorithm.
	The value should be assigned to the member variable: "VerbosityLevel" */
    enum VerbosityLevel { LOW, MEDIUM, HIGH, COMPLETE, VerbosityLevelLength};

    // /*The algorithm is stopped when a value (specified by ther parameter) is less than the "Tolerance" (a member variable)
	// The value should be assigned to the member variable: "StopCriterion" and the applicable values are
	// FUNC_REL: |f_k - f_{k+1}| / max(|f_k|, 1)
    // FUNC_ABS: |f_k - f_{k+1}| 
	// GRAD_ABS: \|gf_k\|
	// GRAD_REL: \|gf_k\| / \|gf_0\|*/
    enum StopCriterion { FUNC_REL, FUNC_ABS, GRAD_ABS, GRAD_REL,  StopCriterionLength};


    /*  Condition type for line search
    ARMIJO: The Armijo-Goldstein condition.[DS83 Algorithm A6.3.1] combined with nonmontone line search
    WOLFE: The weak Wolfe condition [DS83 Algorithm A6.3.1mod]
    STRONG_WOLFE: The strong Wolfe condition [NW06 Algorithm 3.5]
    */
    enum ConditionType {ARMIJO, WOLFE,  STRONG_WOLFE, ConditionTypeLength};
                    
/**
 * @brief Possible statuses for the optimizer to indicate why
 * optimization exited.
 *
 * Negative value indicate errors while positive values indicate convergence.
 */
enum class ExitStatus {
    /** Optimizer did not run */
    did_not_run = 0,
    /** Reached maximum number of allowed iterations */
    max_iter = -1,
    /** Expected to reach maximum allowed time in next iteration */
    max_time = -2,
    /** Encountered non-finite fval/grad/hess */
    not_finite = -3,
    /** Exceeded specified boundaries */
    exceeded_boundary = -4,
    /** Trust region radius too small to proceed */
    delta_too_small = -5,
    /** Converged according to fval difference */
    ftol = 1,
    /** Converged according to x difference */
    xtol = 2,
    /** Converged according to gradient norm */
    gtol = 3
};

const std::map<ExitStatus, std::string> exit_status_to_str{
    {ExitStatus::did_not_run, "did_not_run"},
    {ExitStatus::max_iter, "max_iter"},
    {ExitStatus::max_time, "max_time"},
    {ExitStatus::not_finite, "not_finite"},
    {ExitStatus::exceeded_boundary, "exceeded_boundary"},
    {ExitStatus::delta_too_small, "delta_too_small"},
    {ExitStatus::ftol, "ftol"},
    {ExitStatus::xtol, "xtol"},
    {ExitStatus::gtol, "gtol"},
};


} // namespace OptimLight





#endif // TYPES_H