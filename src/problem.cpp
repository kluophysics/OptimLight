#include "problem.hpp"

namespace OptimLight
{

double Problem::objective_function(const ManifoldPoint &x) const
{
    // Provide a default implementation or leave it pure virtual
    return 0.0;
}

ManifoldVector Problem::gradient(const ManifoldPoint &x) const
{
    // Provide a default implementation or leave it pure virtual
    return ManifoldVector();
}

void Problem::evaluate_obj_and_grad(const ManifoldPoint &x) const
{
    // Provide a default implementation or leave it pure virtual
}

ManifoldVector Problem::riemannian_gradient(const ManifoldPoint &x) const
{
    // Provide a default implementation or leave it pure virtual
    return ManifoldVector();
}
}