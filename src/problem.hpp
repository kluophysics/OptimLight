#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include "manifolds/manifold.hpp"

namespace OptimLight
{
class Problem
{
public:
        // function value evaluated on x, which has to be defined for each problem
        // virtual double objective(const ManifoldPoint & x) = 0;
        virtual double objective_function(const ManifoldPoint & x)  const=0;

        // Euclidean gradient of f, which has to be defined for each problem
        // virtual ManifoldPoint gradient(const ManifoldPoint& x) =0;
        virtual ManifoldVector  gradient(const ManifoldPoint & x) const=0;

        // evaluate objective function and gradient at the same time.
        virtual void evaluate_obj_and_grad(const ManifoldPoint & x) const=0 ;

        // Riemannian Gradient defined on the manifold 
        virtual ManifoldVector  riemannian_gradient(const ManifoldPoint & x ) const=0 ;

        // set the manifold of the objective function
        virtual void set_manifold(Manifold * mani_in) { manifold_ = mani_in; }

        // get the manifold of the objective function
        inline Manifold * get_manifold() const { return manifold_; }  

        virtual ManifoldVector conditioner(const ManifoldPoint & x, 
                                        const ManifoldVector & eta) const 
        {
            return eta;
        };

        // variable
        Manifold * manifold_; // pointer to hold the manifold of the objective function 

};

} // OptimLight

#endif // PROBLEM_HPP