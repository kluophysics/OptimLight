#ifndef MANIFOLD_H
#define MANIFOLD_H

#include "array.hpp"
#include <tuple>
#include <type_traits>
#include <utility>
#include <string>

namespace OptimLight
{
    using ManifoldPoint = Array;
    using ManifoldVector = Array;
    
    class Manifold
    {
    public:
        virtual ~Manifold() = default;

        virtual double metric(const ManifoldPoint& x, 
                            const ManifoldVector& etax, 
                            const ManifoldVector& xix) const = 0;

        virtual ManifoldVector projection(const ManifoldPoint& x, 
                                        const ManifoldVector& etax) const = 0;

        virtual ManifoldPoint retraction(const ManifoldPoint& x, 
                                       const ManifoldVector& etax) const = 0;

        virtual ManifoldVector vector_transport(const ManifoldPoint& x, 
                                              const ManifoldVector& etax,
                                              const ManifoldPoint& y, 
                                              const ManifoldVector& xix) const = 0;

        // virtual int intrinsic_dimension() const = 0;
        // virtual int dimension() const = 0;

        virtual int intrinsic_dimension() const =0;
        virtual int dimension() const =0;
        std::string name; // name of the manifold
        
        ManifoldVector empty; // empty tangent vector
    };

} // namespace OptimLight

#endif // MANIFOLD_H