#ifndef EUCLIDEAN_HPP
#define EUCLIDEAN_HPP

#include "manifold.hpp"

namespace OptimLight {

class Euclidean : public Manifold {
public:
    // Constructor for n-dimensional Euclidean space
    explicit Euclidean(int n_, int p_, bool is_complex = false) 
        : n(n_), p(p_), is_complex_(is_complex) {
        if (n <= 0) {
            throw std::runtime_error("Dimension must be positive");
        }
        name = "Euclidean(" + std::to_string(n) + ")";
        empty = ManifoldVector(n, p, is_complex);
    }

    // Copy constructor
    Euclidean(const Euclidean& other)
        : p(other.p), n(other.n), is_complex_(other.is_complex_) 
    {
        name = other.name;
        empty = other.empty;
    }

    // Manifold operations
    double metric(const ManifoldPoint& x,
                 const ManifoldVector& etax,
                 const ManifoldVector& xix) const override;

    ManifoldVector projection(const ManifoldPoint& x,
                            const ManifoldVector& etax) const override;

    ManifoldPoint retraction(const ManifoldPoint& x,
                            const ManifoldVector& etax) const override;

    ManifoldVector vector_transport(const ManifoldPoint& x,
                                  const ManifoldVector& etax,
                                  const ManifoldPoint& y,
                                  const ManifoldVector& xix) const override;

    // Dimension getters
    int dimension() const override 
    { 
        return n*p;
     }
    int intrinsic_dimension() const override 
    { 
        return n*p;
    }

    bool is_complex() const { return is_complex_; }

    int p;  // Number of rows
    int n;  // Number of columns
    bool is_complex_;

private:
    void check_dimensions(const ManifoldPoint& x, const std::string& name) const;

};

} // namespace OptimLight

#endif // EUCLIDEAN_HPP