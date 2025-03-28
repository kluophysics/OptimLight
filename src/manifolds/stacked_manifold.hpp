#ifndef STACKED_MANIFOLD_HPP
#define STACKED_MANIFOLD_HPP

#include "manifold.hpp"
#include "array.hpp"
#include <string>
#include <vector>
#include <memory>

namespace OptimLight {

class StackedManifold : public Manifold {
public:
    // Default constructor
    StackedManifold() : base_(nullptr), num_copies_(0) {
        name = "StackedManifold";
        empty = ManifoldVector(0, 0, false); // Initialize empty vector
    }

    // Constructor taking base manifold and number of copies
    StackedManifold(const Manifold* base_manifold, int num_copies) 
        : base_(base_manifold), num_copies_(num_copies) {
        if (num_copies < 1) {
            throw std::runtime_error("Number of copies must be positive");
        }
        name = "StackedManifold(" + base_manifold->name + " x " + 
               std::to_string(num_copies) + ")";
        
        // Initialize empty vector with proper dimensions
        int base_dim = base_manifold->empty.n_rows(); // vertically stacked 
        empty = ManifoldVector(base_dim * num_copies, 
                             base_manifold->empty.n_cols(),
                             base_manifold->empty.is_complex());
    }

    // Copy constructor
    StackedManifold(const StackedManifold& other)
        : base_(other.base_),
          num_copies_(other.num_copies_) {
        name = other.name;
        empty = other.empty;
    }

    // Assignment operator
    StackedManifold& operator=(const StackedManifold& other) {
        if (this != &other) {
            base_ = other.base_;
            num_copies_ = other.num_copies_;
            name = other.name;
            empty = other.empty;
        }
        return *this;
    }

    // Override base class virtual functions
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

    int dimension() const  {
        return base_->dimension() * num_copies_;
    }

    int intrinsic_dimension() const  {
        return base_->intrinsic_dimension() * num_copies_;
    }

    // Getter for number of copies
    int get_num_copies() const { return num_copies_; }

private:
    const Manifold* base_;    // Base manifold to be stacked
    int num_copies_;          // Number of copies
    
    // Helper functions
    void check_dimensions(const ManifoldPoint& x, const std::string& name) const;
    ManifoldPoint extract_submanifold(const ManifoldPoint& x, int index) const;
    ManifoldVector extract_subvector(const ManifoldVector& v, int index) const;
};

} // namespace OptimLight
#endif // STACKED_MANIFOLD_HPP