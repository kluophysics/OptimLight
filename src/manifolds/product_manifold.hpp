#ifndef PRODUCT_MANIFOLD_H
#define PRODUCT_MANIFOLD_H

#include "manifold.hpp"
#include <vector>

namespace OptimLight {

class ProductManifold : public Manifold {
public:
    // Default constructor
    ProductManifold();
    
    // Constructor with manifolds and their powers
    ProductManifold(const std::vector<const Manifold*>& base_manifolds,
                   const std::vector<int>& powers);
    
    // Copy constructor and assignment
    ProductManifold(const ProductManifold& other);
    ProductManifold& operator=(const ProductManifold& other);
    
    // Destructor
    ~ProductManifold() override;

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

    // Dimension calculations
    virtual int dimension() const ;
    virtual int intrinsic_dimension() const ;

protected:
    Manifold** manifolds;        // Store all kinds of manifolds
    int numoftypes;              // Number of kinds of manifolds
    std::vector<int> powsinterval; // Manifold intervals
    int numoftotalmani;          // Total number of manifolds

private:
    // Helper functions
    void initialize(const std::vector<const Manifold*>& base_manifolds,
                   const std::vector<int>& powers);
    void cleanup();
    void copy_from(const ProductManifold& other);
    
    // Manifold extraction helpers
    ManifoldPoint extract_submanifold(const ManifoldPoint& x, int type_idx, int copy_idx) const;
    std::pair<int, int> get_manifold_dimensions(int type_idx) const;
    void check_dimensions(const ManifoldPoint& x, const std::string& name) const;
};

} // namespace OptimLight
#endif // PRODUCT_MANIFOLD_H