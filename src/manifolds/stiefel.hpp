#ifndef STIEFEL_H
#define STIEFEL_H

#include "manifold.hpp"
#include <complex>
#include <armadillo>

namespace OptimLight
{

    enum MetricType { 
        EUCLIDEAN, 
        CANONICAL, 
        MetricTypeLength
    };

    enum RetractionType {
        RT_QF,      // QR-based retraction
        RT_EXP,     // Exponential mapping
        RT_CAYLEY,  // Cayley transform
        RT_POLAR,   // Polar decomposition
        RetractionTypeLength
    };

    enum VectorTransportType {
        VT_PROJECTION,
        VT_PARALLELTRANSLATION,
        VT_DIFFERENTIATED,
        VT_CAYLEY,
        VT_RIGGING,
        VT_PARALLELIZATION,
        VectorTransportTypeLength
    };

    class Stiefel : public Manifold
    {
    public:
        // Constructor with full parameters
        // Stiefel();
        Stiefel(int n_, int p_, bool is_complex = false,
                MetricType metric_type = CANONICAL,
                RetractionType retraction_type = RT_QF,
                VectorTransportType vector_transport_type = VT_PROJECTION)
            : n(n_), p(p_),  is_complex_(is_complex),
              metric_type_(metric_type),
              retraction_type_(retraction_type),
              vector_transport_type_(vector_transport_type) {
            if (p_ <= 0 || n_ <= 0 || p_ > n_) {
                throw std::runtime_error("Invalid Stiefel manifold dimensions p=" 
                    + std::to_string(p) + ", n=" + std::to_string(n));
            }
            name = "Stiefel(" + std::to_string(p) + "," + std::to_string(n) + ")";
            // Initialize empty tangent vector with correct dimensions
            empty = ManifoldVector(n, p, is_complex);
        }

        // Copy constructor
        Stiefel(const Stiefel& other)
            : n(other.n),  p(other.p), is_complex_(other.is_complex_),
              metric_type_(other.metric_type_),
              retraction_type_(other.retraction_type_),
              vector_transport_type_(other.vector_transport_type_) {
            name = other.name;
            empty = other.empty;
        }

        // Assignment operator
        Stiefel& operator=(const Stiefel& other) {
            if (this != &other) {
                p = other.p;
                n = other.n;
                is_complex_ = other.is_complex_;
                metric_type_ = other.metric_type_;
                retraction_type_ = other.retraction_type_;
                vector_transport_type_ = other.vector_transport_type_;
                name = other.name;
                empty = other.empty;
            }
            return *this;
        }

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
            return p * n;
        }

        int intrinsic_dimension() const  {
            return p * n - p * (p + 1) / 2;
        }

        // Setters for different types
        void set_metric_type(MetricType type) { metric_type_ = type; }
        void set_retraction_type(RetractionType type) { retraction_type_ = type; }
        void set_vector_transport_type(VectorTransportType type) { vector_transport_type_ = type; }
        
        // Getters for manifold properties
        bool is_complex() const { return is_complex_; }

        int n;  // Number of rows 
        int p;  // Number of columns

        bool is_complex_;

    private:
        void check_dimensions(const ManifoldPoint& x, const std::string& name) const;
        void check_orthogonality(const ManifoldPoint& x, const std::string& name, double tol = 1e-10) const;
        void symmatu(arma::mat& result, const arma::mat& X) const;

        MetricType metric_type_;
        RetractionType retraction_type_;
        VectorTransportType vector_transport_type_;
    };
}

#endif // STIEFEL_H