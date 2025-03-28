#include "stacked_manifold.hpp"
#include <stdexcept>

namespace OptimLight {

void StackedManifold::check_dimensions(const ManifoldPoint& x, 
                                     const std::string& name) const {
    int expected_rows = base_->empty.n_rows() * num_copies_;
    if (x.n_rows() != expected_rows) {
        throw std::runtime_error(name + " dimension mismatch. Expected " + 
                               std::to_string(expected_rows) + " rows, got " + 
                               std::to_string(x.n_rows()));
    }
}

ManifoldPoint StackedManifold::extract_submanifold(const ManifoldPoint& x, 
                                                  int index) const {
    int dim = base_->empty.n_rows();
    int start = index * dim;
    return x.submat(start, 0, start + dim - 1, x.n_cols() - 1);
}

ManifoldVector StackedManifold::extract_subvector(const ManifoldVector& v, 
                                                 int index) const {
    int dim = base_->empty.n_rows();
    int start = index * dim;
    return v.submat(start, 0, start + dim - 1, v.n_cols() - 1);
}

double StackedManifold::metric(const ManifoldPoint& x, 
                             const ManifoldVector& etax, 
                             const ManifoldVector& xix) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    check_dimensions(xix, "xix");

    double sum = 0.0;
    for (int i = 0; i < num_copies_; ++i) {
        ManifoldPoint x_i = extract_submanifold(x, i);
        ManifoldVector etax_i = extract_subvector(etax, i);
        ManifoldVector xix_i = extract_subvector(xix, i);
        
        sum += base_->metric(x_i, etax_i, xix_i);
    }
    return sum;
}

ManifoldVector StackedManifold::projection(const ManifoldPoint& x, 
                                         const ManifoldVector& etax) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");

    ManifoldVector result(x.n_rows(), x.n_cols(), x.is_complex());
    int dim = base_->empty.n_rows();

    for (int i = 0; i < num_copies_; ++i) {
        ManifoldPoint x_i = extract_submanifold(x, i);
        ManifoldVector etax_i = extract_subvector(etax, i);
        
        ManifoldVector proj_i = base_->projection(x_i, etax_i);
        int start = i * dim;
        result.submat(start, 0, start + dim - 1, result.n_cols() - 1) = proj_i;
    }
    return result;
}

ManifoldPoint StackedManifold::retraction(const ManifoldPoint& x, 
                                        const ManifoldVector& etax) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");

    ManifoldPoint result(x.n_rows(), x.n_cols(), x.is_complex());
    int dim = base_->empty.n_rows();

    for (int i = 0; i < num_copies_; ++i) {
        ManifoldPoint x_i = extract_submanifold(x, i);
        ManifoldVector etax_i = extract_subvector(etax, i);
        
        ManifoldPoint retr_i = base_->retraction(x_i, etax_i);
        int start = i * dim;
        result.submat(start, 0, start + dim - 1, result.n_cols() - 1) = retr_i;
    }
    return result;
}

ManifoldVector StackedManifold::vector_transport(const ManifoldPoint& x, 
                                               const ManifoldVector& etax,
                                               const ManifoldPoint& y, 
                                               const ManifoldVector& xix) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    check_dimensions(y, "y");
    check_dimensions(xix, "xix");

    ManifoldVector result(x.n_rows(), x.n_cols(), x.is_complex());
    int dim = base_->empty.n_rows();

    for (int i = 0; i < num_copies_; ++i) {
        ManifoldPoint x_i = extract_submanifold(x, i);
        ManifoldVector etax_i = extract_subvector(etax, i);
        ManifoldPoint y_i = extract_submanifold(y, i);
        ManifoldVector xix_i = extract_subvector(xix, i);
        
        ManifoldVector trans_i = base_->vector_transport(x_i, etax_i, y_i, xix_i);
        int start = i * dim;
        result.submat(start, 0, start + dim - 1, result.n_cols() - 1) = trans_i;
    }
    return result;
}

} // namespace OptimLight