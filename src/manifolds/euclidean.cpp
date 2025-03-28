#include "euclidean.hpp"
#include <stdexcept>

namespace OptimLight {

void Euclidean::check_dimensions(const ManifoldPoint& x, const std::string& name) const {
    if (x.n_rows() != n || x.n_cols() != p) {
        throw std::runtime_error(name + " has wrong dimensions. Expected " + 
                               std::to_string(n) + "x" + std::to_string(p) + 
                               ", got " + std::to_string(x.n_rows()) + "x" + 
                               std::to_string(x.n_cols()));
    }
}

double Euclidean::metric(const ManifoldPoint& x,
                        const ManifoldVector& etax,
                        const ManifoldVector& xix) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    check_dimensions(xix, "xix");

    if (is_complex_) {
        arma::cx_mat etax_cx = etax.as_complex();
        arma::cx_mat xix_cx = xix.as_complex();
        return std::real(arma::accu(arma::conj(etax_cx) % xix_cx));
    } else {
        return arma::accu(etax.real() % xix.real());
    }
}

ManifoldVector Euclidean::projection(const ManifoldPoint& x,
                                   const ManifoldVector& etax) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    // In Euclidean space, projection is identity
    return etax;
}

ManifoldPoint Euclidean::retraction(const ManifoldPoint& x,
                                   const ManifoldVector& etax) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    // Simple addition in Euclidean space
    return x + etax;
}

ManifoldVector Euclidean::vector_transport(const ManifoldPoint& x,
                                         const ManifoldVector& etax,
                                         const ManifoldPoint& y,
                                         const ManifoldVector& xix) const {
    check_dimensions(x, "x");
    check_dimensions(y, "y");
    check_dimensions(etax, "etax");
    check_dimensions(xix, "xix");
    // In Euclidean space, vector transport is identity
    return xix;
}

} // namespace OptimLight