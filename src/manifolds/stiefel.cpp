#include "stiefel.hpp"
#include <cassert>
#include <armadillo>
#include <stdexcept>

namespace OptimLight {

void Stiefel::check_dimensions(const ManifoldPoint& x, const std::string& name) const {
    if (x.n_rows() != n || x.n_cols() !=p) {
        throw std::runtime_error(name + " has wrong dimensions. Expected " + 
                               std::to_string(n) + "x" + std::to_string(p) + 
                               ", got " + std::to_string(x.n_rows()) + "x" + 
                               std::to_string(x.n_cols()));
    }
}

void Stiefel::check_orthogonality(const ManifoldPoint& x, const std::string& name, double tol) const {
    arma::mat I = arma::eye(p, p);
    if (x.is_complex()) {
        arma::cx_mat X = x.as_complex();
        double res = arma::norm(X.t() * X- arma::cx_mat(I, arma::mat(p, p, arma::fill::zeros)), "fro");
        if ( res > tol) {
            throw std::runtime_error(name + " is not orthogonal (complex) res =" + std::to_string(res) );
        }
    } else {
        arma::mat X = x.real();
        double res = arma::norm( X.t()  * X  - I, "fro")  ;
        if (res > tol) {
            throw std::runtime_error(name + " is not orthogonal (real), res = " + std::to_string(res));
        }
    }
}

void Stiefel::symmatu(arma::mat& result, const arma::mat& X) const {
    result = (X + X.t()) / 2.0;
    result = arma::trimatu(result);
}

double Stiefel::metric(const ManifoldPoint& x, 
                      const ManifoldVector& etax, 
                      const ManifoldVector& xix) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    check_dimensions(xix, "xix");

    if (x.is_complex() || etax.is_complex() || xix.is_complex()) {
        arma::cx_mat X = x.as_complex();
        arma::cx_mat Z1 = etax.as_complex();
        arma::cx_mat Z2 = xix.as_complex();

        switch (metric_type_) {
            case EUCLIDEAN:
                return std::real(arma::accu(arma::conj(Z1) % Z2));
            case CANONICAL:
                return std::real(arma::accu(arma::conj(Z1) % Z2) - 
                       0.5 * arma::accu(arma::conj(Z1 * X.t()) % (Z2 * X.t())));
            default:
                throw std::runtime_error("Unknown metric type");
        }
    } else {
        arma::mat X = x.real();
        arma::mat Z1 = etax.real();
        arma::mat Z2 = xix.real();

        switch (metric_type_) {
            case EUCLIDEAN:
                return arma::accu(Z1 % Z2);
            case CANONICAL:
                return arma::accu(Z1 % Z2) - 0.5 * arma::accu((Z1 * X.t()) % (Z2 * X.t()));
            default:
                throw std::runtime_error("Unknown metric type");
        }
    }
}

ManifoldVector Stiefel::projection(const ManifoldPoint& x, 
                                 const ManifoldVector& etax) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    check_orthogonality(x, "x");

    if (is_complex_) {
        arma::cx_mat X = x.as_complex();
        arma::cx_mat Z= etax.as_complex();
        arma::cx_mat XZ = X.t() * Z;
        // arma::cx_mat P = Z-   X * arma::symmatu(X.t() * Z); // this is not correct!
        arma::cx_mat P = Z - 0.5* X * ( XZ + XZ.t() );
        return ManifoldVector(P);
    } else {
        arma::mat X = x.as_mat();
        arma::mat Z= etax.as_mat();
        arma::mat XZ = X.t() * Z;
        // arma::mat P = Z- X * arma::symmatu(X.t() * Z); // this is not correct!
        arma::mat P = Z - 0.5* X * ( XZ + XZ.t() );
        return ManifoldVector(P);
    }
}

ManifoldPoint Stiefel::retraction(const ManifoldPoint& x, 
                                const ManifoldVector& etax) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    check_orthogonality(x, "x");

    if (is_complex_) {
        arma::cx_mat X = x.as_complex();
        arma::cx_mat Z= etax.as_complex();

        switch (retraction_type_) {
            case RT_QF: {
                arma::cx_mat Q, R;
                arma::qr_econ(Q, R, X + Z);
                return ManifoldPoint(Q);  // Added complex flag
            }
            case RT_POLAR: {
                arma::cx_mat U, V;
                arma::vec s;
                arma::svd_econ(U, s, V, X + Z);
                arma::cx_mat result = U * V.t();
                // result.brief_print("result");
                return ManifoldPoint(result);
            }
            case RT_CAYLEY: {
                arma::cx_mat W = Z* X.t() - X * Z.t();
                arma::cx_mat I = arma::eye<arma::cx_mat>(n, n);  // Changed p to n
                arma::cx_mat cayley = arma::solve(I + W/2.0, I - W/2.0) * X;
                return ManifoldPoint(cayley);  // Added complex flag
            }
            case RT_EXP: {
                // Get matrices as Armadillo complex objects
                arma::cx_mat X = x.as_complex();
                arma::cx_mat Z = etax.as_complex();
                
                // Step 1: Compute W = X^H Z
                arma::cx_mat W = X.t() * Z;
                
                // Step 2: Compute A = Z - X*W
                arma::cx_mat A = Z - X * W;
                
                // Step 3: Compute QR decomposition of A
                arma::cx_mat Q, R;
                arma::qr_econ(Q, R, A);
                
                // Step 4: Form the block matrix M
                arma::cx_mat M = arma::zeros<arma::cx_mat>(2*p, 2*p);
                M.submat(0, 0, p-1, p-1) = W;
                M.submat(0, p, p-1, 2*p-1) = -R.t() * R;
                M.submat(p, 0, 2*p-1, p-1) = arma::eye<arma::cx_mat>(p, p);
                
                // Step 5: Compute matrix exponential of M
                arma::cx_mat expM = arma::expmat(M);
                
                // Step 6: Extract the top blocks
                arma::cx_mat expM11 = expM.submat(0, 0, p-1, p-1);
                arma::cx_mat expM21 = expM.submat(p, 0, 2*p-1, p-1);
                
                // Step 7: Compute final result Y = X*expM11 + Q*expM21
                arma::cx_mat Y = X * expM11 + Q * expM21;
                
                return ManifoldPoint(Y);
            }
            default:
                throw std::runtime_error("Unsupported retraction type");
        }
    } else {
        arma::mat X = x.as_mat();
        arma::mat Z= etax.as_mat();

        switch (retraction_type_) {
            case RT_QF: {
                arma::mat Q, R;
                arma::qr_econ(Q, R, X + Z);
                return ManifoldPoint(Q);
            }
            case RT_POLAR: {
                arma::mat U, V;
                arma::vec s;
                arma::svd_econ(U, s, V, X + Z);
                arma::mat result = U * V.t();

                // result.brief_print("result");

                return ManifoldPoint(result, false);
            }
            case RT_CAYLEY: {
                arma::mat W = Z* X.t() - X * Z.t();
                arma::mat I = arma::eye(p, p);
                arma::mat cayley = arma::solve(I + W/2.0, I - W/2.0) * X;
                return ManifoldPoint(cayley);
            }
            case RT_EXP: {
                // Get matrices as Armadillo objects
                arma::mat X = x.as_mat();
                arma::mat Z = etax.as_mat();
                
                // Step 1: Compute W = X^T Z (p × p matrix)
                arma::mat W = X.t() * Z;
                
                // Step 2: Compute A = Z - X*W (n × p matrix)
                arma::mat A = Z - X * W;
                
                // Step 3: Compute QR decomposition of A
                arma::mat Q, R;
                arma::qr_econ(Q, R, A);
                
                // Step 4: Form the block matrix M
                arma::mat M = arma::zeros(2*p, 2*p);
                M.submat(0, 0, p-1, p-1) = W;
                M.submat(0, p, p-1, 2*p-1) = -R.t() * R;
                M.submat(p, 0, 2*p-1, p-1) = arma::eye(p, p);
                
                // Step 5: Compute matrix exponential of M
                arma::mat expM = arma::expmat(M);
                
                // Step 6: Extract the top blocks
                arma::mat expM11 = expM.submat(0, 0, p-1, p-1);
                arma::mat expM21 = expM.submat(p, 0, 2*p-1, p-1);
                
                // Step 7: Compute final result Y = X*expM11 + Q*expM21
                arma::mat Y = X * expM11 + Q * expM21;
                
                return ManifoldPoint(Y);
            }
            default:
                throw std::runtime_error("Unsupported retraction type");
        }
    }
}

ManifoldVector Stiefel::vector_transport(const ManifoldPoint& x, 
                                       const ManifoldVector& etax,
                                       const ManifoldPoint& y, 
                                       const ManifoldVector& xix) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    check_dimensions(y, "y");
    check_dimensions(xix, "xix");
    check_orthogonality(x, "x");
    check_orthogonality(y, "y");

    switch (vector_transport_type_) {
        case VT_PROJECTION:
            return projection(y, xix);
        case VT_PARALLELTRANSLATION:
            return projection(y, xix);
        case VT_CAYLEY: {
            if (is_complex_) {
                arma::cx_mat X = x.as_complex();
                arma::cx_mat Z= etax.as_complex();
                arma::cx_mat W = Z* X.t() - X * Z.t();
                arma::cx_mat I = arma::eye<arma::cx_mat>(p, p);
                arma::cx_mat cayley = arma::solve(I + W/2.0, I - W/2.0);
                arma::cx_mat result = cayley * xix.as_complex();
                return ManifoldVector(result);
            } else {
                arma::mat X = x.as_mat();
                arma::mat Z= etax.as_mat();
                arma::mat W = Z* X.t() - X * Z.t();
                arma::mat I = arma::eye(p, p);
                arma::mat cayley = arma::solve(I + W/2.0, I - W/2.0);
                arma::mat result = cayley * xix.as_mat();
                return ManifoldVector(result);
            }
        }
        default:
            throw std::runtime_error("Unsupported vector transport type");
    }
}

} // namespace OptimLight