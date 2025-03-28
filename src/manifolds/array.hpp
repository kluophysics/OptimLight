#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <complex>
#include <memory>
#include <stdexcept>
#include <armadillo>

namespace OptimLight {

class Array {
private:
    arma::mat real_;          // Real part storage
    arma::mat imag_;          // Imaginary part storage
    const bool is_complex_;   // Flag for complex data

public:
    // Constructors
    Array() : real_(), imag_(), is_complex_(false) {}
    Array(size_t rows, size_t cols, bool complex = false) 
        : real_(rows, cols), imag_(complex ? arma::mat(rows, cols) : arma::mat()), 
          is_complex_(complex) {}
    
    // Copy constructors from armadillo types
    explicit Array(const arma::mat& m, bool complex = false) 
        : real_(m), imag_(complex ? arma::mat(m.n_rows, m.n_cols) : arma::mat()), 
          is_complex_(complex) {}
    
    explicit Array(const arma::cx_mat& m) 
        : real_(arma::real(m)), imag_(arma::imag(m)), is_complex_(true) {}
    
    // Constructor for complex matrix from real and imaginary parts
    Array(const arma::mat& real_part, const arma::mat& imag_part) 
        : real_(real_part), imag_(imag_part), is_complex_(true) {
        if (real_part.n_rows != imag_part.n_rows || 
            real_part.n_cols != imag_part.n_cols) {
            throw std::runtime_error("Dimension mismatch between real and imaginary parts");
        }
    }

    // Copy constructor
    Array(const Array& other) 
        : real_(other.real_), imag_(other.imag_), is_complex_(other.is_complex_) {}
    
    // Dimension info
    size_t n_rows() const { return real_.n_rows; }
    size_t n_cols() const { return real_.n_cols; }
    size_t n_elem() const { return real_.n_elem; }
    
    // Data access
    const arma::mat& real() const { return real_; }
    const arma::mat& imag() const { return imag_; }
    bool is_complex() const { return is_complex_; }
    
    // Complex conversions and access
    arma::cx_mat as_complex() const {
        if (is_complex_) {
            return arma::cx_mat(real_, imag_);
        }
        return arma::cx_mat(real_, arma::zeros(n_rows(), n_cols()));
    }
    
    // Matrix access
    arma::mat& as_mat() { 
        if (is_complex_) throw std::runtime_error("Cannot access complex matrix as real");
        return real_; 
    }
    const arma::mat& as_mat() const { 
        if (is_complex_) throw std::runtime_error("Cannot access complex matrix as real");
        return real_; 
    }

    // Submatrix extraction
    Array submat(size_t first_row, size_t first_col, 
                size_t last_row, size_t last_col) const {
        if (is_complex_) {
            return Array(real_.submat(first_row, first_col, last_row, last_col),
                        imag_.submat(first_row, first_col, last_row, last_col));
        }
        return Array(real_.submat(first_row, first_col, last_row, last_col), false);
    }

    // Submatrix assignment
    void submat(size_t first_row, size_t first_col, 
               size_t last_row, size_t last_col, 
               const Array& X) {
        if (is_complex_ != X.is_complex_) {
            throw std::runtime_error("Complex type mismatch in submatrix assignment");
        }
        real_.submat(first_row, first_col, last_row, last_col) = X.real_;
        if (is_complex_) {
            imag_.submat(first_row, first_col, last_row, last_col) = X.imag_;
        }
    }

    // Operators (declarations)
    Array& operator=(const Array& other);
    Array operator+(const Array& other) const;
    Array operator-(const Array& other) const;
    Array operator*(const Array& other) const;
    Array operator/(const Array& other) const;
    Array& operator+=(const Array& other);
    Array& operator-=(const Array& other);
    Array& operator*=(const Array& other);
    Array& operator/=(const Array& other);
    bool operator==(const Array& other) const;
    bool operator!=(const Array& other) const;

    Array operator-() const;  // Unary negation
    Array operator*(double scalar) const;  // Scalar multiplication
    Array operator*(std::complex<double> scalar) const;  // Complex scalar multiplication
    friend Array operator*(double scalar, const Array& arr);  // Left scalar multiplication
    friend Array operator*(std::complex<double> scalar, const Array& arr);  // Left complex scalar multiplication

    // Add scalar compound assignments
    Array& operator*=(double scalar);
    Array& operator*=(std::complex<double> scalar);

private:
    void check_dimensions(const Array& other, const std::string& op) const;
    void check_mult_dimensions(const Array& other) const;
};

} // namespace OptimLight
#endif // ARRAY_HPP