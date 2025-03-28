#include "array.hpp"

namespace OptimLight {

void Array::check_dimensions(const Array& other, const std::string& op) const {
    if (n_rows() != other.n_rows() || n_cols() != other.n_cols()) {
        throw std::runtime_error("Matrix dimension mismatch in operator " + op);
    }
}

void Array::check_mult_dimensions(const Array& other) const {
    if (n_cols() != other.n_rows()) {
        throw std::runtime_error("Matrix multiplication dimension mismatch");
    }
}

Array& Array::operator=(const Array& other) {
    if (this != &other) {
        real_ = other.real_;
        imag_ = other.imag_;
        const_cast<bool&>(is_complex_) = other.is_complex_;
    }
    return *this;
}

Array Array::operator+(const Array& other) const {
    check_dimensions(other, "+");
    if (is_complex_ || other.is_complex_) {
        arma::cx_mat result = as_complex() + other.as_complex();
        return Array(result);  // Use explicit complex constructor
    }
    return Array(real_ + other.real_, false);
}

Array Array::operator-(const Array& other) const {
    check_dimensions(other, "-");
    if (is_complex_ || other.is_complex_) {
        arma::cx_mat result = as_complex() - other.as_complex();
        return Array(result);  // Use explicit complex constructor
    }
    return Array(real_ - other.real_, false);
}

Array Array::operator*(const Array& other) const {
    check_mult_dimensions(other);
    if (is_complex_ || other.is_complex_) {
        arma::cx_mat result = as_complex() * other.as_complex();
        return Array(result);  // Use explicit complex constructor
    }
    return Array(real_ * other.real_, false);
}

Array Array::operator/(const Array& other) const {
    check_dimensions(other, "/");
    if (is_complex_ || other.is_complex_) {
        arma::cx_mat result = as_complex() / other.as_complex();
        return Array(result);  // Use explicit complex constructor
    }
    return Array(real_ / other.real_, false);
}

Array& Array::operator+=(const Array& other) {
    check_dimensions(other, "+=");
    if (is_complex_ || other.is_complex_) {
        arma::cx_mat result = as_complex() + other.as_complex();
        real_ = arma::real(result);
        imag_ = arma::imag(result);
        const_cast<bool&>(is_complex_) = true;
    } else {
        real_ += other.real_;
    }
    return *this;
}

bool Array::operator==(const Array& other) const {
    if (n_rows() != other.n_rows() || n_cols() != other.n_cols() || 
        is_complex_ != other.is_complex_) {
        return false;
    }
    if (is_complex_) {
        return arma::approx_equal(as_complex(), other.as_complex(), "absdiff", 1e-8);
    }
    return arma::approx_equal(real_, other.real_, "absdiff", 1e-8);
}

bool Array::operator!=(const Array& other) const {
    return !(*this == other);
}

Array Array::operator-() const {
    if (is_complex_) {
        return Array(-real_, -imag_);
    }
    return Array(-real_, false);
}

Array Array::operator*(double scalar) const {
    if (is_complex_) {
        return Array(scalar * real_, scalar * imag_);
    }
    return Array(scalar * real_, false);
}

Array Array::operator*(std::complex<double> scalar) const {
    arma::cx_mat result = scalar * as_complex();
    return Array(result);  // This constructor handles complex conversion
}

Array& Array::operator*=(double scalar) {
    real_ *= scalar;
    if (is_complex_) {
        imag_ *= scalar;
    }
    return *this;
}

Array& Array::operator*=(std::complex<double> scalar) {
    if (!is_complex_) {
        const_cast<bool&>(is_complex_) = true;
        imag_ = arma::zeros(real_.n_rows, real_.n_cols);
    }
    arma::cx_mat result = scalar * as_complex();
    real_ = arma::real(result);
    imag_ = arma::imag(result);
    return *this;
}

Array operator*(double scalar, const Array& arr) {
    return arr * scalar;
}

Array operator*(std::complex<double> scalar, const Array& arr) {
    return arr * scalar;
}

} // namespace OptimLight