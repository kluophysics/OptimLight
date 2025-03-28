#include "product_manifold.hpp"
#include <stdexcept>
#include <numeric>
#include <sstream>

namespace OptimLight {

ProductManifold::ProductManifold()
    : manifolds(nullptr), numoftypes(0), numoftotalmani(0) {
    name = "ProductManifold";
    empty = ManifoldVector(0, 0, false);
}

ProductManifold::ProductManifold(const std::vector<const Manifold*>& base_manifolds,
                               const std::vector<int>& powers) 
    : manifolds(nullptr), numoftypes(0), numoftotalmani(0) {
    initialize(base_manifolds, powers);
}

ProductManifold::ProductManifold(const ProductManifold& other)
    : Manifold(other), manifolds(nullptr) {
    copy_from(other);
}

ProductManifold::~ProductManifold() {
    cleanup();
}

ProductManifold& ProductManifold::operator=(const ProductManifold& other) {
    if (this != &other) {
        cleanup();
        copy_from(other);
    }
    return *this;
}

void ProductManifold::initialize(const std::vector<const Manifold*>& base_manifolds,
                               const std::vector<int>& powers) {
    if (base_manifolds.empty() || powers.empty() || 
        base_manifolds.size() != powers.size()) {
        throw std::runtime_error("Invalid manifold configuration");
    }

    for (int p : powers) {
        if (p <= 0) {
            throw std::runtime_error("Powers must be positive");
        }
    }

    numoftypes = static_cast<int>(base_manifolds.size());
    manifolds = new Manifold*[numoftypes];
    
    // Calculate intervals and total manifolds
    powsinterval.resize(numoftypes + 1);
    powsinterval[0] = 0;
    numoftotalmani = 0;

    for (int i = 0; i < numoftypes; ++i) {
        if (!base_manifolds[i]) {
            throw std::runtime_error("Null manifold pointer at index " + std::to_string(i));
        }
        manifolds[i] = const_cast<Manifold*>(base_manifolds[i]);
        numoftotalmani += powers[i];
        powsinterval[i + 1] = numoftotalmani;
    }

    // Set up empty vector
    size_t total_rows = 0;
    // allow mixed types of manifolds
    // 
    bool is_complex = false;
    size_t cols = base_manifolds[0]->empty.n_cols();

    for (int i = 0; i < numoftypes; ++i) {
        if (base_manifolds[i]->empty.n_cols() != cols) {
            throw std::runtime_error("All manifolds must have same number of columns");
        }
        total_rows += base_manifolds[i]->empty.n_rows() * powers[i];
        is_complex = is_complex || base_manifolds[i]->empty.is_complex();
    }

    empty = ManifoldVector(total_rows, cols, is_complex);

    // Set name
    std::ostringstream oss;
    oss << "ProductManifold(";
    for (int i = 0; i < numoftypes; ++i) {
        oss << manifolds[i]->name;
        if (powers[i] > 1) {
            oss << "^" << powers[i];
        }
        if (i < numoftypes - 1) {
            oss << " × ";
        }
    }
    oss << ")";
    name = oss.str();
}

void ProductManifold::cleanup() {
    delete[] manifolds;
    manifolds = nullptr;
    numoftypes = 0;
    powsinterval.clear();
    numoftotalmani = 0;
}

void ProductManifold::copy_from(const ProductManifold& other) {
    numoftypes = other.numoftypes;
    numoftotalmani = other.numoftotalmani;
    powsinterval = other.powsinterval;
    name = other.name;
    empty = other.empty;

    manifolds = new Manifold*[numoftypes];
    for (int i = 0; i < numoftypes; ++i) {
        manifolds[i] = other.manifolds[i];
    }
}

// Change from auto to explicit pair
std::pair<int, int> ProductManifold::get_manifold_dimensions(int type_idx) const {
    if (type_idx < 0 || type_idx >= numoftypes) {
        throw std::out_of_range("Invalid manifold type index");
    }
    return std::make_pair(manifolds[type_idx]->empty.n_rows(),
                         manifolds[type_idx]->empty.n_cols());
}

ManifoldPoint ProductManifold::extract_submanifold(const ManifoldPoint& x,
                                                 int type_idx,
                                                 int copy_idx) const {
    // Explicitly declare the pair components
    std::pair<int, int> dims = get_manifold_dimensions(type_idx);
    int rows = dims.first;
    int cols = dims.second;
    int start_row = 0;
    
    // Calculate starting row for this submanifold
    for (int i = 0; i < type_idx; ++i) {
        std::pair<int, int> mani_dims = get_manifold_dimensions(i);
        start_row += mani_dims.first * (powsinterval[i + 1] - powsinterval[i]);
    }
    start_row += rows * copy_idx;
    
    return x.submat(start_row, 0, start_row + rows - 1, cols - 1);
}

void ProductManifold::check_dimensions(const ManifoldPoint& x,
                                     const std::string& name) const {
    if (x.n_rows() != empty.n_rows() || x.n_cols() != empty.n_cols()) {
        throw std::runtime_error(name + " has wrong dimensions. Expected " + 
                               std::to_string(empty.n_rows()) + "x" + 
                               std::to_string(empty.n_cols()) + ", got " +
                               std::to_string(x.n_rows()) + "x" + 
                               std::to_string(x.n_cols()));
    }
}

double ProductManifold::metric(const ManifoldPoint& x,
                             const ManifoldVector& etax,
                             const ManifoldVector& xix) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");
    check_dimensions(xix, "xix");

    double sum = 0.0;
    for (int i = 0; i < numoftypes; ++i) {
        for (int j = powsinterval[i]; j < powsinterval[i + 1]; ++j) {
            ManifoldPoint x_sub = extract_submanifold(x, i, j - powsinterval[i]);
            ManifoldVector etax_sub = extract_submanifold(etax, i, j - powsinterval[i]);
            ManifoldVector xix_sub = extract_submanifold(xix, i, j - powsinterval[i]);
            sum += manifolds[i]->metric(x_sub, etax_sub, xix_sub);
        }
    }
    return sum;
}

// In other methods, replace auto with explicit types
ManifoldVector ProductManifold::projection(const ManifoldPoint& x,
                                         const ManifoldVector& etax) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");

    ManifoldVector result(empty.n_rows(), empty.n_cols(), empty.is_complex());
    int current_row = 0;

    for (int i = 0; i < numoftypes; ++i) {
        for (int j = powsinterval[i]; j < powsinterval[i + 1]; ++j) {
            ManifoldPoint x_sub = extract_submanifold(x, i, j - powsinterval[i]);
            ManifoldVector etax_sub = extract_submanifold(etax, i, j - powsinterval[i]);
            
            ManifoldVector proj = manifolds[i]->projection(x_sub, etax_sub);
            std::pair<int, int> dims = get_manifold_dimensions(i);
            result.submat(current_row, 0, current_row + dims.first - 1, dims.second - 1) = proj;
            current_row += dims.first;
        }
    }
    return result;
}

ManifoldPoint ProductManifold::retraction(const ManifoldPoint& x,
                                        const ManifoldVector& etax) const {
    check_dimensions(x, "x");
    check_dimensions(etax, "etax");

    ManifoldPoint result(empty.n_rows(), empty.n_cols(), empty.is_complex());
    int current_row = 0;

    for (int i = 0; i < numoftypes; ++i) {
        for (int j = powsinterval[i]; j < powsinterval[i + 1]; ++j) {
            ManifoldPoint x_sub = extract_submanifold(x, i, j - powsinterval[i]);
            ManifoldVector etax_sub = extract_submanifold(etax, i, j - powsinterval[i]);
            
            ManifoldPoint retr = manifolds[i]->retraction(x_sub, etax_sub);
            std::pair<int, int> dims = get_manifold_dimensions(i);
            result.submat(current_row, 0, current_row + dims.first - 1, dims.second - 1) = retr;
            current_row += dims.first;
        }
    }
    return result;
}

ManifoldVector ProductManifold::vector_transport(const ManifoldPoint& x,
                                               const ManifoldVector& etax,
                                               const ManifoldPoint& y,
                                               const ManifoldVector& xix) const {
    check_dimensions(x, "x");
    check_dimensions(y, "y");
    check_dimensions(etax, "etax");
    check_dimensions(xix, "xix");

    ManifoldVector result(empty.n_rows(), empty.n_cols(), empty.is_complex());
    int current_row = 0;

    for (int i = 0; i < numoftypes; ++i) {
        for (int j = powsinterval[i]; j < powsinterval[i + 1]; ++j) {
            ManifoldPoint x_sub = extract_submanifold(x, i, j - powsinterval[i]);
            ManifoldPoint y_sub = extract_submanifold(y, i, j - powsinterval[i]);
            ManifoldVector etax_sub = extract_submanifold(etax, i, j - powsinterval[i]);
            ManifoldVector xix_sub = extract_submanifold(xix, i, j - powsinterval[i]);
            
            ManifoldVector trans = manifolds[i]->vector_transport(x_sub, etax_sub,
                                                                y_sub, xix_sub);
            std::pair<int, int> dims = get_manifold_dimensions(i);
            result.submat(current_row, 0, current_row + dims.first - 1, dims.second - 1) = trans;
            current_row += dims.first;
        }
    }
    return result;
}

int ProductManifold::dimension() const {
    int dim = 0;
    for (int i = 0; i < numoftypes; ++i) {
        dim += manifolds[i]->dimension() * 
               (powsinterval[i + 1] - powsinterval[i]);
    }
    return dim;
}

int ProductManifold::intrinsic_dimension() const {
    int dim = 0;
    for (int i = 0; i < numoftypes; ++i) {
        dim += manifolds[i]->intrinsic_dimension() * 
               (powsinterval[i + 1] - powsinterval[i]);
    }
    return dim;
}

} // namespace OptimLight