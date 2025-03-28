#include "cg.hpp"
#include <iostream>

int main() {
    using Vector = ConjugateGradient::Vector;

    // Define the objective function (e.g., a simple quadratic function)
    auto objective = [](const Vector& x) {
        return 0.5 * (x[0] * x[0] + x[1] * x[1]);
    };

    // Define the gradient of the objective function
    auto gradient = [](const Vector& x) {
        return Vector{x[0], x[1]};
    };

    // Initial guess
    Vector initial_x = {5.0, 5.0};

    // Create an instance of ConjugateGradient
    ConjugateGradient cg;

    // Initialize the state
    auto state = cg.initialState(objective, gradient, initial_x);

    // Perform a few iterations
    for (int i = 0; i < 20; ++i) {
        cg.updateState(state, objective, gradient);
        std::cout << "Iteration " << i << ": x = [" << state.x[0] << ", " << state.x[1] << "]\n";
    }

    return 0;
}