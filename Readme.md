# OptimLight
It is a collaborative project targeting in optimization using C++ language. It is intended to support
- Unconstrained optimization
- Box constraints
- Manifold optimzation (not yet)

# OptimLight

OptimLight is a C++ framework designed for performing optimization tasks, including line search and trust-region optimization methods. It supports various optimization algorithms such as steepest descent, conjugate gradient, BFGS, Newton method, and dogleg. The framework is flexible, allowing users to work with both double and complex types, and it can utilize popular linear algebra libraries like Armadillo or Eigen, or even user-defined implementations.

## Features

- **Optimization Algorithms**: Implements various optimization methods including:
  - Steepest Descent
  - Conjugate Gradient
  - BFGS
  - Newton Method
  - Dogleg

- **Type Support**: Supports both double and complex types for optimization problems.

- **Manifold Support**: Utilizes a stack of manifolds for advanced optimization techniques.

- **Library Integration**: Options to use default libraries like Armadillo or Eigen, or custom implementations.

- **Testing Framework**: Includes setups for Google Unit Tests to ensure code reliability and correctness.

## Installation

To build the OptimLight framework, you will need CMake and a compatible C++ compiler. Follow these steps:

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd optimlight
   ```

2. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```

3. Configure the project using CMake:
   ```bash
   cmake ..
   ```

4. Build the project:
   ```bash
   make
   ```

## Usage

To use the OptimLight framework, include the necessary headers in your C++ files. Here is a simple example of how to set up an optimization problem:

```cpp
#include "core/problem.hpp"
#include "optimizers/line_search/steepest_descent.hpp"

// Define your optimization problem and use the desired optimizer
```

For more detailed usage examples, please refer to the `examples` directory.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue for any suggestions or improvements.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Need your feedback and help!