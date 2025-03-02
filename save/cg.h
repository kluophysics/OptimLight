#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <limits>
#include <algorithm>

class ConjugateGradient {
public:
    using Vector = std::vector<double>;
    using Matrix = std::vector<std::vector<double>>;
    using Function = std::function<double(const Vector&)>;
    using Gradient = std::function<Vector(const Vector&)>;
    using Preconditioner = std::function<void(const Vector&, Vector&)>;
    using LineSearch = std::function<double(const Vector&, const Vector&, double)>;

    ConjugateGradient(double eta = 0.4, 
                      Preconditioner precondprep = [](const Vector&, Vector&) {},
                      LineSearch linesearch = defaultLineSearch)
        : eta(eta), precondprep(precondprep), linesearch(linesearch) {}

    struct State {
        Vector x;
        Vector x_previous;
        Vector g_previous;
        double f_x_previous;
        Vector y;
        Vector py;
        Vector pg;
        Vector s;
        double alpha;
    };

    void reset(State& state, const Function& obj, const Vector& x) {
        state.x = x;
        precondprep(x, state.pg);
        std::transform(state.pg.begin(), state.pg.end(), state.pg.begin(), std::negate<double>());
        state.f_x_previous = std::numeric_limits<double>::quiet_NaN();
    }

    State initialState(const Function& obj, const Gradient& grad, const Vector& initial_x) {
        State state;
        state.x = initial_x;
        state.x_previous = Vector(initial_x.size(), 0.0);
        state.g_previous = Vector(initial_x.size(), 0.0);
        state.f_x_previous = std::numeric_limits<double>::quiet_NaN();
        state.y = Vector(initial_x.size(), 0.0);
        state.py = Vector(initial_x.size(), 0.0);
        state.pg = grad(initial_x);
        state.s = state.pg;
        std::transform(state.s.begin(), state.s.end(), state.s.begin(), std::negate<double>());
        return state;
    }

    void updateState(State& state, const Function& obj, const Gradient& grad) {
        state.g_previous = grad(state.x);
        state.alpha = linesearch(state.x, state.s, 1.0);
        for (size_t i = 0; i < state.x.size(); ++i) {
            state.x[i] += state.alpha * state.s[i];
        }
        state.y = grad(state.x);
        for (size_t i = 0; i < state.y.size(); ++i) {
            state.y[i] -= state.g_previous[i];
        }
        precondprep(state.x, state.pg);
        double dPd = dot(state.s, state.pg);
        double etak = eta * dot(state.s, state.g_previous) / dPd;
        state.py = state.pg;
        precondprep(state.x, state.pg);
        for (size_t i = 0; i < state.pg.size(); ++i) {
            state.py[i] = state.pg[i] - state.py[i];
        }
        double ydots = dot(state.y, state.s);
        double betak = (dot(state.y, state.pg) - dot(state.y, state.py) * dot(grad(state.x), state.s) / ydots) / ydots;
        double beta = std::max(betak, etak);
        for (size_t i = 0; i < state.s.size(); ++i) {
            state.s[i] = beta * state.s[i] - state.pg[i];
        }
    }

    static double defaultLineSearch(const Vector& x, const Vector& s, double alpha) {
        // Implement a simple line search algorithm
        return alpha;
    }

private:
    double eta;
    Preconditioner precondprep;
    LineSearch linesearch;

    double dot(const Vector& a, const Vector& b) const {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }
};