#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdexcept>

// Simpson's Rule
double simpsons_rule(double (*f)(double), double a, double b, int n) {
    if (n % 2 != 0) {
        throw std::invalid_argument("n must be even for Simpson's Rule");
    }

    double h = (b - a) / n;
    double result = f(a) + f(b);

    for (int i = 1; i < n; i += 2) {
        result += 4 * f(a + i * h);
    }
    for (int i = 2; i < n; i += 2) {
        result += 2 * f(a + i * h);
    }

    return result * h / 3.0;
}

// Part (a): x^{-1/4} * sin(x), using substitution x = t^4 => dx = 4t^3 dt
double f_a_substituted(double t) {
    double x = std::pow(t, 4);
    return 4 * std::pow(t, 3) * std::pow(x, -0.25) * std::sin(x);  // dx/dt * f(x)
}

// Part (b): x^{-4} * sin(x), substitution x = 1/t => dx = -1/t^2 dt
double f_b_substituted(double t) {
    double x = 1.0 / t;
    return std::sin(x);  // since dx = -1/t^2, it cancels with x^4
}

int main() {
    std::cout << std::fixed << std::setprecision(10);

    // Part (a)
    double a1 = std::pow(0.01, 0.25);  // lower bound for t
    double b1 = 1.0;
    int n1 = 4;
    double result_a = simpsons_rule(f_a_substituted, a1, b1, n1);
    std::cout << "Approximation of integral (a): " << result_a << std::endl;

    // Part (b)
    double a2 = 0.01;
    double b2 = 1.0;
    int n2 = 4;
    double result_b = simpsons_rule(f_b_substituted, a2, b2, n2);
    std::cout << "Approximation of integral (b): " << result_b << std::endl;

    return 0;
}
