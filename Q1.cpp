#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

// 被積函數
double f(double x) {
    return std::sin(4 * x) * std::exp(x);
}

// 建立等距向量（模擬 np.linspace）
std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> vec(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        vec[i] = start + i * step;
    }
    return vec;
}

// Composite Trapezoidal Rule
double trapezoidal(double a, double b, double h, int n) {
    std::vector<double> x = linspace(a, b, n + 1);
    double sum = 0.0;
    for (int i = 1; i < n; ++i) {
        sum += f(x[i]);
    }
    double T = h * (0.5 * f(x[0]) + sum + 0.5 * f(x[n]));
    return T;
}

// Composite Simpson's Rule
double simpson(double a, double b, double h, int n) {
    if (n % 2 != 0) {
        throw std::invalid_argument("n must be even for Simpson's Rule");
    }
    std::vector<double> x = linspace(a, b, n + 1);
    double sum_even = 0.0;
    double sum_odd = 0.0;
    for (int i = 2; i < n; i += 2) {
        sum_even += f(x[i]);
    }
    for (int i = 1; i < n; i += 2) {
        sum_odd += f(x[i]);
    }
    double S = h / 3 * (f(x[0]) + 2 * sum_even + 4 * sum_odd + f(x[n]));
    return S;
}

// Composite Midpoint Rule
double midpoint(double a, double b, double h, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double midpoint = a + h * i + h / 2.0;
        sum += f(midpoint);
    }
    return h * sum;
}

int main() {
    double a = 1.0;
    double b = 2.0;
    double h = 0.1;
    int n = static_cast<int>((b - a) / h);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Composite Trapezoidal Rule: " << trapezoidal(a, b, h, n) << std::endl;
    std::cout << "Composite Simpson's Rule:   " << simpson(a, b, h, n) << std::endl;
    std::cout << "Composite Midpoint Rule:    " << midpoint(a, b, h, n) << std::endl;

    return 0;
}
