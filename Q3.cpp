#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <stdexcept>

// 被積分函數 f(x, y)
double f(double x, double y) {
    return 2 * y * std::sin(x) + std::cos(x) * std::cos(x);
}

// y 的上下限是 x 的函數
double y_lower(double x) {
    return std::sin(x);
}
double y_upper(double x) {
    return std::cos(x);
}

// Composite Simpson's Rule (雙重積分)
double simpsons_double(double a, double b, int n, int m) {
    std::vector<double> x(n + 1);
    double hx = (b - a) / n;
    for (int i = 0; i <= n; ++i)
        x[i] = a + i * hx;

    std::vector<double> outer_vals(n + 1);

    for (int i = 0; i <= n; ++i) {
        double xi = x[i];
        double y0 = y_lower(xi);
        double y1 = y_upper(xi);
        double hy = (y1 - y0) / m;

        std::vector<double> y(m + 1);
        for (int j = 0; j <= m; ++j)
            y[j] = y0 + j * hy;

        double inner = f(xi, y[0]) + f(xi, y[m]);
        for (int j = 1; j < m; j += 2)
            inner += 4 * f(xi, y[j]);
        for (int j = 2; j < m - 1; j += 2)
            inner += 2 * f(xi, y[j]);
        inner *= hy / 3.0;

        outer_vals[i] = inner;
    }

    double result = outer_vals[0] + outer_vals[n];
    for (int i = 1; i < n; i += 2)
        result += 4 * outer_vals[i];
    for (int i = 2; i < n - 1; i += 2)
        result += 2 * outer_vals[i];
    result *= hx / 3.0;

    return result;
}

// Gaussian Quadrature (n=3 for x and y)
double gaussian_double(double a, double b, int n) {
    std::vector<double> nodes, weights;

    if (n == 3) {
        nodes = { -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0) };
        weights = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    } else {
        throw std::invalid_argument("Only n=3 supported.");
    }

    double total = 0.0;
    for (int i = 0; i < n; ++i) {
        double xi = 0.5 * (b - a) * nodes[i] + 0.5 * (b + a);
        double wi = weights[i];

        double y0 = y_lower(xi);
        double y1 = y_upper(xi);

        double inner = 0.0;
        for (int j = 0; j < n; ++j) {
            double yj = 0.5 * (y1 - y0) * nodes[j] + 0.5 * (y1 + y0);
            double wj = weights[j];
            inner += wj * f(xi, yj);
        }

        total += wi * inner * 0.5 * (y1 - y0);
    }

    return 0.5 * (b - a) * total;
}

// 使用高精度 Simpson's Rule 計算 Exact Value
double exact_value_simpson(double a, double b, int nx, int ny) {
    return simpsons_double(a, b, nx, ny);
}

int main() {
    double a = 0.0;
    double b = M_PI / 4.0;

    // Exact value 用更細的格子
    double exact = exact_value_simpson(a, b, 100, 100);
    double simpson = simpsons_double(a, b, 4, 4);
    double gaussian = gaussian_double(a, b, 3);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Exact value:           " << exact << std::endl;
    std::cout << "Simpson's Rule:        " << simpson << std::endl;
    std::cout << "Gaussian Quadrature:   " << gaussian << std::endl;
    std::cout << "Error (Simpson):       " << std::abs(simpson - exact) << std::endl;
    std::cout << "Error (Gaussian):      " << std::abs(gaussian - exact) << std::endl;

    return 0;
}

